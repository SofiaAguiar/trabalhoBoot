if (!require(pacman)) install.packages('pacman')
library(pacman)
pacman::p_load(tidyverse, data.table, EnvStats, progress, gt)

bootstrap_varbasica = function(A){
  thetaPontual = var(A)
  b = 1
  B = 200
  reamostragem = list()
  reamostragem[[1]] = c(1:B)
  reamostragem[[2]] = numeric(B)
  while (b <= B) {
    b_i = sample(A, length(A), replace = T)
    thetaBoot = var(b_i)
    reamostragem[[2]][[b]] = thetaBoot
    b = b + 1
  }
  alfa = 0.05
  thetaInf = quantile(reamostragem[[2]], alfa/2)
  thetaSup = quantile(reamostragem[[2]], 1-alfa/2)
  valorInf = 2*thetaPontual - thetaSup
  valorSup = 2*thetaPontual - thetaInf
  return(c(valorInf, valorSup))
  
}

montecarlo <- function(n, FUN, theo_val, ...) {
  
  int_confianca <- data.table(valor = rep(theo_val, 1000),
                              InfBoot = numeric(1000),
                              SupBoot = numeric(1000),
                              InfParam = numeric(1000),
                              SupParam = numeric(1000))
  pb <- progress_bar$new(total=1000)
  
  for (iter in 1:1000) {
    amostra_iter <- FUN(n, ...)
    
    int_confianca[iter, 2:5] <- data.table(
      reduce(list(bootstrap_varbasica(amostra_iter),
                  varTest(amostra_iter)$conf.int[c(1,2)]), c)) %>% transpose()
    pb$tick()
  }
  
  return(list(.001 * nrow(int_confianca[valor %between% list(InfBoot, SupBoot)]), # Bootstrap
              .001 * nrow(int_confianca[valor %between% list(InfParam, SupParam)]))) # Paramétrico
  
}

# Normal

normal_por_n_df <- tibble(x = c(10, 50, 100, 250, 500, 1000, 2000, 3000, 4000, 5000))

pb <- progress_bar$new(total = length(normal_por_n_df$x))
normal_por_n <- map(normal_por_n_df$x, ~{pb$tick(); montecarlo(.x, rnorm, 1)})

normal_por_n_df <- normal_por_n_df %>% 
  mutate(bootstrap = map_dbl(1:nrow(normal_por_n_df), ~pluck(normal_por_n, .x, 1)),
         param = map_dbl(1:nrow(normal_por_n_df), ~pluck(normal_por_n, .x, 2)))

normal_100_nvezes <- map(1:10, ~montecarlo(100, rnorm, 1))

normal_100_nvezes_df <- tibble(
  resultados_boots = map_dbl(1:10, ~pluck(normal_100_nvezes, .x, 1)),
  resultados_param = map_dbl(1:10, ~pluck(normal_100_nvezes, .x, 2))
)

resultados_gerais <- tibble(
  normal = unlist(montecarlo(1000, rnorm, 1)),
  exponencial = unlist(montecarlo(1000, rexp, .25, rate = 2)),
  binomial = unlist(montecarlo(1000, rbinom, 250, size = 1000, prob = .5))
)

# Binomial

var_binom=function(n,p){
  return(n * p * (1 - p))
}

binomial_por_n_df <- tibble(x = c(10, 50, 100, 250, 500, 1000, 2000, 3000, 4000, 5000))

pb <- progress_bar$new(total = length(binomial_por_n_df$x))
binomial_por_n <- map(binomial_por_n_df$x, ~{pb$tick(); montecarlo(.x,
                                                                   rbinom,
                                                                   theo_val=var_binom(.x, 0.5),
                                                                   prob=0.5,
                                                                   size=.x)})

binomial_por_n_df <- binomial_por_n_df %>% 
  mutate(bootstrap = map_dbl(1:nrow(binomial_por_n_df), ~pluck(binomial_por_n, .x, 1)),
         param = map_dbl(1:nrow(binomial_por_n_df), ~pluck(binomial_por_n, .x, 2)))

binomial_100_nvezes <- map(1:10, ~montecarlo(100,
                                             rbinom,
                                             theo_val=var_binom(100, 0.5),
                                             prob=0.5,
                                             size=100))

binomial_100_nvezes_df <- tibble(
  resultados_boots = map_dbl(1:10, ~pluck(binomial_100_nvezes, .x, 1)),
  resultados_param = map_dbl(1:10, ~pluck(binomial_100_nvezes, .x, 2))
)

# Exponencial

exponencial_por_n_df <- tibble(x = c(10, 50, 100, 250, 500, 1000, 2000, 3000, 4000, 5000))

pb <- progress_bar$new(total = length(exponencial_por_n_df$x))
exponencial_por_n <- map(exponencial_por_n_df$x, ~{pb$tick(); montecarlo(.x,
                                                                         rexp,
                                                                         rate=2,
                                                                         theo_val= .25)})

exponencial_por_n_df <- exponencial_por_n_df %>% 
  mutate(bootstrap = map_dbl(1:nrow(exponencial_por_n_df), ~pluck(exponencial_por_n, .x, 1)),
         param = map_dbl(1:nrow(exponencial_por_n_df), ~pluck(exponencial_por_n, .x, 2)))

exponencial_100_nvezes <- map(1:10, ~montecarlo(100, rexp, rate=2, theo_val= .25))

exponencial_100_nvezes_df <- tibble(
  resultados_boots = map_dbl(1:10, ~pluck(exponencial_100_nvezes, .x, 1)),
  resultados_param = map_dbl(1:10, ~pluck(exponencial_100_nvezes, .x, 2))
)

# Export

save(normal_por_n_df,
     normal_100_nvezes_df,
     resultados_gerais,
     binomial_por_n_df,
     binomial_100_nvezes_df,
     exponencial_por_n_df,
     exponencial_100_nvezes_df,
     file = 'simulacoes.RData')
