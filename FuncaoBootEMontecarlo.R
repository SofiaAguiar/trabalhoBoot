library(tidyverse)

library(data.table)
library(EnvStats)

# 1. Função para ic bootstrap básico para uma amostra A, com 200 reamostragem e ic de 90% ----

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

#exemplo da função rodando para uma amostra normal padrao
amostra = rnorm(50)
bootstrap_varbasica(amostra)

# 2. Monte Carlo ----


montecarlo <- function(n, FUN, theo_val, ...) {
  
  int_confianca <- data.table(valor = rep(theo_val, 1000),
                              InfBoot = numeric(1000),
                              SupBoot = numeric(1000),
                              InfParam = numeric(1000),
                              SupParam = numeric(1000))
  
  for (iter in 1:1000) {
    amostra_iter <- FUN(n, ...)
    
    int_confianca[iter, 2:5] <- data.table(
      reduce(list(bootstrap_varbasica(amostra_iter),
                  varTest(amostra_iter)$conf.int[c(1,2)]), c)) %>% transpose()
     
    
  }

  return(list(.001 * nrow(int_confianca[valor %between% list(InfBoot, SupBoot)]), # Bootstrap
              .001 * nrow(int_confianca[valor %between% list(InfParam, SupParam)]))) # Paramétrico

  
}

# Exemplos:
# montecarlo(1000, rnorm, 16, mean = 15, sd = 4)
# montecarlo(1000, rexp, 0.25, rate = 2)
