library(survival)
library(broom)
library(tidyverse)

n <- 100
n_it <- 1000
beta <- log(2)

hazard_invers <- function(x, lambda){
  return(1/lambda * x)
}

sim_cox <- function(it, n){
  
  set.seed(it)
  U <- runif(n)
  # x <- rnorm(n)
  x <- c(rep(0, n/2), rep(1, n/2))
  y <- hazard_invers(-1*log(U)*exp(-1*beta*x), lambda = 0.5)
  
  
  ols <- lm(y~x)
  cox_mdl <- coxph(Surv(y,rep(1,n))~x)
  data.frame(
    mean_t = mean(y),
    var_t = var(y),
    cox_coef = coef(cox_mdl)[1],
    cox_se = (vcov(cox_mdl)[1,1]) %>% sqrt,
    cox_p = cox_mdl %>% tidy %>% pull(p.value),
    cox_CI = confint(cox_mdl),
    ols_coef = coef(ols)[2],
    ols_se = vcov(ols)[2,2] %>% sqrt,
    ols_p = (ols %>% tidy %>% pull(p.value))[2]
  )
}

res_30 <- purrr::map_dfr(1:n_it, sim_cox, n = 30)
res_100 <- purrr::map_dfr(1:n_it, sim_cox, n = 100)

saveRDS(res_30, "Project_2_Simulation Study/res_30.rds")
saveRDS(res_100, "Project_2_Simulation Study/res_100.rds")
