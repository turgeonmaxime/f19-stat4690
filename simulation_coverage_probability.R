library(tidyverse)
library(mvtnorm)

set.seed(123)

B <- 10000
alpha <- 0.05
p <- 2; n <- 100
mu <- c(0, 0)
Sigma <- matrix(c(1, 0.5, 0.5, 1), ncol = p)

results <- map_df(seq_len(B), function(i) {
  # Generate dataset
  data <- rmvnorm(n, mean = mu, sigma = Sigma)
  # Compute test statistics
  mu_hat <- colMeans(data)
  Sn <- cov(data)
  
  # Four confidence regions
  # 1. Univariate without correction
  univ_ci <- cbind(mu_hat - qt(1-0.5*alpha, n - 1)*sqrt(diag(Sn)/n),
                   mu_hat + qt(1-0.5*alpha, n - 1)*sqrt(diag(Sn)/n))
  # 2. Bonferroni adjustment
  bonf_ci <- cbind(mu_hat - qt(1-0.5*alpha/p, n - 1)*sqrt(diag(Sn)/n),
                   mu_hat + qt(1-0.5*alpha/p, n - 1)*sqrt(diag(Sn)/n))
  
  # 3. Simultaneous confidence intervals
  simul_ci <- cbind(mu_hat - ((n - 1)*p/(n-p))*qf(1-0.5*alpha, df1 = p,
                                                  df2 = n - p)*sqrt(diag(Sn)/n),
                    mu_hat + ((n - 1)*p/(n-p))*qf(1-0.5*alpha, df1 = p,
                                          df2 = n - p)*sqrt(diag(Sn)/n))
  # 4. Confidence region
  test_statistic <- n * t(mu_hat - mu) %*% 
    solve(Sn) %*% (mu_hat - mu)
  
  critical_val <- (n - 1)*p*qf(1-0.5*alpha, df1 = p,
                               df2 = n - p)/(n-p)
  
  cr_res <- (drop(test_statistic) > critical_val)
  
  # Did we reject the null hypothesis?
  data.frame(univ = !all(univ_ci[,1] < mu & univ_ci[,2] > mu),
             bonf = !all(bonf_ci[,1] < mu & bonf_ci[,2] > mu),
             simul = !all(simul_ci[,1] < mu & simul_ci[,2] > mu),
             ellipse = cr_res)
})

results %>% summarise_all(mean)
