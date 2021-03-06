## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ---- message=FALSE------------------------------------------------------
library(dslabs)
library(tidyverse)

dataset <- gapminder %>% 
  filter(year == 2012, 
         !is.na(infant_mortality)) %>% 
  select(infant_mortality, 
         life_expectancy, 
         fertility) %>% 
  as.matrix()


## ---- message=FALSE------------------------------------------------------
# Assume we know Sigma
Sigma <- matrix(c(555, -170, 30, -170, 65, -10, 
                  30, -10, 2), ncol = 3)

mu_hat <- colMeans(dataset) 
mu_hat


## ---- message=FALSE------------------------------------------------------
# Test mu = mu_0
mu_0 <- c(25, 50, 3)
test_statistic <- nrow(dataset) * t(mu_hat - mu_0) %*% 
  solve(Sigma) %*% (mu_hat - mu_0)

drop(test_statistic) > qchisq(0.95, df = 3)


## ---- message=FALSE------------------------------------------------------
n <- nrow(dataset); p <- ncol(dataset)

# Test mu = mu_0
mu_0 <- c(25, 50, 3)
test_statistic <- n * t(mu_hat - mu_0) %*% 
  solve(cov(dataset)) %*% (mu_hat - mu_0)

critical_val <- (n - 1)*p*qf(0.95, df1 = p,
                             df2 = n - p)/(n-p)

drop(test_statistic) > critical_val


## ------------------------------------------------------------------------
n <- nrow(dataset); p <- ncol(dataset)

# Test mu = mu_0
mu_0 <- c(25, 50, 3)
test_statistic <- n * t(mu_hat - mu_0) %*% 
  solve(cov(dataset)) %*% (mu_hat - mu_0)

critical_val <- (n - 1)*p*qf(0.95, df1 = p,
                             df2 = n - p)/(n-p)
sample_cov <- diag(cov(dataset))

cbind(mu_hat - sqrt(critical_val*
                      sample_cov/n),
      mu_hat + sqrt(critical_val*
                      sample_cov/n))


## ------------------------------------------------------------------------
U <- matrix(c(1, 0, 0,
              0, 1, 0), 
            ncol = 2)
R <- n*solve(t(U) %*% cov(dataset) %*% U)
transf <- chol(R)


## ------------------------------------------------------------------------
# First create a circle of radius c
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- sqrt(critical_val) * cbind(cos(theta_vect), sin(theta_vect))
# Then turn into ellipse
ellipse <- circle %*% t(solve(transf)) + 
  matrix(mu_hat[1:2], ncol = 2, 
         nrow = nrow(circle), 
         byrow = TRUE)


## ------------------------------------------------------------------------
# Eigendecomposition
decomp <- eigen(t(U) %*% cov(dataset) %*% U)
first <- sqrt(decomp$values[1]) *
  decomp$vectors[,1] * sqrt(critical_val)
second <- sqrt(decomp$values[2]) * 
  decomp$vectors[,2] * sqrt(critical_val)


## ---- echo = FALSE-------------------------------------------------------
plot(ellipse, type = 'l',
     xlab = colnames(dataset)[1],
     ylab = colnames(dataset)[2])
lines(x = c(0 + mu_hat[1], first[1]/sqrt(n) + mu_hat[1]),
      y = c(0 + mu_hat[2], first[2]/sqrt(n) + mu_hat[2]))
lines(x = c(0 + mu_hat[1], -first[1]/sqrt(n) + mu_hat[1]),
      y = c(0 + mu_hat[2], -first[2]/sqrt(n) + mu_hat[2]))
lines(x = c(0 + mu_hat[1], second[1]/sqrt(n) + mu_hat[1]),
      y = c(0 + mu_hat[2], second[2]/sqrt(n) + mu_hat[2]))
lines(x = c(0 + mu_hat[1], -second[1]/sqrt(n) + mu_hat[1]),
      y = c(0 + mu_hat[2], -second[2]/sqrt(n) + mu_hat[2]))
points(x = mu_hat[1],
       y = mu_hat[2])

# Add 1d projections
axis_proj <- cbind(mu_hat - sqrt(critical_val*
                                   sample_cov/n),
                   mu_hat + sqrt(critical_val*
                                   sample_cov/n))
abline(v = axis_proj[1,], lty = 2)
abline(h = axis_proj[2,], lty = 2)


## ---- message=FALSE------------------------------------------------------
# Let's focus on only two variables
dataset <- gapminder %>% 
  filter(year == 2012, 
         !is.na(infant_mortality)) %>% 
  select(infant_mortality, 
         life_expectancy) %>% 
  as.matrix()

n <- nrow(dataset); p <- ncol(dataset)


## ------------------------------------------------------------------------
alpha <- 0.05
mu_hat <- colMeans(dataset) 
sample_cov <- diag(cov(dataset))

# Simultaneous CIs
critical_val <- (n - 1)*p*qf(1-0.5*alpha, df1 = p,
                             df2 = n - p)/(n-p)

simul_ci <- cbind(mu_hat - sqrt(critical_val*
                                  sample_cov/n),
                  mu_hat + sqrt(critical_val*
                                  sample_cov/n))


## ------------------------------------------------------------------------
# Univariate without correction
univ_ci <- cbind(mu_hat - qt(1-0.5*alpha, n - 1) *
                   sqrt(sample_cov/n),
                 mu_hat + qt(1-0.5*alpha, n - 1) *
                   sqrt(sample_cov/n))


## ------------------------------------------------------------------------
# Bonferroni adjustment
bonf_ci <- cbind(mu_hat - qt(1-0.5*alpha/p, n - 1) *
                   sqrt(sample_cov/n),
                 mu_hat + qt(1-0.5*alpha/p, n - 1) *
                   sqrt(sample_cov/n))


## ------------------------------------------------------------------------
simul_ci
univ_ci
bonf_ci


## ---- echo = FALSE-------------------------------------------------------
transf_mat <- solve(chol(solve(cov(dataset)/n)))
# First create a circle of radius c
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- sqrt(critical_val) * cbind(cos(theta_vect), sin(theta_vect))
ellipse1 <- circle %*% t(transf_mat)
ellipse2 <- t(apply(ellipse1, 1, function(row) row + mu_hat))
colnames(ellipse2) <- c("X", "Y")

data_ellipse <- data.frame(ellipse2)

bind_rows(
  data.frame(t(simul_ci)) %>% mutate(type = 'T2-intervals'),
  data.frame(t(univ_ci)) %>% mutate(type = 'Unadjusted'),
  data.frame(t(bonf_ci)) %>% mutate(type = 'Bonferroni')
  ) %>% 
  ggplot() +
  geom_polygon(data = data_ellipse,
               aes(X, Y),
               fill = 'grey60')+
  geom_vline(aes(xintercept = infant_mortality,
                 linetype = type)) +
  geom_hline(aes(yintercept = life_expectancy,
                 linetype = type)) +
  theme_minimal() +
  geom_point(x = mu_hat[1],
             y = mu_hat[2],
             size = 2) +
  xlab("Infant mortality") +
  ylab("Life Expectancy") +
  theme(legend.position = "top",
        legend.title=element_blank()) +
  scale_linetype_discrete(breaks = c('T2-intervals',
                                     'Bonferroni',
                                     'Unadjusted'))

