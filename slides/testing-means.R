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

