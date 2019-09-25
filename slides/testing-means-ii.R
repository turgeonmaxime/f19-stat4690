## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ----message = FALSE-----------------------------------------------------
library(tidyverse)
library(dslabs)

dataset <- gapminder %>% 
  filter(year %in% 2012:2016, 
         continent == "Africa") %>% 
  select(year, country, life_expectancy)


## ----message = FALSE-----------------------------------------------------
# QQ-plots to assess normality
dataset %>% 
  ggplot(aes(sample = life_expectancy)) + 
  stat_qq() + stat_qq_line() + 
  facet_wrap(~year)


## ------------------------------------------------------------------------
C <- matrix(c(1, -1, 0, 0, 0,
              1, 0, -1, 0, 0,
              1, 0, 0, -1, 0,
              1, 0, 0, 0, -1),
            ncol = 5, byrow = TRUE)
C


## ------------------------------------------------------------------------
# Transform data into wide format
dataset <- dataset %>% 
  spread(year, life_expectancy)

head(dataset)


## ------------------------------------------------------------------------
# Compute test statistic
dataset <- dataset %>% 
  select(-country) %>% 
  as.matrix()
n <- nrow(dataset); p <- ncol(dataset)

mu_hat <- colMeans(dataset)
mu_hat


## ----eval = TRUE---------------------------------------------------------
Sn <- cov(dataset)
test_statistic <- n * t(C %*% mu_hat) %*% 
  solve(C %*% Sn %*% t(C)) %*% (C %*% mu_hat)

const <- (n - 1)*(p - 1)/(n - p + 1)
critical_val <- const * qf(0.95, df1 = p - 1,
                           df2 = n - p + 1)

drop(test_statistic) > critical_val


## ------------------------------------------------------------------------
alpha <- 0.05
mu_contr <- C %*% mu_hat
sample_cov <- diag(C %*% Sn %*% t(C))

mu_contr


## ------------------------------------------------------------------------
# Simultaneous CIs
simul_ci <- cbind(mu_contr - sqrt(critical_val*
                                    sample_cov/n),
                  mu_contr + sqrt(critical_val*
                                    sample_cov/n))


## ------------------------------------------------------------------------
# Bonferroni adjustment
bonf_ci <- cbind(mu_contr - qt(1-0.5*alpha/(p-1), 
                               n - 1) *
                   sqrt(sample_cov/n),
                 mu_contr + qt(1-0.5*alpha/(p-1), 
                               n - 1) *
                   sqrt(sample_cov/n))


## ------------------------------------------------------------------------
simul_ci
bonf_ci


## ----message = FALSE-----------------------------------------------------
dataset1 <- gapminder %>% 
  filter(year == 2012, 
         continent == "Africa",
         !is.na(infant_mortality)) %>% 
  select(life_expectancy, infant_mortality) %>% 
  as.matrix()
dim(dataset1)

dataset2 <- gapminder %>% 
  filter(year == 2012, 
         continent == "Asia",
         !is.na(infant_mortality)) %>% 
  select(life_expectancy, infant_mortality) %>% 
  as.matrix()
dim(dataset2)

n1 <- nrow(dataset1); n2 <- nrow(dataset2)
p <- ncol(dataset1)


## ------------------------------------------------------------------------
(mu_hat1 <- colMeans(dataset1))
(mu_hat2 <- colMeans(dataset2))

(S1 <- cov(dataset1))
(S2 <- cov(dataset2))

# Even though it doesn't look reasonable
# We will assume equal covariance for now


## ------------------------------------------------------------------------
mu_hat_diff <- mu_hat1 - mu_hat2

S_pool <- ((n1 - 1)*S1 + (n2 - 1)*S2)/(n1+n2-2)

test_statistic <- t(mu_hat_diff) %*% 
  solve((n1^-1 + n2^-1)*S_pool) %*% mu_hat_diff

const <- (n1 + n2 - 2)*p/(n1 + n2 - p - 2)
critical_val <- const * qf(0.95, df1 = p,
                           df2 = n1 + n2 - p - 2)

drop(test_statistic) > critical_val


## ----echo = FALSE--------------------------------------------------------
R <- solve((n1^-1 + n2^-1)*S_pool)
transf <- chol(R)

# First create a circle of radius c
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- sqrt(critical_val) * cbind(cos(theta_vect), sin(theta_vect))
# Then turn into ellipse
ellipse1 <- circle %*% t(solve(transf)) + 
  matrix(mu_hat_diff, ncol = p, 
         nrow = nrow(circle), 
         byrow = TRUE)

plot(ellipse1, type = 'l',
     xlab = colnames(dataset1)[1],
     ylab = colnames(dataset1)[2],
     main = "Comparing Africa vs. Asia")
points(x = mu_hat_diff[1],
       y = mu_hat_diff[2])


## ------------------------------------------------------------------------
test_statistic <- t(mu_hat_diff) %*% 
  solve(n1^-1*S1 + n2^-1*S2) %*% mu_hat_diff

critical_val <- qchisq(0.95, df = p)

drop(test_statistic) > critical_val


## ----echo = FALSE--------------------------------------------------------
R <- solve(n1^-1*S1 + n2^-1*S2)
transf <- chol(R)

# First create a circle of radius c
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- sqrt(critical_val) * cbind(cos(theta_vect), sin(theta_vect))
# Then turn into ellipse
ellipse2 <- circle %*% t(solve(transf)) + 
  matrix(mu_hat_diff, ncol = p, 
         nrow = nrow(circle), 
         byrow = TRUE)


## ------------------------------------------------------------------------
W1 <- S1 %*% solve(n1^-1*S1 + n2^-1*S2)/n1
W2 <- S2 %*% solve(n1^-1*S1 + n2^-1*S2)/n2

trace_square <- sum(diag(W1%*%W1))/n1 + 
  sum(diag(W2%*%W2))/n2
square_trace <- sum(diag(W1))^2/n1 + 
  sum(diag(W2))^2/n2

(nu <- (p + p^2)/(trace_square + square_trace))


## ------------------------------------------------------------------------
const <- nu*p/(nu - p - 1)
critical_val <- const * qf(0.95, df1 = p,
                           df2 = nu - p - 1)

drop(test_statistic) > critical_val


## ----echo = FALSE--------------------------------------------------------
# First create a circle of radius c
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- sqrt(critical_val) * cbind(cos(theta_vect), sin(theta_vect))
# Then turn into ellipse
ellipse3 <- circle %*% t(solve(transf)) + 
  matrix(mu_hat_diff, ncol = p, 
         nrow = nrow(circle), 
         byrow = TRUE)

xlim <- range(c(ellipse1[,1], 
                ellipse2[,1],
                ellipse3[,1]))
ylim <- range(c(ellipse1[,2], 
                ellipse2[,2],
                ellipse3[,2]))

plot(ellipse2, type = 'l',
     xlab = colnames(dataset1)[1],
     ylab = colnames(dataset1)[2],
     xlim = xlim, ylim = ylim,
     main = "Comparing Africa vs. Asia")
points(x = mu_hat_diff[1],
       y = mu_hat_diff[2])
lines(ellipse1, lty = 2)
lines(ellipse3, lty = 3)
legend('topright', legend = c("Unequal", "Equal", "Nel-VDM"), lty = 1:3)

