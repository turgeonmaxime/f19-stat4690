## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ------------------------------------------------------------------------
library(mvtnorm)
set.seed(123)

n <- 50; p <- 2

mu <- c(1, 2)
Sigma <- matrix(c(1, 0.5, 0.5, 1), ncol = p)

Y <- rmvnorm(n, mean = mu, sigma = Sigma)


## ------------------------------------------------------------------------
loglik <- function(mu, sigma, data = Y) {
  # Compute quantities
  y_bar <- colMeans(Y)
  Sn <- (n-1)*cov(Y)/n
  Sigma_inv <- solve(sigma)
  
  # Compute quadratic form
  quad_form <- drop(t(y_bar - mu) %*% Sigma_inv %*% 
                      (y_bar - mu))
  
  -0.5*n*(log(det(sigma)) + 
            sum(diag(Sigma_inv %*% Sn)) +
            quad_form)
}


## ------------------------------------------------------------------------
grid_xy <- expand.grid(seq(0.5, 1.5, 
                           length.out = 32), 
                       seq(1, 3, 
                           length.out = 32))

contours <- purrr::map_df(seq_len(nrow(grid_xy)), 
                          function(i) {
  # Where we will evaluate loglik
  mu_obs <- as.numeric(grid_xy[i,])
  # Evaluate at the pop covariance
  z <- loglik(mu_obs, sigma = Sigma)
  # Output data.frame
  data.frame(x = mu_obs[1],
             y = mu_obs[2],
             z = z)
})


## ---- message = FALSE----------------------------------------------------
library(tidyverse)
library(ggrepel)
# Create df with pop and sample means
data_means <- data.frame(x = c(mu[1], mean(Y[,1])),
                         y = c(mu[2], mean(Y[,2])),
                         label = c("Pop.", "Sample"))


## ---- message = FALSE----------------------------------------------------
contours %>% 
  ggplot(aes(x, y)) + 
  geom_contour(aes(z = z)) + 
  geom_point(data = data_means) +
  geom_label_repel(data = data_means,
                   aes(label = label))


## ------------------------------------------------------------------------
library(scatterplot3d)
with(contours, scatterplot3d(x, y, z))


## ---- echo = FALSE, eval = FALSE-----------------------------------------
## x <- seq(0.5, 1.5, length.out = 32)
## y <- seq(1, 3, length.out = 32)
## z <- matrix(NA, ncol = 32, nrow = 32)
## for (i in seq_len(32)) {
##   for (j in seq_len(32))
##     z[i,j] <- loglik(c(x[i], y[j]), sigma = Sigma)
## }
## persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")

