## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ------------------------------------------------------------------------
Sigma_Y <- matrix(c(1, 0.4, 0.4, 1), ncol = 2)
Sigma_X <- matrix(c(1, 0.2, 0.2, 1), ncol = 2)
Sigma_YX <- matrix(c(0.5, 0.3, 0.6, 0.4), ncol = 2)
Sigma_XY <- t(Sigma_YX)

rbind(cbind(Sigma_Y, Sigma_YX),
      cbind(Sigma_XY, Sigma_X))

## ---- message = FALSE----------------------------------------------------
library(expm)
sqrt_Y <- sqrtm(Sigma_Y)
sqrt_X <- sqrtm(Sigma_X)
M1 <- solve(sqrt_Y) %*% Sigma_YX %*% solve(Sigma_X)%*% 
  Sigma_XY %*% solve(sqrt_Y)

(decomp1 <- eigen(M1))


## ------------------------------------------------------------------------
decomp1$vectors[,1] %*% solve(sqrt_Y)


## ------------------------------------------------------------------------
M2 <- solve(sqrt_X) %*% Sigma_XY %*% solve(Sigma_Y)%*% 
  Sigma_YX %*% solve(sqrt_X)

decomp2 <- eigen(M2)
decomp2$vectors[,1] %*% solve(sqrt_X)


## ------------------------------------------------------------------------
sqrt(decomp1$values)


## ------------------------------------------------------------------------
# Let's generate data
library(mvtnorm)
Sigma <- rbind(cbind(Sigma_Y, Sigma_YX),
               cbind(Sigma_XY, Sigma_X))

YX <- rmvnorm(100, sigma = Sigma)
Y <- YX[,1:2]
X <- YX[,3:4]

decomp <- cancor(x = X, y = Y)


## ------------------------------------------------------------------------
U <- Y %*% decomp$ycoef
V <- X %*% decomp$xcoef

diag(cor(U, V))
decomp$cor


## ---- message=FALSE------------------------------------------------------
library(tidyverse)
library(dslabs)

X <- olive %>% 
  select(-area, -region) %>% 
  as.matrix

Y <- olive %>% 
  select(region) %>% 
  model.matrix(~ region - 1, data = .)


## ---- message=FALSE------------------------------------------------------
head(unname(Y))

decomp <- cancor(X, Y)

V <- X %*% decomp$xcoef


## ---- message=FALSE------------------------------------------------------
data.frame(
  V1 = V[,1],
  V2 = V[,2],
  region = olive$region
) %>% 
  ggplot(aes(V1, V2, colour = region)) +
  geom_point() + 
  theme_minimal()


## ------------------------------------------------------------------------
# Let's go back to the olive data
decomp <- cancor(X, Y)
V <- X %*% decomp$xcoef
colnames(V) <- paste0("CC", seq_len(8))

library(lattice)
levelplot(cor(X, V[,1:2]))


## ------------------------------------------------------------------------
levelplot(cor(Y, V[,1:2]))

