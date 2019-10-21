## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ------------------------------------------------------------------------
library(psych)

dim(bfi)
names(bfi)


## ----message = FALSE-----------------------------------------------------
library(tidyverse)

data <- bfi %>% 
  select(-gender, -education, -age) %>% 
  filter(complete.cases(.))


## ------------------------------------------------------------------------
cor.plot(cor(data))


## ------------------------------------------------------------------------
decomp <- prcomp(data)
summary(decomp)$importance[,1:3]

cum_prop <- decomp %>% 
  summary %>% 
  .[["importance"]] %>% 
  .["Cumulative Proportion",]

which(cum_prop > 0.8)


## ------------------------------------------------------------------------
Lhat <- decomp$rotation[,1:14] %*% 
  diag(decomp$sdev[1:14])
Psi_hat <- diag(cov(data) - tcrossprod(Lhat))

sum((cov(data) - tcrossprod(Lhat) - diag(Psi_hat))^2)
sum(diag(cov(data)))

