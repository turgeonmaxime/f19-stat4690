## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ------------------------------------------------------------------------
library(psych)

dim(bfi)
tail(names(bfi), n = 3)


## ----message = FALSE-----------------------------------------------------
library(tidyverse)

# Remove demographic variable and keep complete data
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

which(cum_prop > 0.8)[1]


## ------------------------------------------------------------------------
Lhat <- decomp$rotation[,1:14] %*% 
  diag(decomp$sdev[1:14])
Psi_hat <- diag(cov(data) - tcrossprod(Lhat))


## ------------------------------------------------------------------------
# Sum squared error
sum((cov(data) - tcrossprod(Lhat) - diag(Psi_hat))^2)
# Compare to the total variance
sum(diag(cov(data)))

# Our FA model explains:
sum(colSums(Lhat^2)/sum(diag(cov(data))))


## ------------------------------------------------------------------------
(m <- sum(eigen(cor(data))$values > 1))

Lhat <- decomp$rotation[,seq_len(m)] %*% 
  diag(decomp$sdev[seq_len(m)])
Psi_hat <- diag(cov(data) - tcrossprod(Lhat))

# Sum squared error
sum((cov(data) - tcrossprod(Lhat) - diag(Psi_hat))^2)
# Compare to the total variance
sum(diag(cov(data)))


## ------------------------------------------------------------------------
# Our FA model explains:
sum(colSums(Lhat^2)/sum(diag(cov(data))))


## ------------------------------------------------------------------------
# We can also  visualize the fit
Sn <- cov(data)
Sn_fit <- tcrossprod(Lhat) + diag(Psi_hat)

library(lattice)
levelplot(Sn - Sn_fit)


## ------------------------------------------------------------------------
# Or if you prefer % difference
levelplot((Sn - Sn_fit)/Sn)


## ------------------------------------------------------------------------
# CAREFUL: uses correlation matrix
fa_decomp <- factanal(data, factors = m,
                      rotation = 'none')

# Uniquenesses are the diagonal elements 
# of the matrix Psi
Psi_mle <- fa_decomp$uniquenesses
Lmle <- fa_decomp$loadings


## ------------------------------------------------------------------------
# We get an estimate of the correlation
R_mle <- tcrossprod(Lmle) + diag(Psi_mle)

sum((cor(data) - R_mle)^2)

# Our FA model explains:
sum(colSums(Lmle^2)/ncol(data))


## ------------------------------------------------------------------------
levelplot(cor(data) - R_mle)


## ------------------------------------------------------------------------
# To factor the covariance matrix
# Use psych::fa
fa_decomp <- psych::fa(data, nfactors = m,
                       rotate = "none",
                       covar = TRUE,
                       fm = "ml")

# Extract estimates
Psi_mle <- fa_decomp$uniquenesses
Lmle <- fa_decomp$loadings


## ------------------------------------------------------------------------
# We get an estimate of the covariance
Sn_mle <- tcrossprod(Lmle) + diag(Psi_mle)

sum((Sn - Sn_mle)^2)

# Our FA model explains:
sum(colSums(Lmle^2)/sum(diag(Sn)))


## ------------------------------------------------------------------------
levelplot(Sn - Sn_mle)


## ------------------------------------------------------------------------
# Compare MLE with PC estimate
levelplot(Sn_fit - Sn_mle)


## ------------------------------------------------------------------------
# Let's start with m=2 for visualization
fa_decomp <- factanal(data, factors = 2,
                      rotation = 'none')

initial_loadings <- fa_decomp$loadings
varimax_loadings <- varimax(initial_loadings)


## ----echo = FALSE--------------------------------------------------------
par(mfrow = c(1, 2))

plot(initial_loadings, type = 'n',
     xlim = c(-0.5, 0.8),
     ylim = c(-0.6, 0.6),
     main = "Initial")
text(initial_loadings, label = rownames(initial_loadings))
abline(h = 0, v = 0)

plot(varimax_loadings$loadings, type = 'n',
     xlim = c(-0.5, 0.8),
     ylim = c(-0.6, 0.6),
     main = "Varimax")
text(varimax_loadings$loadings, 
     label = rownames(varimax_loadings$loadings))
abline(h = 0, v = 0)

par(mfrow = c(1,1))


## ------------------------------------------------------------------------
# You can extract the matrix T
varimax_loadings$rotmat

# We can also get the angle of rotation
acos(varimax_loadings$rotmat[1,1])


## ------------------------------------------------------------------------
# In more dimensions
fa_decomp <- factanal(data, factors = m,
                      rotation = 'varimax')

levelplot(unclass(fa_decomp$loadings),
          xlab = "", ylab = "")

