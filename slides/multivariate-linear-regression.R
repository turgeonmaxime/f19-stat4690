## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ---- message = FALSE----------------------------------------------------
# Let's revisit the plastic film data
library(heplots)
library(tidyverse)

Y <- Plastic %>% 
  select(tear, gloss, opacity) %>% 
  as.matrix

X <- model.matrix(~ rate, data = Plastic)
head(X)


## ------------------------------------------------------------------------
(B_hat <- solve(crossprod(X)) %*% t(X) %*% Y)


## ------------------------------------------------------------------------
# Compare with lm output
fit <- lm(cbind(tear, gloss, opacity) ~ rate, 
          data = Plastic)
coef(fit)


## ------------------------------------------------------------------------
Y_hat <- fitted(fit)
residuals <- residuals(fit)

crossprod(Y_hat, residuals)
crossprod(X, residuals)

# Is this really zero?
isZero <- function(mat) {
  all.equal(mat, matrix(0, ncol = ncol(mat), 
                        nrow = nrow(mat)),
            check.attributes = FALSE)
}

isZero(crossprod(Y_hat, residuals))
isZero(crossprod(X, residuals))


## ---- message = FALSE----------------------------------------------------
library(candisc)

dataset <- HSB[,c("math", "sci")]

(corr_est <- cor(dataset)[1,2])

# Choose a number of bootstrap samples
B <- 5000
corr_boot <- replicate(B, {
  samp_boot <- sample(nrow(dataset),
                      replace = TRUE)
  dataset_boot <- dataset[samp_boot,]
  cor(dataset_boot)[1,2]
})

quantile(corr_boot,
         probs = c(0.025, 0.975))


## ---- message = FALSE----------------------------------------------------
hist(corr_boot, breaks = 50)
abline(v = corr_est, col = 'red',
       lty = 2)


## ---- message=FALSE------------------------------------------------------
B_boot <- replicate(B, {
  samp_boot <- sample(nrow(Y),
                      replace = TRUE)
  X_boot <- X[samp_boot,]
  Y_boot <- Y[samp_boot,]
  
  solve(crossprod(X_boot)) %*% t(X_boot) %*% Y_boot
})

# The output is a 3-dim array
dim(B_boot)

B_boot[,,1]

# CI for effect of rate on tear
quantile(B_boot["rateHigh", "tear",],
         probs = c(0.025, 0.975))

# CI for effect of rate on gloss
quantile(B_boot["rateHigh", "gloss",],
         probs = c(0.025, 0.975))

# CI for effect of rate on opacity
quantile(B_boot["rateHigh", "opacity",],
         probs = c(0.025, 0.975))


## ---- message=FALSE, eval = FALSE----------------------------------------
## library(ggforce)
## 
## B_boot["rateHigh",,] %>%
##   t() %>%
##   as.data.frame() %>%
##   ggplot(aes(x = .panel_x, y = .panel_y)) +
##   geom_point() +
##   geom_autodensity() +
##   geom_density2d() +
##   facet_matrix(vars(everything()),
##                layer.diag = 2,
##                layer.upper = 3)


## ---- message=FALSE, echo = FALSE----------------------------------------
# Need to repeat code above with echo = FALSE
# In order to get the plot on a separate slide
library(ggforce)

B_boot["rateHigh",,] %>% 
  t() %>% 
  as.data.frame() %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_autodensity() +
  geom_density2d() +
  facet_matrix(vars(everything()),
               layer.diag = 2,
               layer.upper = 3)


## ------------------------------------------------------------------------
# There is some correlation, but not much
B_boot["rateHigh",,] %>% 
  t() %>% 
  cor()

