## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)
options(knitr.kable.NA = '-')


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


## ------------------------------------------------------------------------
# Recall our model
fit <- lm(cbind(tear, gloss, opacity) ~ rate, 
          data = Plastic)

new_x <- data.frame(rate = factor("High", 
                                  levels = c("Low", "High")))
(prediction <- predict(fit, newdata = new_x))


## ------------------------------------------------------------------------
X <- model.matrix(fit)
S <- crossprod(resid(fit))/(nrow(Plastic) - ncol(X))
new_x <- model.matrix(~rate, new_x)

quad_form <- drop(new_x %*% solve(crossprod(X)) %*% t(new_x))

# Estimation covariance
(est_cov <- S * quad_form) 

# Forecasting covariance
(fct_cov <- S *(1 + quad_form)) 


## ------------------------------------------------------------------------
# Estimation CIs
cbind(drop(prediction) - 1.96*diag(est_cov),
      drop(prediction) + 1.96*diag(est_cov))

# Forecasting CIs
cbind(drop(prediction) - 1.96*diag(fct_cov),
      drop(prediction) + 1.96*diag(fct_cov))


## ---- warning = FALSE----------------------------------------------------
# Going back to our example
full_model <- lm(cbind(tear, gloss, 
                       opacity) ~ rate*additive, 
                 data = Plastic)

anova(full_model, test = "Wilks") %>% 
  broom::tidy()  %>% 
  knitr::kable(digits = 3)

anova(full_model, test = "Roy") %>% 
  broom::tidy()  %>% 
  knitr::kable(digits = 3)


## ---- warning = FALSE----------------------------------------------------
# Fit a model with only rate
rate_model <- lm(cbind(tear, gloss, 
                       opacity) ~ rate, 
                 data = Plastic)

# Removing the dfs from approx
anova(full_model, rate_model,
      test = "Wilks") %>% 
  broom::tidy()  %>% 
  dplyr::select(-num.Df, -den.Df) %>% 
  knitr::kable(digits = 3)

anova(full_model, rate_model,
      test = "Roy") %>% 
  broom::tidy()  %>% 
  dplyr::select(-num.Df, -den.Df) %>% 
  knitr::kable(digits = 3)


## ------------------------------------------------------------------------
# Let's look at the eigenvalues
E <- crossprod(residuals(full_model))
H <- crossprod(residuals(rate_model)) - E

result <- eigen(H %*% solve(E),
                only.values = TRUE)
result$values[seq_len(2)]


## ----eval = -1-----------------------------------------------------------
AIC(full_model)
# Error in logLik.lm(full_model) : 
#   'logLik.lm' does not support multiple responses
class(full_model)


## ------------------------------------------------------------------------
logLik.mlm <- function(object, ...) {
  resids <- residuals(object)
  Sigma_ML <- crossprod(resids)/nrow(resids)
  ans <- sum(mvtnorm::dmvnorm(resids, 
                              sigma = Sigma_ML, 
                              log = TRUE))
  
  df <- prod(dim(coef(object))) + 
    choose(ncol(Sigma_ML) + 1, 2)
  attr(ans, c("nobs", "df")) <- c(nrow(resids), df)
  class(ans) <- "logLik"
  return(ans)
}


## ------------------------------------------------------------------------
logLik(full_model)

AIC(full_model)
AIC(rate_model)


## ---- message = FALSE----------------------------------------------------
library(openintro)
model <- lm(cbind(startPr, totalPr) ~ 
              nBids + cond + sellerRate + 
              wheels + stockPhoto, 
            data = marioKart)

X <- model.matrix(model)
P <- X %*% solve(crossprod(X)) %*% t(X)
lev_values <- diag(P)

hist(lev_values, 50)


## ------------------------------------------------------------------------
n <- nrow(HSB)
resids <- residuals(model)
S <- crossprod(resids)/(n - ncol(X))

S_inv <- solve(S)

const <- lev_values/((1 - lev_values)^2*ncol(X))
cook_values <- const * diag(resids %*% S_inv 
                            %*% t(resids))

hist(cook_values, 50)


## ------------------------------------------------------------------------
# Cut-off value
(cutoff <- qchisq(0.5, ncol(S)*(n - ncol(X))))
which(cook_values > cutoff)

