## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=FALSE)


## ------------------------------------------------------------------------
set.seed(1234)
n <- 100
# Generate uniform data
Y <- cbind(runif(n, -1, 1),
           runif(n, -1, 1))

# Check if it falls inside ellipse
Sigma <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
dists <- sqrt(diag(Y %*% solve(Sigma) %*%
                     t(Y)))
inside <- dists < 0.85


## ------------------------------------------------------------------------
# Plot points
colours <- c("red", "blue")[inside + 1]
plot(Y, col = colours)


## ------------------------------------------------------------------------
# Transform data
# (X, Y) -> (X^2, Y^2, sqrt(2)*X*Y)
Y_transf <- cbind(Y[,1]^2, Y[,2]^2,
                  sqrt(2)*Y[,1]*Y[,2])


## ------------------------------------------------------------------------
library(scatterplot3d)
scatterplot3d(Y_transf, color = colours,
              xlab = "X-axis",
              ylab = "Y-axis",
              zlab = "Z-axis")


## ---- echo = FALSE-------------------------------------------------------
scatterplot3d(Y_transf, color = colours,
              angle = 5,
              xlab = "X-axis",
              ylab = "Y-axis",
              zlab = "Z-axis")


## ------------------------------------------------------------------------
# Linear regression
outcome <- ifelse(inside, 1, -1)
head(outcome)


## ------------------------------------------------------------------------
model1 <- lm(outcome ~ Y)
pred1 <- sign(predict(model1))
table(outcome, pred1) # 67%


## ------------------------------------------------------------------------
model2 <- lm(outcome ~ Y_transf)
pred2 <- sign(predict(model2))
table(outcome, pred2) # 92%


## ---- echo = FALSE-------------------------------------------------------
par(mfrow = c(1, 2))

colours1 <- colours
colours1[outcome != pred1] <- "green"
plot(Y, col = colours1,
     main = "No transformation")

colours2 <- colours
colours2[outcome != pred2] <- "green"
plot(Y, col = colours2,
     main = "With transformation")
legend("bottomright", legend = "Mis.", col = "green", 
       pch = 1, bg = "white")


## ---- cache = TRUE, warning = FALSE, message = FALSE, eval = FALSE, echo = FALSE----
## library(tidytuesdayR)
## data_raw <- tt_load_gh("2019-10-15") %>%
##   tt_read_data("big_epa_cars.csv")


## ---- echo = FALSE, message = FALSE, eval = FALSE, echo = FALSE----------
## library(tidyverse)
## data <- data_raw %>%
##   mutate(
##     fuel         = paste0(fuelType1,",",fuelType2),
##     mpg_city     = paste0(city08 ,",",cityA08),
##     mpg_hw       = paste0(highway08 ,",",highwayA08),
##     c02          = paste0(co2,",",co2A),
##     trany        =
##       gsub("Auto\\(AM-S(\\d+)\\)","Automatic \\1-spd",
##       gsub("4-spd Doubled","4-spd",
##       gsub("(variable gear ratios)","variable_gear_ratios",
##                         trany)),perl=TRUE)
##   ) %>%
##   separate(trany,c("transmission","gears"),sep=" ") %>%
##   mutate(gears = gsub("-spd","",gears)) %>%
##   dplyr::select(
##     make         = make,
##     model        = model,
##     year         = year,
##     type         = VClass,
##     displacement = displ,
##     transmission,
##     gears,
##     cylinders    = cylinders,
##     drive,
##     fuel,
##     mpg_city,
##     mpg_hw,
##     c02
##   ) %>%
##   separate_rows(fuel,mpg_city,mpg_hw,c02,sep=",") %>%
##   filter(fuel     !="NA",
##          mpg_city !=0) %>%
##   mutate(mpg_city  = as.numeric(mpg_city),
##          mpg_hw    = as.numeric(mpg_hw),
##          c02       = as.numeric(c02),
##          c02       = na_if(c02,-1)) %>%
##   arrange(make,model,year)


## ---- message = FALSE----------------------------------------------------
library(ElemStatLearn)
library(tidyverse)

data_train <- prostate %>% 
  filter(train == TRUE) %>% 
  dplyr::select(-train)
data_test <- prostate %>% 
  filter(train == FALSE) %>% 
  dplyr::select(-train)


## ------------------------------------------------------------------------
model1 <- lm(lpsa ~ ., 
             data = data_train)
pred1 <- predict(model1, data_test)

mean((data_test$lpsa - pred1)^2)


## ---- message = FALSE----------------------------------------------------
# glmnet does lasso, elastic-net and
# ridge regression
library(glmnet)
X_train <- model.matrix(lpsa ~ . -  1, 
             data = data_train)
model2 <- glmnet(X_train, data_train$lpsa, 
                 alpha = 0, lambda = 0.7)


## ---- message = FALSE----------------------------------------------------
X_test <- model.matrix(lpsa ~ . -  1, 
             data = data_test)
pred2 <- predict(model2, X_test)

mean((data_test$lpsa - pred2)^2)


## ------------------------------------------------------------------------
# Let's start with the identity map for Phi
# We should get the same results as Ridge regression
X_train <- model.matrix(lpsa ~ ., 
                        data = data_train)
Y_train <- data_train$lpsa


## ------------------------------------------------------------------------
# Ridge regression
beta_hat <- solve(crossprod(X_train) + 
                    0.7*diag(ncol(X_train))) %*%
  t(X_train) %*% Y_train

beta_hat[1:3]


## ------------------------------------------------------------------------
# Dual problem
alpha_hat <- solve(tcrossprod(X_train) + 
                     0.7*diag(nrow(X_train))) %*% 
  Y_train

(t(X_train) %*% alpha_hat)[1:3]

all.equal(beta_hat, t(X_train) %*% alpha_hat)


## ------------------------------------------------------------------------
library(kernlab)
# Let's use the quadratic kernel
poly <- polydot(degree = 2)

Kmat <- kernelMatrix(poly, X_train)
Kmat[1:3, 1:3]


## ------------------------------------------------------------------------
alpha_poly <- solve(Kmat + 0.7*diag(nrow(X_train))) %*% 
  Y_train


## ------------------------------------------------------------------------
# Let's predict the test data
X_test <- model.matrix(lpsa ~ ., 
                       data = data_test)
k_pred <- kernelMatrix(poly, X_train, X_test)

pred_poly <- drop(t(alpha_poly) %*% k_pred)
mean((data_test$lpsa - pred_poly)^2)

# Compare with linear kernel
pred_lin <- drop(t(alpha_hat) %*% 
                   tcrossprod(X_train, X_test))
mean((data_test$lpsa - pred_lin)^2)


## ------------------------------------------------------------------------
# Now let's try a Gaussian kernel
# Note: Look at documentation for
# parametrisation
rbf <- rbfdot(sigma = 0.05)
Kmat <- kernelMatrix(rbf, X_train)

alpha_rbf <- solve(Kmat + 0.7*diag(nrow(X_train))) %*% 
  Y_train


## ------------------------------------------------------------------------
k_pred <- kernelMatrix(rbf, X_train, X_test)

pred_rbf <- drop(t(alpha_rbf) %*% k_pred)
mean((data_test$lpsa - pred_rbf)^2)


## ------------------------------------------------------------------------
# Can we do better by choosing a different sigma?
n <- nrow(X_train)

fit_rbf <- function(sigma) {
  rbf <- rbfdot(sigma = sigma)
  Kmat <- kernelMatrix(rbf, X_train)
  alpha_rbf <- solve(Kmat + 0.7*diag(n)) %*% 
    Y_train
  return(list(alpha = alpha_rbf,
              rbf = rbf))
}


## ------------------------------------------------------------------------
pred_rbf <- function(fit) {
  k_pred <- kernelMatrix(fit$rbf, X_train, 
                         X_test)
  pred_rbf <- drop(t(fit$alpha) %*% k_pred)
  return(pred_rbf)
}


## ------------------------------------------------------------------------
sigma_vect <- 10^seq(0, -2, by = -0.1)
MSE <- sapply(sigma_vect,
              function (sigma) {
                fit_rbf <- fit_rbf(sigma)
                pred_rbf <- pred_rbf(fit_rbf)
                mean((data_test$lpsa - pred_rbf)^2)
                })


## ------------------------------------------------------------------------
data.frame(
  sigma = sigma_vect,
  MSE = MSE
) %>% 
  ggplot(aes(sigma, MSE)) +
  geom_line() +
  theme_minimal() +
  scale_x_log10()


## ------------------------------------------------------------------------
data.frame(
  sigma = sigma_vect,
  MSE = MSE
) %>% 
  filter(MSE == min(MSE))


## ------------------------------------------------------------------------
library(caret)

# Blood-Brain barrier data
data(BloodBrain)
length(logBBB)
dim(bbbDescr)


## ------------------------------------------------------------------------
# 5-fold CV with sigma = 0.05
trainIndex <- createFolds(logBBB, k = 5)
str(trainIndex)


## ------------------------------------------------------------------------
# Let's redefine our functions from earlier
fit_rbf <- function(sigma, data_train, Y_train) {
  rbf <- rbfdot(sigma = sigma)
  Kmat <- kernelMatrix(rbf, data_train)
  alpha_rbf <- solve(Kmat + 
                       0.7*diag(nrow(data_train))) %*% 
    Y_train
  return(list(alpha = alpha_rbf, rbf = rbf))
}


## ------------------------------------------------------------------------
pred_rbf <- function(fit, data_train, data_test) {
  k_pred <- kernelMatrix(fit$rbf, data_train, 
                         data_test)
  pred_rbf <- drop(t(fit$alpha) %*% k_pred)
  return(pred_rbf)
}


## ------------------------------------------------------------------------
sapply(trainIndex, function(index){
         data_train <- bbbDescr[-index,] %>% 
           model.matrix( ~ ., data = .)
         Y_train <- logBBB[-index]
         data_test <- bbbDescr[index,] %>% 
           model.matrix( ~ ., data = .)
         fit_rbf <- fit_rbf(0.05, data_train, Y_train)
         pred_rbf <- pred_rbf(fit_rbf, data_train,
                              data_test)
         mean((logBBB[index] - pred_rbf)^2)
       }) -> MSEs


## ------------------------------------------------------------------------
MSEs
mean(MSEs)


## ------------------------------------------------------------------------
# Now, we can repeat for multiple sigmas
mse_calc <- function(sigma, data_train, data_test, 
                     Y_train, Y_test) {
  fit_rbf <- fit_rbf(sigma, data_train, Y_train)
  pred_rbf <- pred_rbf(fit_rbf, data_train,
                       data_test)
  mean((Y_test - pred_rbf)^2)
}


## ------------------------------------------------------------------------
sapply(trainIndex, function(index){
         data_train <- bbbDescr[-index,] %>% 
           model.matrix( ~ ., data = .)
         Y_train <- logBBB[-index]
         data_test <- bbbDescr[index,] %>% 
           model.matrix( ~ ., data = .)
         sapply(sigma_vect, mse_calc, 
                data_train = data_train, 
                data_test = data_test, 
                Y_train = Y_train,
                Y_test = logBBB[index])}) -> MSEs


## ------------------------------------------------------------------------
head(rowMeans(MSEs), n = 4)


## ------------------------------------------------------------------------
data.frame(
  sigma = sigma_vect,
  MSE = rowMeans(MSEs)
) %>% 
  ggplot(aes(sigma, MSE)) +
  geom_line() +
  theme_minimal() +
  scale_x_log10()


## ---- echo = FALSE, cache = TRUE-----------------------------------------
lambda_vect <- 10^seq(3, -2, length.out = 20)
sigma_vect <- 10^seq(0, -2, length.out = 20)

MSEs <- sapply(lambda_vect, function(lambda) {
  sapply(trainIndex, function(index){
    # Create data sets
    data_train <- bbbDescr[-index,] %>% 
      model.matrix( ~ ., data = .)
    Y_train <- logBBB[-index]
    data_test <- bbbDescr[index,] %>% 
      model.matrix( ~ ., data = .)
    Y_test <- logBBB[index]
    
    # Fit model and do 5-fold CV
    sapply(sigma_vect, function(sigma) {
      rbf <- rbfdot(sigma = sigma)
      Kmat <- kernelMatrix(rbf, data_train)
      alpha_rbf <- solve(Kmat + lambda*diag(nrow(data_train))) %*% Y_train
      k_pred <- kernelMatrix(rbf, data_train, data_test)
      pred_rbf <- drop(t(alpha_rbf) %*% k_pred)
      mean((Y_test - pred_rbf)^2)
      })
  }) -> temp
  rowMeans(temp)
})

rownames(MSEs) <- as.character(sigma_vect)

MSEs <- MSEs %>% 
  t %>% 
  as.data.frame() %>% 
  bind_cols(tibble(lambda = lambda_vect))


## ------------------------------------------------------------------------
# We can also tune lambda (see R code)
MSEs <- MSEs %>% 
  gather(sigma, MSE, -lambda) %>% 
  mutate(sigma = as.numeric(sigma)) 

head(MSEs, n = 3) 

MSEs %>% 
  ggplot(aes(lambda, MSE, group = sigma)) + 
  geom_line() +
  theme_minimal() +
  scale_x_log10() +
  geom_vline(xintercept = 0.7, linetype = 'dashed')


## ------------------------------------------------------------------------
MSEs %>% 
  filter(MSE == min(MSE))


## ---- echo = FALSE, warning = FALSE--------------------------------------
sapply(trainIndex, function(index){
         data_train <- bbbDescr[-index,]
         Y_train <- logBBB[-index]
         data_test <- bbbDescr[index,]
         fit_lm <- lm(Y_train ~ ., data = data_train)
         pred_lm <- predict(fit_lm, data_test)
         mean((logBBB[index] - pred_lm)^2)
       }) -> MSEs_linear

