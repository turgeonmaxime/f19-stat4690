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

