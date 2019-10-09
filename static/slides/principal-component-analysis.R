## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=FALSE)


## ----echo = FALSE--------------------------------------------------------
C <- chol(S <- toeplitz(.9 ^ (0:31))) # Cov.matrix and its root
set.seed(17)
X <- matrix(rnorm(32000), 1000, 32)
Z <- X %*% C  ## ==>  cov(Z) ~=  C'C = S
pZ <- prcomp(Z, tol = 0.1)
screeplot(pZ, type = 'lines', main = '')


## ------------------------------------------------------------------------
library(mvtnorm)
Sigma <- matrix(c(1, 0.5, 0.1,
                  0.5, 1, 0.5,
                  0.1, 0.5, 1),
                ncol = 3) 

set.seed(17)
X <- rmvnorm(100, sigma = Sigma)
pca <- prcomp(X)


## ------------------------------------------------------------------------
summary(pca)


## ------------------------------------------------------------------------
screeplot(pca, type = 'l')


## ------------------------------------------------------------------------
pca <- prcomp(USArrests, scale = TRUE)


## ------------------------------------------------------------------------
summary(pca)


## ------------------------------------------------------------------------
screeplot(pca, type = 'l')


## ----message = FALSE-----------------------------------------------------
library(ElemStatLearn)
library(tidyverse)
train <- subset(prostate, train == TRUE,
                select = -train)
test  <- subset(prostate, train == FALSE,
                select = -train)

# First model: Linear regression
lr_model <- lm(lpsa ~ ., data = train)
lr_pred <- predict(lr_model, newdata = test)
(lr_mse <- mean((test$lpsa - lr_pred)^2))

# PCA
decomp <- train %>%
  subset(select = -lpsa) %>%
  as.matrix() %>%
  prcomp
summary(decomp)$importance[,1:3]

screeplot(decomp, type = 'lines')


## ------------------------------------------------------------------------
# Second model: PCs for predictors
train_pc <- train
train_pc$PC1 <- decomp$x[,1]
pc_model <- lm(lpsa ~ PC1, data = train_pc)


## ------------------------------------------------------------------------
test_pc <- as.data.frame(predict(decomp, test))
pc_pred <- predict(pc_model,
                   newdata = test_pc)
(pc_mse <- mean((test$lpsa - pc_pred)^2))


## ------------------------------------------------------------------------
contribution <- decomp$rotation[,"PC1"]
round(contribution, 3)[1:6]
round(contribution, 3)[7:8]

(keep <- names(which(abs(contribution) > 0.01)))

fs_model <- lm(lpsa ~ ., data = train[,c(keep, "lpsa")])
fs_pred <- predict(fs_model, newdata = test)
(fs_mse <- mean((test$lpsa - fs_pred)^2))


## ---- message = FALSE----------------------------------------------------
model_plot <- data.frame(
  "obs" = test$lpsa,
  "LR" = lr_pred,
  "PC" = pc_pred,
  "FS" = fs_pred
) %>%
  gather(Model, pred, -obs)


## ---- message = FALSE----------------------------------------------------
ggplot(model_plot,
       aes(pred, obs, colour = Model)) +
  geom_point() +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0) +
  theme(legend.position = 'top') +
  xlab("Predicted") + ylab("Observed")


## ----mnist, cache = TRUE-------------------------------------------------
library(dslabs)

mnist <- read_mnist()

dim(mnist$train$images)
dim(mnist$test$images)

head(mnist$train$labels)


## ------------------------------------------------------------------------
matrix(mnist$train$images[1,], ncol = 28) %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE)


## ----cache = TRUE--------------------------------------------------------
decomp <- prcomp(mnist$train$images)


## ------------------------------------------------------------------------
screeplot(decomp, type = 'lines',
          npcs = 20, main = "")


## ------------------------------------------------------------------------
decomp$x[,1:2] %>%
  as.data.frame() %>%
  mutate(label = factor(mnist$train$labels)) %>%
  ggplot(aes(PC1, PC2, colour = label)) +
  geom_point(alpha = 0.5) +
  theme_minimal()


## ------------------------------------------------------------------------
# And on the test set
decomp %>%
  predict(newdata = mnist$test$images) %>%
  as.data.frame() %>%
  mutate(label = factor(mnist$test$labels)) %>%
  ggplot(aes(PC1, PC2, colour = label)) +
  geom_point(alpha = 0.5) +
  theme_minimal()


## ------------------------------------------------------------------------
par(mfrow = c(2, 2))
for (i in seq_len(4)) {
  matrix(decomp$rotation[,i], ncol = 28) %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, main = paste0("PC", i))
}


## ------------------------------------------------------------------------
# Approximation with 90 PCs
approx_mnist <- decomp$rotation[, seq_len(90)] %*% 
  decomp$x[1, seq_len(90)]
par(mfrow = c(1, 2))

matrix(mnist$train$images[1,], ncol = 28) %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, main = "Original")
matrix(approx_mnist, ncol = 28) %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, main = "Approx")


## ------------------------------------------------------------------------
library(mvtnorm)
Sigma <- matrix(c(1, 0.5, 0.1,
                  0.5, 1, 0.5,
                  0.1, 0.5, 1),
                ncol = 3) 
mu <- c(1, 2, 2)


## ------------------------------------------------------------------------
set.seed(17)
X <- rmvnorm(100, mean = mu,
             sigma = Sigma)
pca <- prcomp(X)

colMeans(X)
colMeans(pca$x)


## ------------------------------------------------------------------------
# On the other hand
pca <- prcomp(X, center = FALSE)
colMeans(pca$x)


## ------------------------------------------------------------------------
set.seed(1234)
# Random measurement error
sigma <- 5

# Exact relationship between 
# Celsius and Fahrenheit
temp_c <- seq(-40, 40, by = 1)
temp_f <- 1.8*temp_c + 32


## ------------------------------------------------------------------------
# Add measurement error
temp_c_noise <- temp_c + rnorm(n = length(temp_c), 
                               sd = sigma)
temp_f_noise <- temp_f + rnorm(n = length(temp_f), 
                               sd = sigma)


## ------------------------------------------------------------------------
# Linear model
(fit <- lm(temp_f_noise ~ temp_c_noise))
confint(fit)

# PCA
pca <- prcomp(cbind(temp_c_noise, temp_f_noise))
pca$rotation
pca$rotation[2,"PC1"]/pca$rotation[1,"PC1"]

