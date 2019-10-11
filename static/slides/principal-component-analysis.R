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


## ---- message = FALSE----------------------------------------------------
library(dslabs)
library(ggridges)

# Data on Breast Cancer
as.data.frame(brca$x) %>% 
  gather(variable, measurement) %>% 
  mutate(variable = reorder(variable, measurement, 
                            median)) %>% 
  ggplot(aes(x = measurement, y = variable)) + 
  geom_density_ridges() + theme_ridges() +
  coord_cartesian(xlim = c(0, 250))


## ------------------------------------------------------------------------
# Remove some variables
rem_index <- which(colnames(brca$x) %in%
                     c("area_worst", "area_mean",
                       "perimeter_worst", 
                       "perimeter_mean"))
dataset <- brca$x[,-rem_index]
decomp <- prcomp(dataset)
summary(decomp)$importance[,1:3]


## ------------------------------------------------------------------------
screeplot(decomp, type = 'l')


## ------------------------------------------------------------------------
# Let's put a CI around the first eigenvalue
first_ev <- decomp$sdev[1]^2
n <- nrow(dataset)

# Recall that TV = 2166
c("LB" = first_ev/(1+qnorm(0.975)*sqrt(2/n)),
  "Est." = first_ev,
  "UP" = first_ev/(1-qnorm(0.975)*sqrt(2/n)))


## ------------------------------------------------------------------------
B <- 1000; n <- 100; p <- 3

results <- purrr::map_df(seq_len(B), function(b) {
    X <- matrix(rnorm(p*n, sd = sqrt(c(1, 2, 3))), 
                ncol = p, byrow = TRUE)
    tmp <- eigen(cov(X), symmetric = TRUE, 
                 only.values = TRUE)
    tibble(ev1 = tmp$values[1],
           ev2 = tmp$values[2],
           ev3 = tmp$values[3])
})


## ------------------------------------------------------------------------
results %>% 
  gather(ev, value) %>% 
  ggplot(aes(value, fill = ev)) + 
  geom_density(alpha = 0.5) +
  theme_minimal() +
  geom_vline(xintercept = c(1, 2, 3),
             linetype = 'dashed')


## ------------------------------------------------------------------------
results %>% 
  summarise_all(mean)


## ------------------------------------------------------------------------
p <- 2
results <- purrr::map_df(seq_len(B), function(b) {
  X <- matrix(rnorm(p*n, sd = c(1, 2)), ncol = p,
              byrow = TRUE)
  tmp <- eigen(cov(X), symmetric = TRUE)
  
  tibble(
    xend = tmp$vectors[1,1],
    yend = tmp$vectors[2,1]
  )
})


## ------------------------------------------------------------------------
results %>% 
  ggplot() +
  geom_segment(aes(xend = xend, yend = yend),
               x = 0, y = 0, colour = 'grey60') +
  geom_segment(x = 0, xend = 0,
               y = 0, yend = 1,
               colour = 'blue', size = 2) +
  expand_limits(y = 0, x = c(-1, 1)) +
  theme_minimal()


## ---- message = FALSE----------------------------------------------------
# Or looking at angles
results %>% 
  transmute(theta = atan2(yend, xend)) %>% 
  ggplot(aes(theta)) +
  geom_histogram() +
  theme_minimal() + 
  geom_vline(xintercept = pi/2)

