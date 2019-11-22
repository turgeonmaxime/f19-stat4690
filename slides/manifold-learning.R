## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache = FALSE, message = FALSE)


## -----------------------------------------------------------------------------
n <- 1000
F1 <- runif(n, 0, 10)
F2 <- runif(n, -1, 1)


## -----------------------------------------------------------------------------
X <- F1 * cos(F1)
Y <- F2
Z <- F1 * sin(F1)


## -----------------------------------------------------------------------------
library(scatterplot3d)
library(colorspace)
colours <- cut(F1, breaks = seq(0, 10), 
               labels = diverging_hcl(10))
par(mfrow = c(1, 2))
scatterplot3d(X, Y, Z, pch = 19, asp = 1, 
              color = colours)
scatterplot3d(X, Y, Z, pch = 19, asp = 1, 
              color = colours, angle = 80)


## -----------------------------------------------------------------------------
# Let's see if PCA can unroll the Swiss roll
decomp <- prcomp(cbind(X, Y, Z))

plot(decomp$x[,1:2],
     col = as.character(colours), pch = 19)


## ---- message = FALSE---------------------------------------------------------
library(dslabs)
library(tidyverse)

mnist <- read_mnist()

data <- mnist$train$images[mnist$train$labels == 2, ]


## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
# With crossing
matrix(data[1,], ncol = 28)[ , 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, asp = 1)
# Without crossing
matrix(data[4,], ncol = 28)[ , 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, asp = 1)


## -----------------------------------------------------------------------------
decomp <- prcomp(data)
decomp$x[,1:2] %>%
  as.data.frame() %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(alpha = 0.5) +
  theme_minimal()


## -----------------------------------------------------------------------------
# First PC
par(mfrow = c(1, 2))
index_right <- which.max(decomp$x[,1])
matrix(data[index_right,], ncol = 28)[ , 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, asp = 1)
index_left <- which.min(decomp$x[,1])
matrix(data[index_left,], ncol = 28)[ , 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, asp = 1)


## -----------------------------------------------------------------------------
# Second PC
par(mfrow = c(1, 2))
index_top <- which.max(decomp$x[,2])
matrix(data[index_top,], ncol = 28)[ , 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, asp = 1)
index_bottom <- which.min(decomp$x[,2])
matrix(data[index_bottom,], ncol = 28)[ , 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, asp = 1)


## ---- echo = FALSE------------------------------------------------------------
k <- 15
par(mfrow = c(3, 5))
for (i in seq(0, k - 1)) {
  index <- pmin(round(nrow(data)/(k - 1)) * i + 1, nrow(data))
  pc_val <- sort(decomp$x[, 1])[index]
  matrix(data[decomp$x[, 1] == pc_val,], ncol = 28)[ , 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, asp = 1, main = paste0("PC1=", round(pc_val)))
}


## ---- echo = FALSE------------------------------------------------------------
indices <- pmin(round(nrow(data)/(k - 1)) * seq(0, k - 1) + 1, nrow(data))
pc_vals <- sort(decomp$x[, 1])[indices]

data_samp <- decomp$x[,1:2] %>%
  as.data.frame() %>%
  filter(PC1 %in% pc_vals)

decomp$x[,1:2] %>%
  as.data.frame() %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(alpha = 0.1) +
  geom_point(data = data_samp, colour = "black", 
             fill = "red", pch = 21) +
  theme_minimal()


## ---- cache = FALSE, message = FALSE------------------------------------------
library(dimRed)


## ----message = TRUE-----------------------------------------------------------
isomap_sr <- embed(cbind(X, Y, Z), "Isomap", knn = 10,
                   ndim = 2)


## -----------------------------------------------------------------------------
isomap_sr@data@data %>% 
  plot(col = as.character(colours), pch = 19)


## ----isomap, cache = TRUE, message = TRUE-------------------------------------
isomap_res <- embed(data, "Isomap", knn = 10,
                    ndim = 2)


## -----------------------------------------------------------------------------
isomap_res@data %>% 
  as.data.frame() %>% 
  ggplot(aes(iso.1, iso.2)) +
  geom_point(alpha = 0.5) +
  theme_minimal()


## ---- echo = FALSE, eval = TRUE-----------------------------------------------
k <- 15
par(mfrow = c(3, 5))
for (i in seq(0, k - 1)) {
  index <- pmin(round(nrow(data)/(k - 1)) * i + 1, nrow(data))
  iso_val <- sort(isomap_res@data@data[, 1])[index]
  matrix(data[isomap_res@data@data[, 1] == iso_val,], ncol = 28)[ , 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, asp = 1, main = paste0("ISO1=", round(iso_val)))
}


## -----------------------------------------------------------------------------
n <- 1000
data_knn <- data[seq_len(n),]

# Compute distances
Delta <- dist(data_knn)

# Take 5-NN to first obs.
neighbours <- order(Delta[seq_len(n-1)])[1:5]


## -----------------------------------------------------------------------------
main_obs <- data_knn[1,]
nb_data <- data_knn[neighbours,]


## ----echo = FALSE-------------------------------------------------------------
par(mfrow = c(2,3))
for (i in 1:5) {
  matrix(nb_data[i,], ncol = 28)[, 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, main = glue::glue("NN {i}"), asp = 1)
}


## -----------------------------------------------------------------------------
# Center neighbours around main obs
nb_data_c <- sweep(nb_data, 2, main_obs)
# Local cov matrix
Gmat <- tcrossprod(nb_data_c)


## -----------------------------------------------------------------------------
# Find weights
w_vect <- solve(Gmat, rep(1, 5))
w_vect <- w_vect/sum(w_vect)


## -----------------------------------------------------------------------------
# Compare original with approx.
approx <- drop(w_vect %*% nb_data)
par(mfrow = c(1, 2))

matrix(main_obs, ncol = 28)[, 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, main = "Original", asp = 1)
matrix(approx, ncol = 28)[, 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, main = "Approx. (5 NN)", asp = 1)


## ---- echo = FALSE------------------------------------------------------------
# Take 50-NN to first obs.
K <- 25
neighbours <- order(Delta[1:99])[seq_len(K)]
# Center neighbours around main obs
nb_data <- data_knn[neighbours,]
nb_data_c <- sweep(nb_data, 2, main_obs)
# Local cov matrix
Gmat <- tcrossprod(nb_data_c)
# Find weights
w_vect <- solve(Gmat, rep(1, K))
w_vect <- w_vect/sum(w_vect)
# Compare original with approx.
approx <- drop(w_vect %*% nb_data)
par(mfrow = c(1, 2))

matrix(main_obs, ncol = 28)[, 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, main = "Original", asp = 1)
matrix(approx, ncol = 28)[, 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, main = glue::glue("Approx. ({K} NN)"), asp = 1)


## ----message = TRUE-----------------------------------------------------------
lle_sr <- embed(cbind(X, Y, Z), "LLE", knn = 20,
                ndim = 2)


## -----------------------------------------------------------------------------
lle_sr@data@data %>% 
  plot(col = as.character(colours), pch = 19)


## ----lle, cache = TRUE, message = TRUE----------------------------------------
lle_res <- embed(data, "LLE", knn = 50,
                 ndim = 2)


## -----------------------------------------------------------------------------
lle_res@data %>%
  as.data.frame() %>%
  ggplot(aes(LLE1, LLE2)) +
  geom_point(alpha = 0.5) +
  theme_minimal()


## ---- echo = FALSE, eval = TRUE-----------------------------------------------
k <- 15
par(mfrow = c(3, 5))
for (i in seq(0, k - 1)) {
  index <- pmin(round(nrow(data)/(k - 1)) * i + 1, nrow(data))
  lle_val <- sort(lle_res@data@data[, 1])[index]
  matrix(data[lle_res@data@data[, 1] == lle_val,], ncol = 28)[ , 28:1] %>%
  image(col = gray.colors(12, rev = TRUE),
        axes = FALSE, asp = 1, main = paste0("ISO1=", round(lle_val)))
}

