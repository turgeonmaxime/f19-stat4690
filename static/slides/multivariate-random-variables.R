## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ------------------------------------------------------------------------
A <- matrix(c(5, 4, 4, 5), ncol = 2)

results <- eigen(A, symmetric = TRUE,
                 only.values = TRUE)

# Generalized variance
prod(results$values)

# Total variance
sum(results$values)

# Compare this with the following
B <- matrix(c(5, -4, -4, 5), ncol = 2)

# Generalized variance
# GV(A) = 9
det(B)

# Total variance
# TV(A) = 10
sum(diag(B))


## ------------------------------------------------------------------------
Sigma <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)

# First create a circle
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- cbind(cos(theta_vect), sin(theta_vect))
# Then turn into ellipse
ellipse <- circle %*% Sigma


## ------------------------------------------------------------------------
# Principal axes
result <- eigen(Sigma, symmetric = TRUE)

first <- result$values[1]*result$vectors[,1]
second <- result$values[2]*result$vectors[,2]


## ------------------------------------------------------------------------
# Plot results
plot(ellipse, type = 'l')
lines(x = c(0, first[1]),
      y = c(0, first[2]))
lines(x = c(0, second[1]),
      y = c(0, second[2]))


## ------------------------------------------------------------------------
# Generalized Variance
det(Sigma)

# Predicted volume of the ellipse above
pi/(gamma(1))*sqrt(det(Sigma))


## ------------------------------------------------------------------------
# How can we estimate the area?
# Monte Carlo simulation!
Sigma_inv <- solve(Sigma)

x_1 <- runif(1000, min = min(ellipse[,1]),
             max = max(ellipse[,1]))
x_2 <- runif(1000, min = min(ellipse[,2]),
             max = max(ellipse[,2]))

X <- cbind(x_1, x_2)
distances <- apply(X, 1, function(row) {
  sqrt(t(row) %*% Sigma_inv %*% row)
  })


## ------------------------------------------------------------------------
# Estimate
length_x <- diff(range(ellipse[,1])) 
length_y <- diff(range(ellipse[,2]))
area_rect <- length_x * length_y

area_rect * mean(distances <= 1)

