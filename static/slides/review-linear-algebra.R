## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ------------------------------------------------------------------------
A <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
A

# Determinant
det(A)

# Rank
library(Matrix)
rankMatrix(A)

# Condition number
kappa(A)

# How to compute the trace?
sum(diag(A))

# Transpose
t(A)

# Inverse
solve(A)

A %*% solve(A) # CHECK


## ------------------------------------------------------------------------
A <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
B <- matrix(c(4, 3, 2, 1), nrow = 2, ncol = 2)

# Addition
A + B

# Scalar multiplication
3*A

# Matrix multiplication
A %*% B

# Hadamard product aka entrywise multiplication
A * B

# Matrix-vector product
vect <- c(1, 2)
A %*% vect

# BE CAREFUL: R recycles vectors
A * vect


## ------------------------------------------------------------------------
A <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

result <- eigen(A)

names(result)

result$values

result$vectors

1/sqrt(2)


## ------------------------------------------------------------------------
v1 <- c(1/sqrt(2), 1/sqrt(2))
v2 <- c(1/sqrt(2), -1/sqrt(2))

Lambda <- diag(c(1.5, 0.5))
P <- cbind(v1, v2)

P %*% Lambda %*% t(P)

# Now let's look at a random matrix----
A <- matrix(rnorm(3 * 3), ncol = 3, nrow = 3)
# Let's make it symmetric
A[lower.tri(A)] <- A[upper.tri(A)]
A

result <- eigen(A, symmetric = TRUE)
Lambda <- diag(result$values)
P <- result$vectors

P %*% Lambda %*% t(P)

# How to check if they are equal?
all.equal(A, P %*% Lambda %*% t(P))


## ------------------------------------------------------------------------
A <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

# Eigenvalue method
result <- eigen(A)
Lambda <- diag(result$values)
P <- result$vectors
A_sqrt <- P %*% Lambda^0.5 %*% t(P)

all.equal(A, A_sqrt %*% A_sqrt) # CHECK

# Cholesky method
# It's upper triangular!
(L <- chol(A))

all.equal(A, t(L) %*% L) # CHECK


## ------------------------------------------------------------------------
set.seed(123)

A <- matrix(rnorm(3*3), ncol = 3)
# Make A symmetric
A[lower.tri(A)] <- A[upper.tri(A)]

# Set initial value
v_current <- rnorm(3)
v_current <- v_current/norm(v_current, type = "2")


## ------------------------------------------------------------------------
# We'll perform 100 iterations
for (i in seq_len(100)) {
  # Save result from previous iteration
  v_previous <- v_current
  # Compute matrix product
  numerator <- A %*% v_current
  # Normalize
  v_current <- numerator/norm(numerator, type = "2")
}

v_current

# Corresponding eigenvalue
num <- t(v_current) %*% A %*% v_current
denom <- t(v_current) %*% v_current
num/denom

# CHECK results
result <- eigen(A, symmetric = TRUE)
result$values[which.max(abs(result$values))]
result$vectors[,which.max(abs(result$values))]


## ---- echo=FALSE, fig.height=6, fig.width=6------------------------------
set.seed(123)
p <- 2

A <- matrix(rnorm(p*p), ncol = p)
# Make A symmetric
A[lower.tri(A)] <- A[upper.tri(A)]

# Set initial value
v_current <- rnorm(p)
v_current <- v_current/norm(v_current, type = "2")
results <- matrix(NA, ncol = p, nrow = 100)

for (i in seq_len(100)) {
  results[i,] <- v_previous <- v_current
  # Compute matrix product
  numerator <- A %*% v_current
  # Normalize
  v_current <- numerator/norm(numerator, type = "2")
}

# Plot sequence
decomp <- eigen(A)
leading_vect <- decomp$vectors[,which.max(abs(decomp$values))]

theta_vect <- seq(0, 2*pi, length.out = 100)
x_circle <- cos(theta_vect)
y_circle <- sin(theta_vect)

plot(x_circle, y_circle, col = 'grey60', type = 'l',
     xlab = "", ylab = "")
lines(results[,1], results[,2])
points(results[,1], results[,2], pch = 19,
       col = c("green", rep("black", 98), "red"))
points(x = c(0, leading_vect[1]), y = c(0, leading_vect[2]), 
       col = 'blue', type = 'b', pch = 19, cex = 2)


## ---- cache = FALSE------------------------------------------------------
set.seed(1234)
A <- matrix(rnorm(3 * 2), ncol = 2, nrow = 3)
result <- svd(A)
names(result)

result$d
result$u
result$v

D <- diag(result$d)
all.equal(A, result$u %*% D %*% t(result$v)) #CHECK


## ---- cache = FALSE------------------------------------------------------
# Note: crossprod(A) == t(A) %*% A
#  tcrossprod(A) == A %*% t(A)
U <- eigen(tcrossprod(A))$vectors
V <- eigen(crossprod(A))$vectors

D <- matrix(0, nrow = 3, ncol = 2)
diag(D) <- result$d

all.equal(A, U %*% D %*% t(V)) # CHECK


## ---- cache = FALSE------------------------------------------------------
# What went wrong?
# Recall that eigenvectors are unique 
# only up to a sign!

# These elements should all be positive
diag(t(U) %*% A %*% V)

# Therefore we need to multiply the 
# corresponding columns of U or V 
# (but not both!) by -1
cols_flip <- which(diag(t(U) %*% A %*% V) < 0)
V[,cols_flip] <- -V[,cols_flip]

all.equal(A, U %*% D %*% t(V)) # CHECK

