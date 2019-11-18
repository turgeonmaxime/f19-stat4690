## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=FALSE,
                      message = FALSE)


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


## ------------------------------------------------------------------------
# First, let's look at the average eigenvalue
decomp <- eigen(Sn, symmetric = TRUE,
                only.values = TRUE)

mean(decomp$values)

sum(decomp$values > mean(decomp$values))


## ------------------------------------------------------------------------
# We will go from 1 factor to 15
prop_var_explained <- purrr::map_df(
  seq(1, 15),  function(m) {
    fa_decomp <- factanal(data, factors = m)
    Lmle <- fa_decomp$loadings
    prop <- sum(colSums(Lmle^2))/ncol(data)
    
    data.frame(
      prop = prop,
      m = m
    )})


## ------------------------------------------------------------------------
prop_var_explained %>% 
  ggplot(aes(m, prop)) + 
  geom_point() + 
  geom_line() + 
  theme_minimal() +
  expand_limits(y = 0) +
  geom_vline(xintercept = 6,
             linetype = 'dotted')


## ----warning=FALSE-------------------------------------------------------
inform_crit <- purrr::map_df(
  seq(1, 15),  function(m) {
    fa_decomp <- psych::fa(data, nfactors = m,
                           fm = 'ml')
    data.frame(
      BIC = fa_decomp$BIC, eBIC = fa_decomp$EBIC,
      SABIC = fa_decomp$SABIC, m = m
    )})


## ------------------------------------------------------------------------
inform_crit %>% 
  gather(Criteria, value, -m) %>% 
  ggplot(aes(m, value, colour = Criteria)) + 
  geom_point() + 
  geom_line() + 
  theme_minimal() +
  theme(legend.position = 'top') +
  geom_vline(xintercept = 6,
             linetype = 'dotted')


## ------------------------------------------------------------------------
inform_crit %>% 
  gather(Criteria, value, -m) %>% 
  group_by(Criteria) %>% 
  filter(value == min(value)) %>% 
  select(-value) %>% 
  spread(Criteria, m)


## ------------------------------------------------------------------------
m <- 6
fa_decomp <- psych::fa(data, nfactors = m,
                       covar = TRUE,
                       fm = 'ml')


## ------------------------------------------------------------------------
Sn <- cov(data)
Sn_mle <- tcrossprod(fa_decomp$loadings) + 
  diag(fa_decomp$uniquenesses)

test_stat <- nrow(data) * (
  log(det(Sn_mle)) -
    log(det(Sn))
)


## ------------------------------------------------------------------------
p <- ncol(data)
nu <- 0.5*((p-m)^2 - p - m)

test_stat > qchisq(0.95, df = nu)


## ---- message = FALSE----------------------------------------------------
m <- 6
fa_decomp <- psych::fa(data, nfactors = m,
                       covar = TRUE,
                       fm = 'ml')


## ---- message = FALSE----------------------------------------------------
# Extract estimates
InvPsi <- diag(fa_decomp$uniquenesses^-1)
Lhat <- fa_decomp$loadings

hat_mat <- solve(t(Lhat) %*% InvPsi %*% Lhat) %*%
  t(Lhat) %*% InvPsi

scores <- scale(data, center = TRUE,
                scale = FALSE) %*% t(hat_mat)


## ----message = FALSE-----------------------------------------------------
GGally::ggpairs(as.data.frame(scores))


## ---- message = FALSE----------------------------------------------------
scores_reg <- scale(data, center = TRUE,
                    scale = FALSE) %*% 
  solve(Sn) %*% Lhat


## ----message = FALSE-----------------------------------------------------
GGally::ggpairs(as.data.frame(scores_reg))


## ------------------------------------------------------------------------
# Let's look at agreement
round(cor(scores, scores_reg), 2)


## ------------------------------------------------------------------------
# Or graphically
plot(as.vector(scores),
     as.vector(scores_reg),
     xlab = "WLS",
     ylab = "Regression")
abline(a = 0, b = 1, 
       col='red',
       lty = 2)


## ---- echo = FALSE-------------------------------------------------------
inner_join(
  as.data.frame(scores_reg) %>% 
    mutate(ID = row_number()) %>% 
    gather(Factor, Score, -ID) %>% 
    rename(Regression = Score),
  as.data.frame(scores) %>% 
    mutate(ID = row_number()) %>% 
    gather(Factor, Score, -ID) %>% 
    rename(WLS = Score),
  by = c("ID", "Factor")
) %>% 
  ggplot(aes(WLS, Regression, colour = Factor)) + 
  geom_abline(slope = 1, intercept = 0,
              linetype = 'dashed') +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_minimal() +
  theme(legend.position = 'top')


## ------------------------------------------------------------------------
state.x77[1:5,1:4]


## ------------------------------------------------------------------------
library(GGally)
data <- as.data.frame(state.x77)
ggpairs(data)


## ------------------------------------------------------------------------
# Potential outliers?
data %>% 
  rownames_to_column('State') %>%
  top_n(Population, n = 2)

data %>% 
  rownames_to_column('State') %>%
  top_n(Area, n = 2)


## ------------------------------------------------------------------------
# 1. Principal Component Factor analysis
decomp <- princomp(data, cor = TRUE)
decomp


## ------------------------------------------------------------------------
biplot(decomp)


## ------------------------------------------------------------------------
m <- 3
Lhat <- decomp$loadings[,seq_len(m)] %*% 
  diag(decomp$sdev[seq_len(m)])
colnames(Lhat) <- paste0("PC", 1:3)
Psi_hat <- diag(cor(data) - tcrossprod(Lhat))

# Our FA model explains:
sum(colSums(Lhat^2))/ncol(data)


## ------------------------------------------------------------------------
# We can also  visualize the fit
Rn <- cor(data)
Rn_fit <- tcrossprod(Lhat) + diag(Psi_hat)

levelplot(Rn - Rn_fit)


## ------------------------------------------------------------------------
scores_pc <- scale(data, center = TRUE,
                    scale = TRUE) %*% 
  solve(Rn) %*% Lhat


## ------------------------------------------------------------------------
ggpairs(as.data.frame(scores_pc))


## ------------------------------------------------------------------------
# What state has the outlying values?
scores_pc %>% 
  as.data.frame %>% 
  rownames_to_column('State') %>%
  filter(PC2 > 4 | PC3 < -4)


## ------------------------------------------------------------------------
# 2. Maximum Likelihood Factor analysis
fa_decomp <- factanal(data, factors = m)


## ------------------------------------------------------------------------
# Extract estimates
Psi_mle <- fa_decomp$uniquenesses
Lmle <- fa_decomp$loadings

# Our FA model explains:
sum(colSums(Lmle^2))/ncol(data)


## ------------------------------------------------------------------------
# We can also  visualize the fit
Rn_mle <- tcrossprod(Lmle) + diag(Psi_mle)

levelplot(Rn - Rn_mle)


## ------------------------------------------------------------------------
scores_mle <- scale(data, center = TRUE,
                    scale = TRUE) %*% 
  solve(Rn) %*% Lmle


## ------------------------------------------------------------------------
ggpairs(as.data.frame(scores_mle))


## ------------------------------------------------------------------------
# 3. Compare both loadings
round(cor(scores_pc, scores_mle), 2)


## ------------------------------------------------------------------------
levelplot(Lhat,
          xlab = "", ylab = "")


## ------------------------------------------------------------------------
# Let's rotate the PC loadings
Lhat <- varimax(Lhat)$loadings
scores_pc <- scale(data, center = TRUE,
                    scale = TRUE) %*% 
  solve(Rn) %*% Lhat


## ------------------------------------------------------------------------
# Compare both loadings again
round(cor(scores_pc, scores_mle), 2)


## ------------------------------------------------------------------------
levelplot(unclass(Lhat),
          xlab = "", ylab = "")


## ------------------------------------------------------------------------
levelplot(unclass(Lmle),
          xlab = "", ylab = "")


## ---- echo = FALSE, eval = FALSE-----------------------------------------
## inner_join(
##   as.data.frame(scores_pc) %>%
##     mutate(ID = row_number()) %>%
##     gather(Factor, Score, -ID) %>%
##     rename(PC = Score) %>%
##     mutate(Factor = stringr::str_replace(Factor, "PC", "Factor")),
##   as.data.frame(scores_mle) %>%
##     mutate(ID = row_number()) %>%
##     gather(Factor, Score, -ID) %>%
##     rename(MLE = Score),
##   by = c("ID", "Factor")
## ) %>%
##   ggplot(aes(PC, MLE, colour = Factor)) +
##   geom_abline(slope = 1, intercept = 0,
##               linetype = 'dashed') +
##   geom_point(alpha = 0.5) +
##   geom_smooth(method = 'lm', se = FALSE) +
##   theme_minimal() +
##   theme(legend.position = 'top')


## ------------------------------------------------------------------------
# 4. Compare different m
# Let's start with scree plot
screeplot(decomp, type = 'l')


## ----warning = FALSE-----------------------------------------------------
# Then let's look at information criteria
inform_crit <- purrr::map_df(
  seq_len(4),  function(m) {
    fa_decomp <- psych::fa(data, nfactors = m,
                           fm = 'ml')
    data.frame(
      BIC = fa_decomp$BIC, eBIC = fa_decomp$EBIC,
      SABIC = fa_decomp$SABIC, m = m
    )})


## ------------------------------------------------------------------------
inform_crit %>% 
  gather(Criteria, value, -m) %>% 
  ggplot(aes(m, value, colour = Criteria)) + 
  geom_point() + 
  geom_line() + 
  theme_minimal() +
  theme(legend.position = 'top')


## ------------------------------------------------------------------------
inform_crit %>% 
  gather(Criteria, value, -m) %>% 
  group_by(Criteria) %>% 
  filter(value == min(value)) %>% 
  select(-value) %>% 
  spread(Criteria, m)


## ------------------------------------------------------------------------
m <- 2
Lhat <- decomp$loadings[,seq_len(m)] %*% 
  diag(decomp$sdev[seq_len(m)])
colnames(Lhat) <- paste0("PC", 1:2)
Psi_hat <- diag(cov(data) - tcrossprod(Lhat))

Lhat <- varimax(Lhat)$loadings
scores_pc <- scale(data, center = TRUE,
                    scale = TRUE) %*% 
  solve(Rn) %*% Lhat


## ------------------------------------------------------------------------
# MLE
fa_decomp <- factanal(data, factors = m)
Psi_mle <- fa_decomp$uniquenesses
Lmle <- fa_decomp$loadings

scores_mle <- scale(data, center = TRUE,
                    scale = TRUE) %*% 
  solve(Rn) %*% Lmle


## ------------------------------------------------------------------------
# Compare both loadings again
round(cor(scores_pc, scores_mle), 2)


## ------------------------------------------------------------------------
levelplot(unclass(Lhat),
          xlab = "", ylab = "")


## ------------------------------------------------------------------------
levelplot(unclass(Lmle),
          xlab = "", ylab = "")


## ------------------------------------------------------------------------
plot(scores_mle, type = 'n')
text(scores_mle, labels = rownames(scores_mle))
abline(h = 0, v = 0, lty = 2)


## ------------------------------------------------------------------------
# Let's plot the factors on the map
library(maps)
states <- map_data("state")

data_plot <- scores_mle %>%
  as.data.frame() %>% 
  rownames_to_column("region") %>% 
  mutate(region = tolower(region)) %>% 
  inner_join(states, by = "region")


## ------------------------------------------------------------------------
ggplot(data = data_plot) + 
  geom_polygon(aes(x = long, y = lat, 
                   fill = Factor1, 
                   group = group)) +
  coord_fixed(1.3) +
  ggthemes::theme_map() + 
  ggtitle("First Factor")


## ------------------------------------------------------------------------
ggplot(data = data_plot) + 
  geom_polygon(aes(x = long, y = lat, 
                   fill = Factor2, 
                   group = group)) +
  coord_fixed(1.3) +
  ggthemes::theme_map() + 
  ggtitle("Second Factor")

