## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ------------------------------------------------------------------------
Delta <- dist(swiss)
D <- Delta^2

# Center columns
B <- scale(D, center = TRUE, scale = FALSE)
# Center rows
B <- t(scale(t(B), center = TRUE, scale = FALSE))
B <- -0.5 * B


## ------------------------------------------------------------------------
decomp <- eigen(B)
Lambda <- diag(pmax(decomp$values, 0))
X_tilde <-  decomp$vectors %*% Lambda^0.5


## ------------------------------------------------------------------------
plot(X_tilde)


## ---- message = FALSE----------------------------------------------------
mds <- cmdscale(Delta, k = 2)

all.equal(X_tilde[,1:2], mds,
          check.attributes = FALSE)


## ---- message = FALSE----------------------------------------------------
library(tidyverse)
# Let's add annotations
dimnames(X_tilde) <- list(rownames(swiss), 
                          paste0("MDS", seq_len(ncol(X_tilde))))
X_tilde <- as.data.frame(X_tilde) %>% 
  rownames_to_column("District")


## ------------------------------------------------------------------------
X_tilde <- X_tilde %>% 
  mutate(Canton = case_when(
    District %in% c("Courtelary", "Moutier", 
                    "Neuveville") ~ "Bern",
    District %in% c("Broye", "Glane", "Gruyere",
                    "Sarine", "Veveyse") ~ "Fribourg",
    District %in% c("Conthey", "Entremont", "Herens", 
                    "Martigwy", "Monthey", 
                    "St Maurice", "Sierre", 
                    "Sion") ~ "Valais"))


## ------------------------------------------------------------------------
X_tilde <- X_tilde %>% 
  mutate(Canton = case_when(!is.na(Canton) ~ Canton,
    District %in% c("Boudry", "La Chauxdfnd", 
                    "Le Locle", "Neuchatel", 
                    "ValdeTravers", 
                    "Val de Ruz") ~ "Neuchatel",
    District %in% c("V. De Geneve", "Rive Droite", 
                    "Rive Gauche") ~ "Geneva",
    District %in% c("Delemont", "Franches-Mnt", 
                    "Porrentruy") ~ "Jura",
    TRUE ~ "Vaud"))


## ---- message = FALSE----------------------------------------------------
library(ggrepel)
X_tilde %>% 
  ggplot(aes(MDS1, MDS2)) + 
  geom_point(aes(colour = Canton)) + 
  geom_label_repel(aes(label = District)) +
  theme_minimal() +
  theme(legend.position = "top")


## ------------------------------------------------------------------------
library(psych)
cities[1:5, 1:5]


## ------------------------------------------------------------------------
mds <- cmdscale(cities, k = 2)
colnames(mds) <- c("MDS1", "MDS2")

mds <- mds %>% 
  as.data.frame %>% 
  rownames_to_column("Cities")


## ------------------------------------------------------------------------
mds %>% 
  ggplot(aes(MDS1, MDS2)) + 
  geom_point() +
  geom_label_repel(aes(label = Cities)) +
  theme_minimal() 


## ------------------------------------------------------------------------
mds %>% 
  mutate(MDS1 = -MDS1, MDS2 = -MDS2) %>% 
  ggplot(aes(MDS1, MDS2)) + 
  geom_point() +
  geom_label_repel(aes(label = Cities)) +
  theme_minimal() 


## ---- message=FALSE------------------------------------------------------
library(MASS)

Delta <- dist(swiss)
mds <- sammon(Delta, k = 2)

plot(mds$points)


## ------------------------------------------------------------------------
# Fit for different values of k
stresses <- sapply(seq(2, 10),
                   function(k) {
                     sammon(Delta, k = k,
                            trace = FALSE)$stress
                     })
plot(seq(2, 10), stresses, type = 'b')


## ------------------------------------------------------------------------
library(scatterplot3d)
mds <- sammon(Delta, k = 3)
scatterplot3d(mds$points, 
              xlab = "MDS1", ylab = "MDS2", 
              zlab = "MDS3")


## ---- message = FALSE----------------------------------------------------
mds_s <- sammon(Delta, k = 2)
mds_k <- isoMDS(Delta, k = 2)

par(mfrow = c(1, 2))
plot(mds_s$points, main = "Sammon",
     xlab = "MDS1", ylab = "MDS2")
plot(mds_k$points, main = "Kruskal",
     xlab = "MDS1", ylab = "MDS2")


## ------------------------------------------------------------------------
# Sammon and Kruskal have different
# optimal k
stresses <- sapply(seq(2, 10),
                   function(k) {
                     isoMDS(Delta, k = k,
                            trace = FALSE)$stress
                     })
plot(seq(2, 10), stresses, type = 'b')


## ------------------------------------------------------------------------
mds_opt_s <- sammon(Delta, k = 3,
                    trace = FALSE)
mds_opt_k <- isoMDS(Delta, k = 6,
                    trace = FALSE)

# Let's cluster in the MDS space
cluster_s <- kmeans(mds_opt_s$points, centers = 2)
cluster_k <- kmeans(mds_opt_k$points, centers = 2)


## ------------------------------------------------------------------------
par(mfrow = c(1, 2))
plot(mds_s$points, main = "Sammon",
     xlab = "MDS1", ylab = "MDS2",
     col = cluster_s$cluster)
plot(mds_k$points, main = "Kruskal",
     xlab = "MDS1", ylab = "MDS2",
     col = cluster_k$cluster)


## ------------------------------------------------------------------------
# More interestingly, you can use MDS to
# cluster data where you only have distances
stresses <- sapply(seq(2, 6),
                   function(k) {
                     isoMDS(as.matrix(cities), 
                            k = k,
                            trace = FALSE)$stress
                     })
plot(seq(2, 6), stresses, type = 'b')


## ------------------------------------------------------------------------
mds_cities <- isoMDS(as.matrix(cities), k = 6,
                     trace = FALSE)
cluster_cities <- kmeans(mds_cities$points, 
                         centers = 2)

plot(mds_cities$points, main = "Kruskal",
     xlab = "MDS1", ylab = "MDS2",
     type = 'n')
text(mds_cities$points, colnames(cities),
     col = cluster_cities$cluster)

