## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache = FALSE, message = FALSE)


## -----------------------------------------------------------------------------
(Lmat <- matrix(c(0, 0, 1, 0, 1, 0, 0, 0,
                  1, 1, 0, 1, 0, 0, 0, 0), ncol = 4,
                byrow = TRUE))


## -----------------------------------------------------------------------------
# Alternatively
library(igraph)
edge_list <- c(c("A", "B"), c("A", "C"), c("B", "C"),
               c("C", "A"), c("D", "C"))
graph <- graph_from_edgelist(
  matrix(edge_list, ncol = 2,
         byrow = TRUE)
  )
plot(graph)


## -----------------------------------------------------------------------------
m_vect <- colSums(Lmat)
# Initialize
pr_vect <- rep(0.25, 4)
for (i in seq_len(4)) {
  pr_vect[i] <- sum(Lmat[i,]*pr_vect/m_vect)
}
pr_vect


## -----------------------------------------------------------------------------
# Repeat 5 more times
for(count in seq_len(5)) {
  for (i in seq_len(4)) {
    pr_vect[i] <- sum(Lmat[i,]*pr_vect/m_vect)
    }
}
pr_vect


## -----------------------------------------------------------------------------
Amat <- Lmat %*% diag(1/m_vect)
decomp <- eigen(Amat)
decomp$values


## -----------------------------------------------------------------------------
Re(decomp$vectors[,1])

# Which is proportional to what we found
pr_vect
Re(decomp$vectors[,1])/-(4/3)


## -----------------------------------------------------------------------------
d <- 0.85
n <- length(m_vect)
Amat <- (1 - d) * matrix(1, ncol = n, nrow = n)/n + 
  d * Lmat %*% diag(1/m_vect)
decomp <- eigen(Amat)
decomp$values

Re(decomp$vectors[,1])


## -----------------------------------------------------------------------------
page_rank(graph, damping = 0.85)$vector
Re(decomp$vectors[,1])/(sum(Re(decomp$vectors[,1])))


## ---- message=FALSE-----------------------------------------------------------
library(tidyverse)
url <- paste0("https://raw.githubusercontent.com/",
        "turgeonmaxime/twitter_pagerank_STAT4690/", 
              "master/edge_list_subset.csv")
edge_sub <- read_csv(url) %>% 
  mutate(Node_Id_1 = as.character(Node_Id_1),
         Node_Id_2 = as.character(Node_Id_2)) %>% 
  as.matrix()


## -----------------------------------------------------------------------------
graph_twitter <- graph_from_edgelist(edge_sub)

plot(graph_twitter, vertex.label = NA, 
     vertex.size = 2, edge.arrow.size = 0.5)


## -----------------------------------------------------------------------------
# PageRank algorithm
PR_vect <- page_rank(graph_twitter)$vector

# Visualize scores
library(colorspace)
quartiles <- quantile(PR_vect,
                      probs = seq(0, 1, length.out=5))
colours <- cut(PR_vect, breaks = quartiles, 
               labels = sequential_hcl(4), 
               include.lowest = TRUE)


## -----------------------------------------------------------------------------
plot(graph_twitter, vertex.label = NA, 
     vertex.size = 2, edge.arrow.size = 0.5,
     vertex.color = as.character(colours))

