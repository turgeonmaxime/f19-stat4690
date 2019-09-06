################################
# Data Visualization in Base R #
#           STAT 4690          #
################################
library(dslabs)

# Univariate Plots----
# Histograms
hist(olive$oleic)

# Density plots
plot(density(olive$oleic))

# ECDFs
plot(ecdf(olive$oleic))

# Boxplots
boxplot(olive$oleic)

# Bivariate Plots----
# Scatter plots
plot(stars$magnitude, 
     stars$temp)

# 2D Density plots
library(MASS)

image(kde2d(stars$magnitude, 
            stars$temp))

# Bagplots
library(mrfDepth)

stars_bg <- stars[,c("magnitude", "temp")]
bagplot(compBagplot(stars_bg))

# Different approach
# Be careful: both packages use the same function name
library(aplpack)

bagplot(stars$magnitude,
        stars$temp)

# Pairs plots----
olive_acids <- subset(olive, select = c(-region, -area))

plot(olive_acids)

## put histograms on the diagonal
## From the help page for graphics::pairs
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

plot(olive_acids, 
     diag.panel = panel.hist)
