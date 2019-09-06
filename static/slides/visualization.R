## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ---- message = FALSE, eval = FALSE--------------------------------------
## library(tidyverse)
## 
## count(mtcars, cyl)
## # Or with the pipe
## mtcars %>% count(cyl)


## ----eval = FALSE--------------------------------------------------------
## mutate(mtcars, liters_per_100km = mpg/235.215)


## ----eval = FALSE--------------------------------------------------------
## filter(mtcars, cyl %in% c(6, 8))


## ----eval = FALSE--------------------------------------------------------
## mtcars %>% group_by(cyl) %>%
##   summarise(avg_mpg = mean(mpg))


## ---- message=FALSE------------------------------------------------------
library(tidyverse)
library(dslabs)

dim(olive)

olive %>% 
  ggplot(aes(oleic)) + 
  geom_histogram()

olive %>% 
  ggplot(aes(oleic, fill = region)) + 
  geom_histogram() + 
  theme(legend.position = 'top')

# Or with facets
olive_bg <- olive %>% dplyr::select(-region)
olive %>% 
  ggplot(aes(oleic, fill = region)) + 
  geom_histogram(data = olive_bg, 
                 fill = 'grey') +
  geom_histogram() +
  facet_grid(. ~ region) +
  theme(legend.position = 'top')


## ---- message=FALSE------------------------------------------------------
olive %>% 
  ggplot(aes(oleic)) + 
  geom_density()

olive %>% 
  ggplot(aes(oleic, fill = region)) + 
  geom_density(alpha = 0.5) + 
  theme(legend.position = 'top')


## ---- message=FALSE------------------------------------------------------
olive %>% 
  ggplot(aes(oleic)) + 
  stat_ecdf() + 
  ylab("Cumulative Probability")

# You can add a "rug"
olive %>% 
  ggplot(aes(oleic)) + 
  stat_ecdf() + 
  geom_rug(sides = "b") + 
  ylab("Cumulative Probability")

olive %>% 
  ggplot(aes(oleic, colour = region)) + 
  stat_ecdf() + 
  ylab("Cumulative Probability") + 
  theme(legend.position = 'top')


## ---- message=FALSE------------------------------------------------------
olive %>% 
  ggplot(aes(y = oleic)) + 
  geom_boxplot(x = 0)

olive %>% 
  ggplot(aes(x = region, y = oleic)) + 
  geom_boxplot()

# Add all points on top of boxplots
# Note: need to remove outliers or you will get 
#       duplicates
olive %>% 
  ggplot(aes(x = region, y = oleic)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.25, height = 0)


## ------------------------------------------------------------------------
stars %>% 
  ggplot(aes(magnitude, temp)) + 
  geom_point()

stars %>% 
  ggplot(aes(magnitude, temp)) + 
  geom_point(aes(colour = type))


## ------------------------------------------------------------------------
library(scatterplot3d)

greenhouse_gases %>% 
  spread(gas, concentration) %>% 
  with(scatterplot3d(CH4,   # x axis
                     CO2,   # y axis
                     N2O    # z axis
))


## ------------------------------------------------------------------------
stars %>% 
  ggplot(aes(magnitude, temp)) + 
  geom_point(aes(colour = type)) +
  geom_density_2d()


## ---- cache = TRUE, warning = FALSE, message = FALSE---------------------
devtools::source_gist("00772ccea2dd0b0f1745", 
                      filename = "000_geom_bag.r")
devtools::source_gist("00772ccea2dd0b0f1745", 
                      filename = "001_bag_functions.r")

stars %>% 
  ggplot(aes(magnitude, temp)) + 
  geom_bag() + 
  theme_bw()

stars %>% 
  ggplot(aes(magnitude, temp)) + 
  geom_bag() + 
  geom_point(aes(colour = type)) +
  theme_bw()

gapminder %>%
  filter(year == 2012,
         !is.na(infant_mortality)) %>% 
  ggplot(aes(infant_mortality, life_expectancy)) + 
  geom_bag(aes(fill = continent)) + 
  geom_point(aes(colour = continent)) +
  theme_bw()

gapminder %>%
  filter(year == 2012,
         !is.na(infant_mortality)) %>% 
  ggplot(aes(infant_mortality, life_expectancy)) + 
  geom_bag(aes(fill = continent)) + 
  geom_point(aes(colour = continent)) +
  facet_wrap(~continent) +
  theme_bw()


## ---- message = FALSE----------------------------------------------------
library(GGally)

olive %>% 
  dplyr::select(-region, -area) %>% 
  ggpairs


## ---- cache = TRUE-------------------------------------------------------
library(ggforce)

olive %>% 
  dplyr::select(-region, -area) %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) + 
  geom_point() + 
  facet_matrix(vars(everything()))

olive %>% 
  dplyr::select(-region, -area) %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) + 
  geom_point() + 
  geom_autodensity() +
  facet_matrix(vars(everything()),
               layer.diag = 2)

olive %>% 
  dplyr::select(-region, -area) %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) + 
  geom_point() + 
  geom_autodensity() +
  geom_density2d() +
  facet_matrix(vars(everything()), 
               layer.diag = 2, 
               layer.upper = 3)

