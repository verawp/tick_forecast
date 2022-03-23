library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggrepel)
library(corrplot)


# Load data ---------------------------------------------------------------

# State polygons in sf format
state_shapes <- ne_states(country = "United States of America",
                          returnclass = "sf")

# Tick data
# https://projects.ecoforecast.org/neon4cast-docs/theme-tick-populations.html
# tick_data <- read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz",
#                       guess_max = 1e6)

# Tick data with weather attached
tick_data <- read_csv(file = "data/ticks_with_filled_weather.csv")

# Preview the data
head(tick_data)


# Variable exploration ----------------------------------------------------

correlations <- tick_data %>%
  select(amblyomma_americanum, mean_temp, min_temp, max_temp, mean_var_temp,
         mean_rh_pct, min_rh_pct, max_rh_pct) %>%
  cor(use = "pairwise.complete.obs")

corrplot(corr = correlations, type = "upper")































