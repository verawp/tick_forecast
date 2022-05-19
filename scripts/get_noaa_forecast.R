library(tidyverse)
library(lubridate)
library(neon4cast)

tick_sites <- read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz") %>% 
 pull(siteID) %>%
  unique()

download_noaa(siteID = tick_sites, dir = "data/weather_forecasts/")
noaa_fc <- stack_noaa(dir = "data/weather_forecasts/",
                      model = "NOAAGEFS_6hr")
noaa_fc
