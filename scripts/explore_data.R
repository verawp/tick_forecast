library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggrepel)
library(riem)


# Load data ---------------------------------------------------------------

# State polygons in sf format
state_shapes <- ne_states(country = "United States of America",
                          returnclass = "sf")

# Tick data
# https://projects.ecoforecast.org/neon4cast-docs/theme-tick-populations.html
tick_data <- read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz",
                      guess_max = 1e6)

# Preview the data
head(tick_data)

# Convert tick data to spatial format
tick_sf <- st_as_sf(tick_data,
                    coords = c("decimalLongitude", "decimalLatitude")) %>%
  st_set_crs(st_crs(state_shapes))


# Create reference map of sites -------------------------------------------

# Manually summarize the locations to a single point each
plot_counts <- tribble(
  ~label, ~lat, ~long,
  "BLAN: 2 plots", 39.086224, -77.972736,
  "KONZ: 1 plot", 39.103714, -96.597835, 
  "ORNL: 6 plots", 35.983081, -84.211655,
  "SCBI: 2 plots", 38.893863, -78.1427,
  "SERC: 5 plots", 38.868184, -76.535301,
  "TALL: 3 plots", 32.976366, -87.432192,
  "UKFS: 3 plots", 39.037238, -95.19863
) %>%
  st_as_sf(coords = c("long", "lat")) %>%
  st_set_crs(st_crs(state_shapes))

# Create map
site_map <- ggplot() +
  geom_sf(data = state_shapes) +
  geom_sf(data = plot_counts, color = "black", fill = "red", size = 2, pch = 21) +
  geom_label_repel(
    data = plot_counts,
    aes(label = label, geometry = geometry),
    stat = "sf_coordinates",
    min.segment.length = 0,
    nudge_y = 2,
    nudge_x = 2,
    color = "black") +
  coord_sf(xlim = c(-100, -65), ylim = c(30, 50)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

ggsave(filename = "figures/site_overview_map.png", plot = site_map,
       device = "png", width = 6, height = 4, units = "in")


# Explore some options for daily weather data -----------------------------

# Locations in KS
tick_data %>%
  select(siteID, decimalLatitude, decimalLongitude) %>%
  distinct() %>%
  filter(siteID %in% c("KONZ", "UKFS"))

# KONZ, distance: ~5.4 mi
riem_measures(station = "MHK", date_start = "2015-01-01", date_end = "2020-01-01")

# UKFS, distance: ~0 mi
riem_measures(station = "LWC", date_start = "2015-01-01", date_end = "2020-01-01")


# Locations in AL
tick_data %>%
  select(siteID, decimalLatitude, decimalLongitude) %>%
  distinct() %>%
  filter(siteID %in% c("TALL"))

# TALL, distance: ~21.5 mi
riem_measures(station = "TCL", date_start = "2015-01-01", date_end = "2020-01-01")


# Locations in TN
tick_data %>%
  select(siteID, decimalLatitude, decimalLongitude) %>%
  distinct() %>%
  filter(siteID %in% c("ORNL"))

# ORNL, distance: ~21.5 mi
riem_measures(station = "OQT", date_start = "2014-01-01", date_end = "2020-01-01")

# Locations in VA
tick_data %>%
  select(siteID, decimalLatitude, decimalLongitude) %>%
  distinct() %>%
  filter(siteID %in% c("BLAN", "SCBI"))

# SCBI, distance: ~10.85 mi, doesn't capture first couple months of data though
riem_measures(station = "FRR", date_start = "2014-01-01", date_end = "2020-01-01")
# SCBI / BLAN, distance: ~13.5 mi / ~5.3 mi
riem_measures(station = "OKV", date_start = "2014-01-01", date_end = "2020-01-01")


# Locations in MD
tick_data %>%
  select(siteID, decimalLatitude, decimalLongitude) %>%
  distinct() %>%
  filter(siteID %in% c("SERC"))

# SERC, distance: ~8.9 mi
riem_measures(station = "NAK", date_start = "2015-01-01", date_end = "2020-01-01")


