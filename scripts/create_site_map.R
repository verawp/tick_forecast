library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggrepel)


# Load data ---------------------------------------------------------------

# State polygons in sf format
state_shapes <- ne_states(country = "United States of America",
                          returnclass = "sf")

# Site metadata
neon_sites <- read_csv(file = "data/Ticks_NEON_Field_Site_Metadata_20210928.csv")


# Prep and plot -----------------------------------------------------------

neon_sites_sf <- neon_sites %>%
  select(field_domain_id, field_site_id, field_site_name, field_site_state,
         field_latitude, field_longitude) %>%
  st_as_sf(x = .,
           coords = c("field_longitude", "field_latitude"),
           # WGS84
           crs = 4326)

# Large scale map
site_map <- ggplot() +
  geom_sf(data = state_shapes) +
  geom_sf(data = neon_sites_sf, color = "red") +
  geom_label_repel(
    data = neon_sites_sf,
    aes(label = field_site_id, geometry = geometry),
    stat = "sf_coordinates",
    min.segment.length = 0,
    nudge_y = 2,
    nudge_x = 2,
    color = "black") +
  coord_sf(xlim = c(-105, -70), ylim = c(28, 48)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

# Export map
ggsave(filename = "figures/site_overview_map.png", plot = site_map,
       device = "png", width = 6, height = 4, units = "in")
