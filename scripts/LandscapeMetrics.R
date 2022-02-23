



##################################################################
##################      Landscape metrics for ticks


library(tidyverse)
library(sf)
library(terra)
library(dplyr)
library(spData)
library(spDataLarge)


# Get lat longs

Sites <- read.csv("Ticks_NEON_Field_Site_Metadata_20210928.txt")
coords <- Sites[,c("field_longitude", "field_latitude")]
colnames(coords) <- c("X", "Y")
SiteNames<-Sites$field_site_name
SiteCodes<-Sites$field_site_id


# Project and supply spatial metadata

NEONsites = coords %>% 
  st_as_sf(coords = 1:2, crs = "EPSG:4326")

st_is_longlat(NEONsites)
st_crs(NEONsites)  # get CRS
st_crs(NEONsites)$IsGeographic #TRUE <- unprojected

target_crs <- "ESRI:102003"

# These coordinates are in the USA Contiguous Albers Equal Area Conic projected CRS 
# and the EPSG code is 102003. It's important to use an equal area projection for area relevant calculations.

SitesProj = sf::st_transform(NEONsites, target_crs)
st_crs(SitesProj)  # get CRS
st_crs(SitesProj)$IsGeographic #FALSE <- now it's projected

st_distance(SitesProj$geometry[1], SitesProj$geometry[2]) #note the projected coordinate system is in meters


# Make Buffers

SiteBuffer_s2 = st_buffer(NEONsites, dist = 1e4) # silent use of s2, like "projecting on the fly"
SiteBuffer_projected = st_buffer(SitesProj, 1e4) # explicit projected 10km buffer, better option


# Get croplands data
LandSectors <- as.list(rep(1,1,nrow(SiteBuffer_projected$geometry)))

for (i in 1:length(SiteBuffer_projected$geometry)) {
data <- GetCDLData(aoi = SiteBuffer_projected$geometry[i], year = 2018, type = 'b')
LandSectors[[i]] = data
raster::plot(LandSectors[[i]])
}


# Reclassify raster

library(raster)
library(landscapemetrics)

plot(LandSectors[[1]])
unique(LandSectors[[1]])

SiteCodes[1]
BLAN = LandSectors[[1]]

check_landscape(BLAN)

CDL<-read.csv("CDLcodes.csv")
RclTable<-CDL["Codes", "NewValues"]

BLANsimp <- reclassify(BLAN, as.matrix(RclTable))
plot(BLANsimp)


# Visualization functions

BLANsimp  %>%
show_cores() #show cores


# Check options for metrics

list_lsm(level = "patch") #see options
list_lsm(level = "class")   %>%  #see options
  print(n=Inf)
list_lsm(level = "landscape")  %>%  #see options
print(n=Inf)


# Calculate a single metric
lsm_c_clumpy(BLANsimp)


# Calculate multiple metrics
metrics = calculate_lsm(BLANsimp, 
              what = c("lsm_c_pland", "lsm_l_ta", "lsm_l_te"))

