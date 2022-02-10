# Download and clean tick data: Adapted from https://github.com/eco4cast/neon4cast-ticks/blob/master/02_ticks_targets.R
# on 2022-02-10 by MRB

# Check for neonstore package and install if needed
if(!"neonstore" %in% installed.packages()){
  remotes::install_github("cboettig/neonstore", ref = "patch/api-updates")
}

library(tidyverse) # for data wrangling and piping (dplyr probably ok)
library(lubridate) # for finding year from dates
library(stringr) # for searching within character strings 
library(parallel) # for using more than one core in download
library(MMWRweek) # for converting from date to MMWR week
library(neonstore)

# Select target species and life stage
target_species <- "Amblyomma americanum"
target_lifestage <- "Nymph"

target_sites <- c("BLAN", "KONZ", "LENO", "ORNL", "OSBS", "SCBI",
                  "SERC", "TALL", "UKFS")


Sys.setenv("NEONSTORE_HOME" = "/neonstore")
Sys.setenv("NEONSTORE_DB" = "/neonstore")

# Get data from NEON
product <- "DP1.10093.001"
neon_download(product = product,
              site = target_sites)

tick_field_raw <- neon_read("tck_fielddata-basic", keep_filename = TRUE)
tick_taxon_raw <- neon_read("tck_taxonomyProcessed-basic")

# There are lots of reasons why sampling didn't occur (logistics, too wet, too cold, etc.)
# so, keep records when sampling occurred
tick_field <- tick_field_raw %>% 
  filter(totalSampledArea > 0) %>% 
  mutate(time = floor_date(collectDate, unit = "day")) %>% 
  unite(namedLocation, time, col = "occasionID", sep = "_")

# Combine adults into single category and make wide to get zero counts
tick_taxon_wide <- tick_taxon_raw %>% 
  filter(sampleCondition == "OK") %>% # remove taxonomy samples with quality issues
  mutate(sexOrAge = if_else(sexOrAge == "Female" | sexOrAge == "Male", 
                            "Adult",     # convert to Adult
                            sexOrAge),
         time = floor_date(collectDate, unit = "day")) %>% 
  unite(namedLocation, time, col = "occasionID", sep = "_") %>% 
  pivot_wider(id_cols = occasionID, # make wide by species and life stage
              names_from = c(acceptedTaxonID, sexOrAge),
              values_from = individualCount, 
              names_sep = "_",
              # duplicates occur because of Adults that where F/M - add them 
              values_fn = {sum}, 
              values_fill = 0)

# Join taxonomy and field data
tick_joined <- left_join(tick_taxon_wide, tick_field, by = "occasionID") %>% 
  select(-NA_NA, -geodeticDatum, -samplingImpractical, -targetTaxaPresent,
         -adultCount, -nymphCount, -larvaCount, -samplingProtocolVersion, 
         -measuredBy, -sampleCode, -biophysicalCriteria, -plotType)

# All the species column names
spp_cols <- tick_joined %>% 
  select(contains("Larva"), contains("Nymph"), contains("Adult")) %>% 
  colnames()

# Get matching taxon ids
taxon_ids <- tick_taxon_raw %>%
  filter(!is.na(acceptedTaxonID)) %>% 
  select(acceptedTaxonID, scientificName, taxonRank) %>% 
  distinct() 

# Pivot to long version
tick_long <- tick_joined %>% 
  pivot_longer(cols = all_of(spp_cols), 
               names_to = "taxonAge",
               values_to = "processedCount",
               values_drop_na = TRUE) %>% 
  separate(col = taxonAge, into = c("acceptedTaxonID", "lifeStage"), sep = "_")

# Add taxon ids
tick_long <- left_join(tick_long, taxon_ids, by = "acceptedTaxonID") 

# Standardize the data and subset to targets
tick_standard <- tick_long %>% 
  filter(siteID %in% target_sites, # sites we want
         lifeStage == target_lifestage, # life stage we want
         scientificName == target_species, # species we want
         grepl("Forest", nlcdClass)) %>%  # forest plots
  mutate(date = floor_date(collectDate, unit = "day"),
         date = ymd(date),
         year = year(date),
         mmwrWeek = MMWRweek(date)$MMWRweek,
         time = MMWRweek2Date(year, mmwrWeek)) %>% 
  select(time, processedCount, totalSampledArea, siteID) %>%
  mutate(totalSampledArea = as.numeric(totalSampledArea)) %>% 
  group_by(siteID, time) %>%
  summarise(totalCount = sum(processedCount), # all counts in a week
            totalArea = sum(totalSampledArea),# total area surveyed in a week
            amblyomma_americanum = totalCount / totalArea * 1600) %>% # scale to the size of a plot
  # Add MMWR week (CDC's Morbidity and Mortality Weekly Report)
  mutate(mmwrWeek = MMWRweek(time)$MMWRweek) %>% 
  arrange(siteID, time) %>% 
  filter() %>% 
  select(time, mmwrWeek, siteID, amblyomma_americanum)

# In case NEON makes a provisional data release during the challenge
# need to make sure that we filter out 2021 data that is in the "future"
run_date <- today()
year(run_date) <- year(run_date) - 1
run_mmwrWeek <- MMWRweek(run_date)$MMWRweek
challenge_time <- MMWRweek2Date(year(run_date), run_mmwrWeek)

tick_targets <- tick_standard %>% 
  filter(time < challenge_time) %>%
  mutate(year = year(time),
         day_num = yday(time)) %>%
  select(site_id = siteID, mmwr_week = mmwrWeek, time, year, day_num, amblyomma_americanum)

# Write targets to csv
write_csv(tick_targets,
          file = "data/ticks_target.csv")
