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


# Weather data from NEON --------------------------------------------------

# Download steps take an hour or more total

# SAAT: Single aspirated air temperature
neon_download(product = "DP1.00002.001",
              site = target_sites)

# See what was downloaded
neon_index()


# Combine and export
walk(.x = target_sites,
     .f = ~ {
       
       temp_df <- neon_read(table = "SAAT_30min-basic", site = .x)
       
       write_rds(x = temp_df,
                 file = paste0("data/high_temporal_res/", tolower(.x), "_30min_temps.rds"))
       
     })

# Take a look at a single dataset for insight:
# There are 3 vertical positions, with their own vals for each 30 min
read_rds("data/high_temporal_res/blan_30min_temps.rds") %>%
  count(year(startDateTime), yday(startDateTime), verticalPosition)

# For now I'm going to stick with the lowest one, 010

# Compile all temperature datasets into one data frame
temp_compiled <- map_df(.x = list.files(path = "data/high_temporal_res/",
                                        pattern = "30min_temps.rds",
                                        full.names = TRUE),
                        .f = ~ {
                          
                          # Keep the site ID for later use
                          site_id <- gsub(pattern = "data/high_temporal_res/|_30min_temps.rds",
                                          replacement = "",
                                          x = .x)
                          
                          # Read in data and filter for the lowest tier of the tower
                          temp_data <- read_rds(.x) %>%
                            filter(verticalPosition == "010")
                          
                          # Get proportion of each day's mean records that are NA values
                          prop_na <- temp_data %>%
                            add_count(year = year(startDateTime), day = yday(startDateTime)) %>%
                            group_by(year, day) %>%
                            summarize(total_n = unique(n),
                                      non_na_count = sum(!is.na(tempSingleMean)),
                                      prop_na = round(((total_n - non_na_count) / total_n),
                                                      digits = 2)) %>%
                            select(year, day, total_n, prop_na)
                          
                          # Aggregate data to day level
                          agg_data <- temp_data %>%
                            group_by(year = year(startDateTime), day = yday(startDateTime)) %>%
                            summarize(mean_temp_c = mean(tempSingleMean, na.rm = TRUE),
                                      min_temp_c = min(tempSingleMinimum, na.rm = TRUE),
                                      max_temp_c = max(tempSingleMaximum, na.rm = TRUE),
                                      mean_var_temp = mean(tempSingleVariance, na.rm = TRUE)) %>%
                            ungroup()
                          
                          # Join aggregated data with NA proportions
                          data_output <- left_join(
                            x = agg_data,
                            y = prop_na,
                            by = c("year", "day")) %>%
                            mutate(
                              # Tag with site ID
                              site = site_id,
                              # Replace Infinite values and NaN with NAs
                              across(.cols = where(is.numeric),
                                     .fns =  ~na_if(., Inf)),
                              across(.cols = where(is.numeric),
                                     .fns =  ~na_if(., -Inf)),
                              across(.cols = where(is.double),
                                     .fns = ~if_else(condition = is.nan(.),
                                                     true = NA_real_,
                                                     false = .)),
                              date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
                            # Rearrange for easier reading
                            select(site, year, day, date, total_n, prop_na, everything())
                          
                          # Return to the outer function
                          return(data_output)
                          
                        }) 

# Export the final temp dataset
write_csv(x = temp_compiled,
          file = "data/neon_temp_compiled.csv")


# Relative humidity
neon_download(product = "DP1.00098.001",
              site = target_sites)

# Combine and export
walk(.x = target_sites,
     .f = ~ {
       
       temp_df <- neon_read(table = "RH_30min-basic", site = .x)
       
       write_rds(x = temp_df,
                 file = paste0("data/high_temporal_res/", tolower(.x), "_30min_rh.rds"))
       
     })

# Take a look at a single dataset for insight:
# There are 2 vertical positions, with their own vals for each 30 min
read_rds("data/high_temporal_res/blan_30min_rh.rds") %>%
  count(year(startDateTime), yday(startDateTime), verticalPosition)

# For now I'm going to stick with the lowest one, 000

# Compile all RH datasets into one data frame
rh_compiled <- map_df(.x = list.files(path = "data/high_temporal_res/",
                                      pattern = "30min_rh.rds",
                                      full.names = TRUE),
                      .f = ~ {
                        
                        # Keep the site ID for later use
                        site_id <- gsub(pattern = "data/high_temporal_res/|_30min_rh.rds",
                                        replacement = "",
                                        x = .x)
                        
                        # Read in data and filter for the lowest tier of the tower
                        rh_data <- read_rds(.x) %>%
                          filter(verticalPosition == "000")
                        
                        # Get proportion of each day's mean records that are NA values
                        prop_na <- rh_data %>%
                          add_count(year = year(startDateTime), day = yday(startDateTime)) %>%
                          group_by(year, day) %>%
                          summarize(total_n = unique(n),
                                    non_na_count = sum(!is.na(RHMean)),
                                    prop_na = round(((total_n - non_na_count) / total_n),
                                                    digits = 2)) %>%
                          select(year, day, total_n, prop_na)
                        
                        # Aggregate data to day level
                        agg_data <- rh_data %>%
                          group_by(year = year(startDateTime), day = yday(startDateTime)) %>%
                          summarize(mean_rh_pct = mean(RHMean, na.rm = TRUE),
                                    min_rh_pct = min(RHMinimum, na.rm = TRUE),
                                    max_rh_pct = max(RHMaximum, na.rm = TRUE),
                                    mean_var_rh = mean(RHVariance, na.rm = TRUE)) %>%
                          ungroup()
                        
                        # Join aggregated data with NA proportions
                        data_output <- left_join(
                          x = agg_data,
                          y = prop_na,
                          by = c("year", "day")) %>%
                          mutate(
                            # Tag with site ID
                            site = site_id,
                            # Replace Infinite values and NaN with NAs
                            across(.cols = where(is.numeric),
                                   .fns =  ~na_if(., Inf)),
                            across(.cols = where(is.numeric),
                                   .fns =  ~na_if(., -Inf)),
                            across(.cols = where(is.double),
                                   .fns = ~if_else(condition = is.nan(.),
                                                   true = NA_real_,
                                                   false = .)),
                            date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
                          # Rearrange for easier reading
                          select(site, year, day, date, total_n, prop_na, everything())
                        
                        # Return to the outer function
                        return(data_output)
                        
                      }) 

# Export the final temp dataset
write_csv(x = rh_compiled,
          file = "data/neon_rh_compiled.csv")


# Precip
neon_download(product = "DP1.00006.001",
              site = target_sites)

# Combine and export
walk(.x = target_sites,
     .f = ~ {
       
       temp_df <- neon_read(table = "THRPRE_30min-basic", site = .x)
       
       write_rds(x = temp_df,
                 file = paste0("data/high_temporal_res/", tolower(.x), "_30min_thrpre.rds"))
       
     })

# Take a look at a single dataset for insight:
# There's only one vertical position
read_rds("data/high_temporal_res/blan_30min_thrpre.rds") %>%
  pull(verticalPosition) %>%
  unique()

# Compile all precip datasets into one data frame
precip_compiled <- map_df(.x = list.files(path = "data/high_temporal_res/",
                                          pattern = "30min_thrpre.rds",
                                          full.names = TRUE),
                          .f = ~ {
                            
                            # Keep the site ID for later use
                            site_id <- gsub(pattern = "data/high_temporal_res/|_30min_thrpre.rds",
                                            replacement = "",
                                            x = .x)
                            
                            # Read in data and filter for the lowest tier of the tower
                            precip_data <- read_rds(.x) %>%
                              filter(verticalPosition == "000")
                            
                            # Get proportion of each day's bulk records that are NA values
                            prop_na <- precip_data %>%
                              add_count(year = year(startDateTime), day = yday(startDateTime)) %>%
                              group_by(year, day) %>%
                              summarize(total_n = unique(n),
                                        non_na_count = sum(!is.na(TFPrecipBulk)),
                                        prop_na = round(((total_n - non_na_count) / total_n),
                                                        digits = 2)) %>%
                              select(year, day, total_n, prop_na)
                            
                            # Aggregate data to day level
                            agg_data <- precip_data %>%
                              group_by(year = year(startDateTime), day = yday(startDateTime)) %>%
                              summarize(sum_precip_mm= sum(TFPrecipBulk, na.rm = TRUE)) %>%
                              ungroup()
                            
                            # Join aggregated data with NA proportions
                            data_output <- left_join(
                              x = agg_data,
                              y = prop_na,
                              by = c("year", "day")) %>%
                              mutate(
                                # Tag with site ID
                                site = site_id,
                                # Replace Infinite values and NaN with NAs
                                across(.cols = where(is.numeric),
                                       .fns =  ~na_if(., Inf)),
                                across(.cols = where(is.numeric),
                                       .fns =  ~na_if(., -Inf)),
                                across(.cols = where(is.double),
                                       .fns = ~if_else(condition = is.nan(.),
                                                       true = NA_real_,
                                                       false = .)),
                                ,
                                date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
                              # Rearrange for easier reading
                              select(site, year, day, date, total_n, prop_na, everything())
                            
                            # Return to the outer function
                            return(data_output)
                            
                          }) 

# Export the final temp dataset
write_csv(x = precip_compiled,
          file = "data/neon_precip_compiled.csv")


# Combine the weather datasets with the tick target data

tick_targets
temp_compiled
rh_compiled
precip_compiled

# Vector of final columns that will need some values turned into NAs again
cols_to_fix <- c("mean_temp", "min_temp", "max_temp", "mean_var_temp",
                 "mean_rh_pct", "min_rh_pct", "max_rh_pct", "mean_var_rh",
                 "sum_precip_mm")

full_dataset <- reduce(.x = list(
  tick_targets,
  temp_compiled %>%
    rename(prop_na_temp = prop_na) %>%
    group_by(site_id = site, year, mmwr_week = epiweek(date)) %>%
    summarize(site_id = toupper(site_id),
              mean_temp = mean(mean_temp_c, na.rm = TRUE),
              min_temp = min(min_temp_c, na.rm = TRUE),
              max_temp = max(max_temp_c, na.rm = TRUE),
              mean_var_temp = mean(mean_var_temp, na.rm = TRUE)) %>%
    ungroup(),
  rh_compiled %>%
    rename(prop_na_rh = prop_na) %>%
    group_by(site_id = site, year, mmwr_week = epiweek(date)) %>%
    summarize(site_id = toupper(site_id),
              mean_rh_pct = mean(mean_rh_pct, na.rm = TRUE),
              min_rh_pct = min(min_rh_pct, na.rm = TRUE),
              max_rh_pct = max(max_rh_pct, na.rm = TRUE),
              mean_var_rh = mean(mean_var_rh, na.rm = TRUE)),
  precip_compiled %>%
    rename(prop_na_precip = prop_na) %>%
    group_by(site_id = site, year, mmwr_week = epiweek(date)) %>%
    summarize(site_id = toupper(site_id),
              sum_precip_mm = sum(sum_precip_mm, na.rm = TRUE))),
  .f = left_join,
  by = c("site_id", "year", "mmwr_week")) %>%
  mutate(# Replace Infinite values and NaN with NAs
    across(.cols = cols_to_fix,
           .fns =  ~na_if(., Inf)),
    across(.cols = cols_to_fix,
           .fns =  ~na_if(., -Inf)),
    across(.cols = cols_to_fix,
           .fns = ~if_else(condition = is.nan(.),
                           true = NA_real_,
                           false = .))) %>%
  distinct()

write_csv(x = full_dataset,
          file = "data/ticks_with_weather.csv")




