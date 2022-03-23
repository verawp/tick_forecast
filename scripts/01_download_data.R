# Download tick and weather data: Tick section adapted from https://github.com/eco4cast/neon4cast-ticks/blob/master/02_ticks_targets.R
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
library(riem)
library(measurements)
library(prism)
library(sf)
library(raster)


# 1. Download and clean tick data -----------------------------------------

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


# Get info on non-sampling dates that were dropped earlier:
# There are lots of reasons why sampling didn't occur (logistics, too wet, too cold, etc.)
# so, keep records when sampling occurred
tick_no_field <- tick_field_raw %>% 
  filter(totalSampledArea <= 0 | is.na(totalSampledArea)) %>% 
  mutate(time = floor_date(collectDate, unit = "day")) %>% 
  unite(namedLocation, time, col = "occasionID", sep = "_") %>%
  select(siteID, nlcdClass, decimalLatitude, decimalLongitude, elevation, totalSampledArea, remarks)

write_csv(x = tick_no_field,
          file = "data/ticks_missed_sample_dates.csv")


# 2. Download and process NEON  data --------------------------------------

# Note: download steps take an hour or more total


# 2.1 SAAT: Single aspirated air temperature ------------------------------

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


# 2.2 RH: Relative humidity -----------------------------------------------

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


# 2.3 PRIPRE: Primary precipitation ---------------------------------------

# Precip
neon_download(product = "DP1.00006.001",
              site = target_sites)

# Combine and export. Different sites may have different available precipitation
# types
walk(.x = target_sites,
     .f = ~ {
       
       # First try retrieving primary precip
       pri_df <- neon_read(table = "PRIPRE_30min-basic", site = .x)
       
       # If primary is NULL, try secondary
       if (is.null(pri_df)) {
         
         sec_df <- neon_read(table = "SECPRE_30min-basic", site = .x)
         
         # If secondary is NULL, try throughfall
         if (is.null(sec_df)) {
           
           # Keep throughfall data and mark it
           precip_df <- neon_read(table = "THRPRE_30min-basic", site = .x) %>%
             mutate(precip_type = "THRPRE")
           
           # If secondary is not NULL, keep it and mark it
         } else {
           
           precip_df <- sec_df %>%
             mutate(precip_type = "SECPRE")
           
         }
         # If primary is not NULL, keep it and mark it
       } else {
         
         precip_df <- pri_df %>%
           mutate(precip_type = "PRIPRE")
         
       }
       
       # Write the resulting data out to csv
       write_rds(x = precip_df,
                 file = paste0("data/high_temporal_res/", tolower(.x), "_30min_precip.rds"))
     })

# Take a look at a single dataset for insight:
read_rds("data/high_temporal_res/blan_30min_precip.rds")

# Compile all precip datasets into one data frame
precip_compiled <- map_df(.x = list.files(path = "data/high_temporal_res/",
                                          pattern = "30min_precip.rds",
                                          full.names = TRUE),
                          .f = ~ {
                            
                            # Keep the site ID for later use
                            site_id <- gsub(pattern = "data/high_temporal_res/|_30min_precip.rds",
                                            replacement = "",
                                            x = .x)
                            
                            # Read in data
                            precip_data <- read_rds(.x) %>%
                              # Different measurement column names based on site;
                              # rename so that they can be used agnostically
                              rename(precip_bulk = matches("PrecipBulk"))
                            
                            
                            # Get proportion of each day's bulk records that are NA values
                            prop_na <- precip_data %>%
                              add_count(year = year(startDateTime), day = yday(startDateTime)) %>%
                              group_by(year, day) %>%
                              summarize(total_n = unique(n),
                                        non_na_count = sum(!is.na(precip_bulk)),
                                        prop_na = round(((total_n - non_na_count) / total_n),
                                                        digits = 2)) %>%
                              select(year, day, total_n, prop_na)
                            
                            # Aggregate data to day level
                            agg_data <- precip_data %>%
                              group_by(year = year(startDateTime), day = yday(startDateTime)) %>%
                              summarize(sum_precip_mm = sum(precip_bulk, na.rm = TRUE),
                                        precip_type = unique(precip_type)) %>%
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


# 2.4 Small mammal trapping -----------------------------------------------

# Small mammal box trapping
neon_download(product = "DP1.10072.001",
              site = target_sites)

# Combine and export
walk(.x = target_sites,
     .f = ~ {
       temp_df <- neon_read(table = "mam_pertrapnight", site = .x)
       
       write_rds(x = temp_df,
                 file = paste0("data/high_temporal_res/", tolower(.x), "_small_mammal.rds"))
     })

# Take a look at a single dataset for insight:
read_rds("data/high_temporal_res/blan_small_mammal.rds")


# Compile all mammal datasets into one data frame
mammal_daily_compiled <- map_df(.x = list.files(path = "data/high_temporal_res/",
                                                pattern = "small_mammal.rds",
                                                full.names = TRUE),
                                .f = ~ {
                                  
                                  # Keep the site ID for later use
                                  site_id <- gsub(pattern = "data/high_temporal_res/|_small_mammal.rds",
                                                  replacement = "",
                                                  x = .x)
                                  
                                  # Read in data and select the most helpful cols
                                  small_mammal <- read_rds(.x) %>%
                                    select(siteID, plotID, trapCoordinate, plotType,
                                           nlcdClass, trapStatus, collectDate,
                                           lifeStage, taxonID, scientificName,
                                           taxonRank, contains("ticks"))
                                  
                                  mammal_success <- small_mammal %>%
                                    # Something captured
                                    filter(trapStatus %in% c("4 - more than 1 capture in one trap",
                                                             "5 - capture"))
                                  
                                  # Return to the outer function
                                  return(mammal_success)
                                  
                                }) 

# Export the final temp dataset
write_csv(x = mammal_daily_compiled,
          file = "data/neon_mammals_daily_compiled.csv")

mammal_larval_ticks <- mammal_daily_compiled %>%
  add_count(year = year(collectDate), mmwr_week = epiweek(collectDate), siteID,
            lifeStage, scientificName,
            larvalTicksAttached, name = "ind_w_larval_ticks") %>%
  filter(larvalTicksAttached == "Y") %>%
  select(year, mmwr_week, siteID, lifeStage, scientificName, ind_w_larval_ticks) %>%
  distinct()

mammal_nymphal_ticks <- mammal_daily_compiled %>%
  add_count(year = year(collectDate), mmwr_week = epiweek(collectDate), siteID,
            lifeStage, scientificName,
            nymphalTicksAttached, name = "ind_w_nymphal_ticks") %>%
  filter(nymphalTicksAttached == "Y") %>%
  select(year, mmwr_week, siteID, lifeStage, scientificName, ind_w_nymphal_ticks) %>%
  distinct()

mammal_adult_ticks <- mammal_daily_compiled %>%
  add_count(year = year(collectDate), mmwr_week = epiweek(collectDate), siteID,
            lifeStage, scientificName,
            adultTicksAttached, name = "ind_w_adult_ticks") %>%
  filter(adultTicksAttached == "Y") %>%
  select(year, mmwr_week, siteID, lifeStage, scientificName, ind_w_adult_ticks) %>%
  distinct()

mammal_ticks <- reduce(.x = list(mammal_larval_ticks, mammal_nymphal_ticks, mammal_adult_ticks),
                       .f = full_join,
                       by = c("year", "mmwr_week", "siteID", "lifeStage", "scientificName"))

mammal_counts <- mammal_daily_compiled %>%
  count(year = year(collectDate), mmwr_week = epiweek(collectDate),
        siteID, lifeStage, scientificName)

mammals_weekly <- full_join(x = mammal_counts,
                            y = mammal_ticks,
                            by = c("year", "mmwr_week", "siteID", "lifeStage", "scientificName")) %>%
  mutate(across(.cols = c(n, ind_w_larval_ticks, ind_w_nymphal_ticks, ind_w_adult_ticks),
                .fns = ~replace_na(data = ., replace = 0)))

# Export the final temp dataset
write_csv(x = mammals_weekly,
          file = "data/neon_mammals_weekly_compiled.csv")


# 3. External (non-NEON) weather data downloads ---------------------------

# A review of the NEON weather data shows that there's a fair amount of gaps.
# See script 02_clean_and_manage_data.R to see how many gaps there are in the
# combined NEON tick & weather output.
# Because of the gaps we'll want some external weather data sources, which I
# pull below

# Our sites:
neon_sites <- read_csv(file = "data/Ticks_NEON_Field_Site_Metadata_20210928.csv") %>%
  select(field_domain_id, field_site_id, field_site_name, field_site_state,
         field_latitude, field_longitude) 

# Note: I went ahead and looked up the start/end dates of samples for each of
# the sites in the dataset. I'm going to (for now) retrieve weather data
# as early as the September before sampling started so that we have some
# info covering the winter before sampling started

temp_compiled <- read_csv(file = "data/neon_temp_compiled.csv")
precip_compiled <- read_csv(file = "data/neon_precip_compiled.csv")


# 3.1 RIEM temperatures ---------------------------------------------------

# MHK near KONZ
konz_additional <- riem_measures(station = "MHK",
                                 date_start = "2014-09-01",
                                 date_end = "2020-09-20")

konz_prop_na <- konz_additional %>%
  add_count(year = year(valid), day = yday(valid)) %>%
  group_by(year, day) %>%
  summarize(total_n = unique(n),
            non_na_count = sum(!is.na(tmpf)),
            prop_na = round(((total_n - non_na_count) / total_n),
                            digits = 2)) %>%
  select(year, day, total_n, prop_na)

# Aggregate data to day level
konz_agg <- konz_additional %>%
  group_by(year = year(valid), day = yday(valid)) %>%
  summarize(mean_tmpf = mean(tmpf, na.rm = TRUE),
            min_tmpf = min(tmpf, na.rm = TRUE),
            max_tmpf = max(tmpf, na.rm = TRUE)) %>%
  transmute(year,
            day,
            mean_temp_c = conv_unit(mean_tmpf, from = "F", to = "C"),
            min_temp_c = conv_unit(min_tmpf, from = "F", to = "C"),
            max_temp_c = conv_unit(max_tmpf, from = "F", to = "C")) %>%
  ungroup()

# Join aggregated data with NA proportions
konz_output <- left_join(
  x = konz_agg,
  y = konz_prop_na,
  by = c("year", "day")) %>%
  mutate(
    # Tag with site ID
    site_id = "konz",
    date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
  # Rearrange for easier reading
  select(site_id, year, day, date, total_n, prop_na, everything()) %>%
  arrange(date)


konz_output

temp_compiled

inner_join(x = temp_compiled,
           y = konz_output,
           by = c("site" = "site_id", "date", "year", "day"),
           suffix = c("_full", "_konz")) %>%
  ggplot() +
  geom_point(
    aes(x = mean_temp_c_full, y = mean_temp_c_konz, fill = as.factor(year)),
    color = "black", shape = 21) +
  xlab("NEON temp, C") +
  ylab("ASOS temp, C") +
  facet_wrap(vars(year)) +
  scale_fill_viridis_d("Year") +
  theme_bw() +
  ggtitle("Site: KONZ")


# LWC near UKFS
ukfs_additional <- riem_measures(station = "LWC",
                                 date_start = "2014-09-01",
                                 date_end = "2020-09-27")

ukfs_prop_na <- ukfs_additional %>%
  add_count(year = year(valid), day = yday(valid)) %>%
  group_by(year, day) %>%
  summarize(total_n = unique(n),
            non_na_count = sum(!is.na(tmpf)),
            prop_na = round(((total_n - non_na_count) / total_n),
                            digits = 2)) %>%
  select(year, day, total_n, prop_na)

# Aggregate data to day level
ukfs_agg <- ukfs_additional %>%
  group_by(year = year(valid), day = yday(valid)) %>%
  summarize(mean_tmpf = mean(tmpf, na.rm = TRUE),
            min_tmpf = min(tmpf, na.rm = TRUE),
            max_tmpf = max(tmpf, na.rm = TRUE)) %>%
  transmute(year,
            day,
            mean_temp_c = conv_unit(mean_tmpf, from = "F", to = "C"),
            min_temp_c = conv_unit(min_tmpf, from = "F", to = "C"),
            max_temp_c = conv_unit(max_tmpf, from = "F", to = "C")) %>%
  ungroup()

# Join aggregated data with NA proportions
ukfs_output <- left_join(
  x = ukfs_agg,
  y = ukfs_prop_na,
  by = c("year", "day")) %>%
  mutate(
    # Tag with site ID
    site_id = "ukfs",
    date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
  # Rearrange for easier reading
  select(site_id, year, day, date, total_n, prop_na, everything()) %>%
  arrange(date)

inner_join(x = temp_compiled,
           y = ukfs_output,
           by = c("site" = "site_id", "date", "year", "day"),
           suffix = c("_full", "_ukfs")) %>%
  ggplot() +
  geom_point(
    aes(x = mean_temp_c_full, y = mean_temp_c_ukfs, fill = as.factor(year)),
    color = "black", shape = 21) +
  xlab("NEON temp, C") +
  ylab("ASOS temp, C") +
  facet_wrap(vars(year)) +
  scale_fill_viridis_d("Year") +
  theme_bw() +
  ggtitle("Site: UKFS")



# OQT near ORNL
ornl_additional <- riem_measures(station = "OQT",
                                 date_start = "2013-09-01",
                                 date_end = "2020-09-06")

ornl_prop_na <- ornl_additional %>%
  add_count(year = year(valid), day = yday(valid)) %>%
  group_by(year, day) %>%
  summarize(total_n = unique(n),
            non_na_count = sum(!is.na(tmpf)),
            prop_na = round(((total_n - non_na_count) / total_n),
                            digits = 2)) %>%
  select(year, day, total_n, prop_na)

# Aggregate data to day level
ornl_agg <- ornl_additional %>%
  group_by(year = year(valid), day = yday(valid)) %>%
  summarize(mean_tmpf = mean(tmpf, na.rm = TRUE),
            min_tmpf = min(tmpf, na.rm = TRUE),
            max_tmpf = max(tmpf, na.rm = TRUE)) %>%
  transmute(year,
            day,
            mean_temp_c = conv_unit(mean_tmpf, from = "F", to = "C"),
            min_temp_c = conv_unit(min_tmpf, from = "F", to = "C"),
            max_temp_c = conv_unit(max_tmpf, from = "F", to = "C")) %>%
  ungroup()

# Join aggregated data with NA proportions
ornl_output <- left_join(
  x = ornl_agg,
  y = ornl_prop_na,
  by = c("year", "day")) %>%
  mutate(
    # Tag with site ID
    site_id = "ornl",
    date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
  # Rearrange for easier reading
  select(site_id, year, day, date, total_n, prop_na, everything()) %>%
  arrange(date)

inner_join(x = temp_compiled,
           y = ornl_output,
           by = c("site" = "site_id", "date", "year", "day"),
           suffix = c("_full", "_ornl")) %>%
  ggplot() +
  geom_point(
    aes(x = mean_temp_c_full, y = mean_temp_c_ornl, fill = as.factor(year)),
    color = "black", shape = 21) +
  xlab("NEON temp, C") +
  ylab("ASOS temp, C") +
  facet_wrap(vars(year)) +
  scale_fill_viridis_d("Year") +
  theme_bw() +
  ggtitle("Site: ORNL")



# TCL near TALL
tall_additional <- riem_measures(station = "TCL",
                                 date_start = "2013-09-01",
                                 date_end = "2020-09-20")

tall_prop_na <- tall_additional %>%
  add_count(year = year(valid), day = yday(valid)) %>%
  group_by(year, day) %>%
  summarize(total_n = unique(n),
            non_na_count = sum(!is.na(tmpf)),
            prop_na = round(((total_n - non_na_count) / total_n),
                            digits = 2)) %>%
  select(year, day, total_n, prop_na)

# Aggregate data to day level
tall_agg <- tall_additional %>%
  group_by(year = year(valid), day = yday(valid)) %>%
  summarize(mean_tmpf = mean(tmpf, na.rm = TRUE),
            min_tmpf = min(tmpf, na.rm = TRUE),
            max_tmpf = max(tmpf, na.rm = TRUE)) %>%
  transmute(year,
            day,
            mean_temp_c = conv_unit(mean_tmpf, from = "F", to = "C"),
            min_temp_c = conv_unit(min_tmpf, from = "F", to = "C"),
            max_temp_c = conv_unit(max_tmpf, from = "F", to = "C")) %>%
  ungroup()

# Join aggregated data with NA proportions
tall_output <- left_join(
  x = tall_agg,
  y = tall_prop_na,
  by = c("year", "day")) %>%
  mutate(
    # Tag with site ID
    site_id = "tall",
    date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
  # Rearrange for easier reading
  select(site_id, year, day, date, total_n, prop_na, everything()) %>%
  arrange(date)

inner_join(x = temp_compiled,
           y = tall_output,
           by = c("site" = "site_id", "date", "year", "day"),
           suffix = c("_full", "_tall")) %>%
  ggplot() +
  geom_point(
    aes(x = mean_temp_c_full, y = mean_temp_c_tall, fill = as.factor(year)),
    color = "black", shape = 21) +
  xlab("NEON temp, C") +
  ylab("ASOS temp, C") +
  facet_wrap(vars(year)) +
  scale_fill_viridis_d("Year") +
  theme_bw() +
  ggtitle("Site: TALL")



# DYA kind of near LENO
leno_additional <- riem_measures(station = "DYA",
                                 date_start = "2015-09-01",
                                 date_end = "2019-08-04")

leno_prop_na <- leno_additional %>%
  add_count(year = year(valid), day = yday(valid)) %>%
  group_by(year, day) %>%
  summarize(total_n = unique(n),
            non_na_count = sum(!is.na(tmpf)),
            prop_na = round(((total_n - non_na_count) / total_n),
                            digits = 2)) %>%
  select(year, day, total_n, prop_na)

# Aggregate data to day level
leno_agg <- leno_additional %>%
  group_by(year = year(valid), day = yday(valid)) %>%
  summarize(mean_tmpf = mean(tmpf, na.rm = TRUE),
            min_tmpf = min(tmpf, na.rm = TRUE),
            max_tmpf = max(tmpf, na.rm = TRUE)) %>%
  transmute(year,
            day,
            mean_temp_c = conv_unit(mean_tmpf, from = "F", to = "C"),
            min_temp_c = conv_unit(min_tmpf, from = "F", to = "C"),
            max_temp_c = conv_unit(max_tmpf, from = "F", to = "C")) %>%
  ungroup()

# Join aggregated data with NA proportions
leno_output <- left_join(
  x = leno_agg,
  y = leno_prop_na,
  by = c("year", "day")) %>%
  mutate(
    # Tag with site ID
    site_id = "leno",
    date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
  # Rearrange for easier reading
  select(site_id, year, day, date, total_n, prop_na, everything()) %>%
  arrange(date)

inner_join(x = temp_compiled,
           y = leno_output,
           by = c("site" = "site_id", "date", "year", "day"),
           suffix = c("_full", "_leno")) %>%
  ggplot() +
  geom_point(
    aes(x = mean_temp_c_full, y = mean_temp_c_leno, fill = as.factor(year)),
    color = "black", shape = 21) +
  xlab("NEON temp, C") +
  ylab("ASOS temp, C") +
  facet_wrap(vars(year)) +
  scale_fill_viridis_d("Year") +
  theme_bw() +
  ggtitle("Site: LENO")



# 42J near OSBS
osbs_additional <- riem_measures(station = "42J",
                                 date_start = "2013-09-01",
                                 date_end = "2020-07-12")

osbs_prop_na <- osbs_additional %>%
  add_count(year = year(valid), day = yday(valid)) %>%
  group_by(year, day) %>%
  summarize(total_n = unique(n),
            non_na_count = sum(!is.na(tmpf)),
            prop_na = round(((total_n - non_na_count) / total_n),
                            digits = 2)) %>%
  select(year, day, total_n, prop_na)

# Aggregate data to day level
osbs_agg <- osbs_additional %>%
  group_by(year = year(valid), day = yday(valid)) %>%
  summarize(mean_tmpf = mean(tmpf, na.rm = TRUE),
            min_tmpf = min(tmpf, na.rm = TRUE),
            max_tmpf = max(tmpf, na.rm = TRUE)) %>%
  transmute(year,
            day,
            mean_temp_c = conv_unit(mean_tmpf, from = "F", to = "C"),
            min_temp_c = conv_unit(min_tmpf, from = "F", to = "C"),
            max_temp_c = conv_unit(max_tmpf, from = "F", to = "C")) %>%
  ungroup()

# Join aggregated data with NA proportions
osbs_output <- left_join(
  x = osbs_agg,
  y = osbs_prop_na,
  by = c("year", "day")) %>%
  mutate(
    # Tag with site ID
    site_id = "osbs",
    date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
  # Rearrange for easier reading
  select(site_id, year, day, date, total_n, prop_na, everything()) %>%
  arrange(date)

inner_join(x = temp_compiled,
           y = osbs_output,
           by = c("site" = "site_id", "date", "year", "day"),
           suffix = c("_full", "_osbs")) %>%
  ggplot() +
  geom_point(
    aes(x = mean_temp_c_full, y = mean_temp_c_osbs, fill = as.factor(year)),
    color = "black", shape = 21) +
  xlab("NEON temp, C") +
  ylab("ASOS temp, C") +
  facet_wrap(vars(year)) +
  scale_fill_viridis_d("Year") +
  theme_bw() +
  ggtitle("Site: OSBS")


# FRR near SCBI
scbi_additional <- riem_measures(station = "FRR",
                                 date_start = "2013-09-01",
                                 date_end = "2020-10-04")

scbi_prop_na <- scbi_additional %>%
  add_count(year = year(valid), day = yday(valid)) %>%
  group_by(year, day) %>%
  summarize(total_n = unique(n),
            non_na_count = sum(!is.na(tmpf)),
            prop_na = round(((total_n - non_na_count) / total_n),
                            digits = 2)) %>%
  select(year, day, total_n, prop_na)

# Aggregate data to day level
scbi_agg <- scbi_additional %>%
  group_by(year = year(valid), day = yday(valid)) %>%
  summarize(mean_tmpf = mean(tmpf, na.rm = TRUE),
            min_tmpf = min(tmpf, na.rm = TRUE),
            max_tmpf = max(tmpf, na.rm = TRUE)) %>%
  transmute(year,
            day,
            mean_temp_c = conv_unit(mean_tmpf, from = "F", to = "C"),
            min_temp_c = conv_unit(min_tmpf, from = "F", to = "C"),
            max_temp_c = conv_unit(max_tmpf, from = "F", to = "C")) %>%
  ungroup()

# Join aggregated data with NA proportions
scbi_output <- left_join(
  x = scbi_agg,
  y = scbi_prop_na,
  by = c("year", "day")) %>%
  mutate(
    # Tag with site ID
    site_id = "scbi",
    date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
  # Rearrange for easier reading
  select(site_id, year, day, date, total_n, prop_na, everything()) %>%
  arrange(date)

inner_join(x = temp_compiled,
           y = scbi_output,
           by = c("site" = "site_id", "date", "year", "day"),
           suffix = c("_full", "_scbi")) %>%
  ggplot() +
  geom_point(
    aes(x = mean_temp_c_full, y = mean_temp_c_scbi, fill = as.factor(year)),
    color = "black", shape = 21) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("NEON temp, C") +
  ylab("ASOS temp, C") +
  facet_wrap(vars(year)) +
  scale_fill_viridis_d("Year") +
  theme_bw() +
  ggtitle("Site: SCBI")


# FRR near BLAN
blan_additional <- riem_measures(station = "FRR",
                                 date_start = "2014-09-01",
                                 date_end = "2020-09-13")

blan_prop_na <- blan_additional %>%
  add_count(year = year(valid), day = yday(valid)) %>%
  group_by(year, day) %>%
  summarize(total_n = unique(n),
            non_na_count = sum(!is.na(tmpf)),
            prop_na = round(((total_n - non_na_count) / total_n),
                            digits = 2)) %>%
  select(year, day, total_n, prop_na)

# Aggregate data to day level
blan_agg <- blan_additional %>%
  group_by(year = year(valid), day = yday(valid)) %>%
  summarize(mean_tmpf = mean(tmpf, na.rm = TRUE),
            min_tmpf = min(tmpf, na.rm = TRUE),
            max_tmpf = max(tmpf, na.rm = TRUE)) %>%
  transmute(year,
            day,
            mean_temp_c = conv_unit(mean_tmpf, from = "F", to = "C"),
            min_temp_c = conv_unit(min_tmpf, from = "F", to = "C"),
            max_temp_c = conv_unit(max_tmpf, from = "F", to = "C")) %>%
  ungroup()

# Join aggregated data with NA proportions
blan_output <- left_join(
  x = blan_agg,
  y = blan_prop_na,
  by = c("year", "day")) %>%
  mutate(
    # Tag with site ID
    site_id = "blan",
    date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
  # Rearrange for easier reading
  select(site_id, year, day, date, total_n, prop_na, everything()) %>%
  arrange(date)

inner_join(x = temp_compiled,
           y = blan_output,
           by = c("site" = "site_id", "date", "year", "day"),
           suffix = c("_full", "_blan")) %>%
  ggplot() +
  geom_point(
    aes(x = mean_temp_c_full, y = mean_temp_c_blan, fill = as.factor(year)),
    color = "black", shape = 21) +
  xlab("NEON temp, C") +
  ylab("ASOS temp, C") +
  facet_wrap(vars(year)) +
  scale_fill_viridis_d("Year") +
  theme_bw() +
  ggtitle("Site: BLAN")


# NAK kind of near SERC (Not very close)
serc_additional <- riem_measures(station = "NAK",
                                 date_start = "2014-09-01",
                                 date_end = "2020-09-20")

serc_prop_na <- serc_additional %>%
  add_count(year = year(valid), day = yday(valid)) %>%
  group_by(year, day) %>%
  summarize(total_n = unique(n),
            non_na_count = sum(!is.na(tmpf)),
            prop_na = round(((total_n - non_na_count) / total_n),
                            digits = 2)) %>%
  select(year, day, total_n, prop_na)

# Aggregate data to day level
serc_agg <- serc_additional %>%
  group_by(year = year(valid), day = yday(valid)) %>%
  summarize(mean_tmpf = mean(tmpf, na.rm = TRUE),
            min_tmpf = min(tmpf, na.rm = TRUE),
            max_tmpf = max(tmpf, na.rm = TRUE)) %>%
  transmute(year,
            day,
            mean_temp_c = conv_unit(mean_tmpf, from = "F", to = "C"),
            min_temp_c = conv_unit(min_tmpf, from = "F", to = "C"),
            max_temp_c = conv_unit(max_tmpf, from = "F", to = "C")) %>%
  ungroup()

# Join aggregated data with NA proportions
serc_output <- left_join(
  x = serc_agg,
  y = serc_prop_na,
  by = c("year", "day")) %>%
  mutate(
    # Tag with site ID
    site_id = "serc",
    date = (ymd(paste0(year, "-01-01")) + day) - 1) %>%
  # Rearrange for easier reading
  select(site_id, year, day, date, total_n, prop_na, everything()) %>%
  arrange(date)

inner_join(x = temp_compiled,
           y = serc_output,
           by = c("site" = "site_id", "date", "year", "day"),
           suffix = c("_full", "_serc")) %>%
  ggplot() +
  geom_point(
    aes(x = mean_temp_c_full, y = mean_temp_c_serc, fill = as.factor(year)),
    color = "black", shape = 21) +
  xlab("NEON temp, C") +
  ylab("ASOS temp, C") +
  facet_wrap(vars(year)) +
  scale_fill_viridis_d("Year") +
  theme_bw() +
  ggtitle("Site: SERC")


# Combine the RIEM weather data pulls for all plots
full_riem <- reduce(.x = list(blan_output, konz_output, leno_output, ornl_output,
                              osbs_output, scbi_output, serc_output, tall_output,
                              ukfs_output),
                    .f = bind_rows)

write_csv(x = full_riem,
          file = "data/temperature_from_riem.csv")


# 3.2 RIEM Precipitation --------------------------------------------------

# Exploring this idea for now, but not pursuing in full

# MHK near KONZ

# The ASOS data are a little more difficult to aggregate
# Quoting riem::riem_measures() -- "One hour precipitation for the period from
# the observation time to the time of the previous hourly precipitation reset.
# This varies slightly by site."
# Help with this from https://github.com/rmcd1024/daily_precipitation_totals
konz_precip <- konz_additional %>% 
  filter(!is.na(p01i)) %>% 
  select(station, valid, p01i) %>% 
  mutate(date = with_tz(as.POSIXct(valid, tz ='UCT'),
                        "America/New_York"),
         date = date - dst(date)*3600,
         year = year(date),
         month = month(date),
         day = day(date),
         hour = hour(date),
         minute = minute(date)) %>% 
  group_by(year, month, day, hour) 

modalminute <- konz_precip %>%
  arrange(desc(minute), .by_group = TRUE) %>% 
  mutate(maxminute = minute[which.max(p01i)]) %>% 
  {which.max(tabulate(.$maxminute))}

konz_precip_day <- konz_precip %>% 
  filter(minute <= modalminute) %>% 
  summarize(hourlyprecip = max(p01i, na.rm = TRUE)) %>% 
  group_by(year, month, day) %>% 
  summarize(sum_precip_mm = sum(hourlyprecip, na.rm = TRUE) * 25.4) %>%
  ungroup() %>%
  mutate(site_id = "konz",
         date = ymd(paste(year, month, day, sep = "-")))

# Hmm...could be worse I guess
inner_join(x = precip_compiled,
           y = konz_precip_day,
           by = c("site" = "site_id", "date", "year", "day"),
           suffix = c("_full", "_konz")) %>%
  ggplot() +
  geom_point(
    aes(x = sum_precip_mm_full , y = sum_precip_mm_konz, fill = as.factor(year)),
    color = "black", shape = 21) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("NEON precip, mm") +
  ylab("ASOS precip, mm") +
  facet_wrap(vars(year)) +
  scale_fill_viridis_d("Year") +
  theme_bw() +
  ggtitle("Site: KONZ")


# 3.3 PRISM VPD -----------------------------------------------------------

# Download PRISM min/max vapor pressure deficit in bulk

prism_set_dl_dir(path = "data/prism/")

# Minimum Vapor Pressure Deficit
get_prism_dailys(type = "vpdmin",
                 minDate = "2013-01-01",
                 maxDate = "2021-12-31")

# Maximum Vapor Pressure Deficit
get_prism_dailys(type = "vpdmax",
                 minDate = "2013-01-01",
                 maxDate = "2021-12-31")

neon_sites_sf <- neon_sites %>%
  select(field_domain_id, field_site_id, field_site_name, field_site_state,
         field_latitude, field_longitude) %>%
  st_as_sf(x = .,
           coords = c("field_longitude", "field_latitude"),
           # WGS84
           crs = 4326)

# Get raster list
file_archive <- list.files(path = "data/prism", pattern = ".bil$",
                           recursive = T, full.names = TRUE)

vpdmin_archive <- grep(pattern = "vpdmin", x = file_archive, value = TRUE)
vpdmax_archive <- grep(pattern = "vpdmax", x = file_archive, value = TRUE)


vpd_min <- map_df(.x = vpdmin_archive,
                  .f = ~ {
                    
                    date_stamp <- .x %>%
                      str_extract(string = ., pattern = "[0-9]{8}_bil.bil") %>%
                      str_remove(string = ., pattern = "_bil.bil")
                    
                    prism_raster <- raster(x = .x)
                    
                    vpdmin_vals <- extract(x = prism_raster,
                                           y = neon_sites_sf)
                    
                    out_df <- tibble(neon_sites_sf) %>%
                      dplyr::select(-geometry) %>%
                      mutate(date = ymd(date_stamp)) %>%
                      bind_cols(., vpd_min = vpdmin_vals)
                    
                    return(out_df)
                    
                  })

write_rds(x = vpd_min,
          file = "data/vpd_min_prism.rds")

vpd_max <- map_df(.x = vpdmax_archive,
                  .f = ~ {
                    
                    date_stamp <- .x %>%
                      str_extract(string = ., pattern = "[0-9]{8}_bil.bil") %>%
                      str_remove(string = ., pattern = "_bil.bil")
                    
                    prism_raster <- raster(x = .x)
                    
                    vpdmax_vals <- extract(x = prism_raster,
                                           y = neon_sites_sf)
                    
                    out_df <- tibble(neon_sites_sf) %>%
                      dplyr::select(-geometry) %>%
                      mutate(date = ymd(date_stamp)) %>%
                      bind_cols(., vpd_max = vpdmax_vals)
                    
                    return(out_df)
                    
                  })

write_rds(x = vpd_max,
          file = "data/vpd_max_prism.rds")


# 3.4 GridMet RH ----------------------------------------------------------

# We'll use GridMet's RH, see how it matches what we've gotten from NEON,
# then try out our VPD calculations with it. I've also downloaded VPD
# from GridMet in case we want to use that

# I've downloaded GridMet .nc (NetCDF) files, from which these values can
# be pulled

# A function to be used in mapping over NetCDF files to retrieve and format
# past climate data
extract_from_brick <- function(brick_file, sites_sf, origin_date, variable){
  
  # brick_file = A file to be read in as a raster brick of climate data (gridmet intended)
  # sites_sf = An sf object containing NEON site locations
  # origin_date = The origin date of the .nc files being read in, as
  #               a character string. (e.g., when 'description: days
  #               since 1900-01-01', then the value would be "1900-01-01")
  # variable = A character string containing the name of the variable being
  #            extracted. This is not used to retrieve the variable, just to
  #            label its column in the output data frame 
  
  raster_brick <- raster::brick(x = brick_file)
  
  # Ensure that they are in compatible reference systems
  assert_that(st_crs(raster_brick) == st_crs(sites_sf),
              msg = "CRS of brick_file != CRS of sites_sf")  
  
  raster::extract(raster_brick, sites_sf, method = "simple", df = TRUE) %>%
    pivot_longer(names_to = "raw_dates", values_to = variable, cols = -"ID") %>%
    mutate(date = gsub(pattern = "X", replacement = "", x = raw_dates),
           date = as.numeric(date),
           date = ymd(origin_date) + days(date),
           ID = as.character(ID)) %>%
    left_join(x = .,
              y = as_tibble(sites_sf) %>%
                rownames_to_column() %>%
                dplyr::select(rowname, field_site_id),
              by = c("ID" = "rowname")) %>%
    select(field_site_id, date, all_of(variable))
  
}

# Iterate over gridmet files and compile daily estimates for each site's
# location into a data frame
gridmet_rh_max <- map_df(.x = list.files(path = "data/gridmet/",
                                         pattern = "rmax_.{4}\\.nc",
                                         full.names = TRUE),
                         .f = ~ extract_from_brick(brick_file = .x,
                                                   sites_sf = neon_sites_sf,
                                                   origin_date = "1900-01-01",
                                                   variable = "rh_max"))

gridmet_rh_min <- map_df(.x = list.files(path = "data/gridmet/",
                                         pattern = "rmin_.{4}\\.nc",
                                         full.names = TRUE),
                         .f = ~ extract_from_brick(brick_file = .x,
                                                   sites_sf = neon_sites_sf,
                                                   origin_date = "1900-01-01",
                                                   variable = "rh_min"))

gridmet_vpd <- map_df(.x = list.files(path = "data/gridmet/",
                                      pattern = "vpd_.{4}\\.nc",
                                      full.names = TRUE),
                      .f = ~ extract_from_brick(brick_file = .x,
                                                sites_sf = neon_sites_sf,
                                                origin_date = "1900-01-01",
                                                variable = "vpd"))

# Combine all variables and then export
full_gridmet <- reduce(.x = list(gridmet_rh_min, gridmet_rh_max, gridmet_vpd),
                       .f = full_join,
                       by = c("field_site_id", "date"))

write_csv(x = full_gridmet,
          file = "data/gridmet_data_compiled.csv")

