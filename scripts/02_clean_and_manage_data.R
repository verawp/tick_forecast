library(tidyverse)
library(lubridate)
library(MMWRweek)
library(plantecophys)


# 1. Combine NEON tick and NEON weather data ------------------------------

# The datasets we'll use:
tick_targets <- read_csv(file = "data/ticks_target.csv")
temp_compiled <- read_csv(file = "data/neon_temp_compiled.csv")
rh_compiled <- read_csv(file = "data/neon_rh_compiled.csv")
precip_compiled <- read_csv(file = "data/neon_precip_compiled.csv")
vpd_min_prism <- read_rds(file = "data/vpd_min_prism.rds")
vpd_max_prism <- read_rds(file = "data/vpd_max_prism.rds")
full_gridmet <- read_csv(file = "data/gridmet_data_compiled.csv")

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
              sum_precip_mm = sum(sum_precip_mm, na.rm = TRUE),
              precip_type = unique(precip_type))),
  .f = full_join,
  by = c("site_id", "year", "mmwr_week")) %>%
  mutate(
    # Replace Infinite values and NaN with NAs
    across(.cols = all_of(cols_to_fix),
           .fns =  ~na_if(., Inf)),
    across(.cols = all_of(cols_to_fix),
           .fns =  ~na_if(., -Inf)),
    across(.cols = all_of(cols_to_fix),
           .fns = ~if_else(condition = is.nan(.),
                           true = NA_real_,
                           false = .))) %>%
  distinct()

# Note that there will be NAs in the mmwr_week and day_num columns where the
# data goes beyond the limits of the tick data sampling
write_csv(x = full_dataset,
          file = "data/ticks_with_neon_weather.csv")


# 2. Combine NEON data with RIEM ------------------------------------------

# The temperature data from the RIEM package (in script 01)
full_riem <- read_csv(file = "data/temperature_from_riem.csv")

# Summarize the temp data from RIEM to the weekly level
riem_weekly <- full_riem %>%
  rename(prop_na_temp = prop_na) %>%
  group_by(site_id, year, mmwr_week = epiweek(date)) %>%
  summarize(site_id = toupper(site_id),
            mean_temp = mean(mean_temp_c, na.rm = TRUE),
            min_temp = min(min_temp_c, na.rm = TRUE),
            max_temp = max(max_temp_c, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    # Replace Infinite values and NaN with NAs
    across(.cols = c("mean_temp", "min_temp", "max_temp"),
           .fns =  ~na_if(., Inf)),
    across(.cols = c("mean_temp", "min_temp", "max_temp"),
           .fns =  ~na_if(., -Inf)),
    across(.cols = c("mean_temp", "min_temp", "max_temp"),
           .fns = ~if_else(condition = is.nan(.),
                           true = NA_real_,
                           false = .))) %>%
  distinct()

# Merge the RIEM temperature data with the NEON tick/weather dataset. But only
# use the RIEM temperature data if there is no NEON available...
tick_neon_riem <- full_join(
  x = full_dataset,
  y = riem_weekly,
  by = c("site_id", "year", "mmwr_week"),
  # Columns with identical names will be given suffixes based on which dataset
  # they came from:
  suffix = c("_neon", "_riem")) %>%
  # Use the RIEM-sourced data when NEON is NA:
  mutate(
    mean_temp = case_when(
      # If NEON is NA, then use RIEM
      is.na(mean_temp_neon) ~ mean_temp_riem,
      # In all other cases use NEON
      TRUE ~ mean_temp_neon),
    min_temp = case_when(
      is.na(min_temp_neon) ~ min_temp_riem,
      TRUE ~ min_temp_neon),
    max_temp = case_when(
      is.na(max_temp_neon) ~ max_temp_riem,
      TRUE ~ max_temp_neon),
    
    # Round to two decimal places
    across(.cols = c("amblyomma_americanum", "mean_temp", "min_temp",
                     "max_temp",  "mean_rh_pct", "min_rh_pct",
                     "max_rh_pct", "mean_var_rh", "sum_precip_mm"),
           .fns = ~round(x = ., digits = 2)),
    # Record which source the temp data came from
    temp_source = case_when(
      
      # Condition: All NEON columns are NA ...AND...
      (is.na(mean_temp_neon) & is.na(min_temp_neon) & is.na(max_temp_neon)) &
        # ...none of the RIEM columns are NA
        (!is.na(mean_temp_riem) & !is.na(min_temp_riem) & !is.na(max_temp_riem)) ~ "RIEM",
      
      # Condition: All of the NEON columns are NA ...AND...
      (is.na(mean_temp_neon) & is.na(min_temp_neon) & is.na(max_temp_neon)) &
        # ...all of the RIEM columns are NA as well
        (is.na(mean_temp_riem) & is.na(min_temp_riem) & is.na(max_temp_riem)) ~ "All temp. sources NA",
      
      # Condition: None of the NEON columns are NA
      !is.na(mean_temp_neon) & !is.na(min_temp_neon) & !is.na(max_temp_neon) ~ "NEON",
      
      # If one of the above options doesn't apply...
      TRUE ~ "Unexpected values"),
    # Fill in missing dates based on their MMWR week
    time = if_else(condition = is.na(time),
                   true = MMWRweek2Date(MMWRyear = year, MMWRweek = mmwr_week),
                   false = time),
    day_num = if_else(condition = is.na(day_num),
                      true = yday(time),
                      false = day_num)) %>%
  select(site_id, time, year, day_num, mmwr_week, amblyomma_americanum,
         mean_temp, min_temp, max_temp, temp_source, contains("_rh_"),
         contains("precip"))

# There's some unneeded data in here, so I am going to go through and filter
# out the data so that there's not weather data outside of the September before
# the first tick sample at each site (so we can look at winter effects), and none
# after sample records stop.
str(tick_neon_riem)

konz_limits <- interval(start = "2014-09-01",
                        end = "2020-09-20")

ornl_limits <- interval(start = "2013-09-01",
                        end = "2020-09-06")

tall_limits <- interval(start = "2013-09-01",
                        end = "2020-09-20")

leno_limits <- interval(start = "2015-09-01",
                        end = "2019-08-04")

osbs_limits <- interval(start = "2013-09-01",
                        end = "2020-07-12")

scbi_limits <- interval(start = "2013-09-01",
                        end = "2020-10-04")

blan_limits <- interval(start = "2014-09-01",
                        end = "2020-09-13")

serc_limits <- interval(start = "2014-09-01",
                        end = "2020-09-20")

ukfs_limits <- interval(start = "2014-09-01",
                        end = "2020-09-27")


# The filtered version of the dataset:
tick_neon_riem_clean <- tick_neon_riem %>%
  filter(
    (site_id == "KONZ" & time %within% konz_limits) |
      (site_id == "ORNL" & time %within% ornl_limits) |
      (site_id == "TALL" & time %within% tall_limits) |
      (site_id == "LENO" & time %within% leno_limits) |
      (site_id == "OSBS" & time %within% osbs_limits) |
      (site_id == "SCBI" & time %within% scbi_limits) |
      (site_id == "BLAN" & time %within% blan_limits) |
      (site_id == "SERC" & time %within% serc_limits) |
      (site_id == "UKFS" & time %within% ukfs_limits)) %>%
  arrange(site_id, time) 


# Anything weird happening with these different sources?
ggplot(data = tick_neon_riem_clean %>%
         filter(temp_source %in% c("NEON", "RIEM"))) +
  geom_point(aes(x = time, y = mean_temp, fill = temp_source),
             color = "black", alpha = 0.45, shape = 21, size = 2) +
  geom_point(data = tick_neon_riem_clean %>%
               filter(temp_source == "All temp. sources NA"),
             aes(x = time, y = 0, size = ""), color = "red",
             shape = 4) +
  scale_fill_viridis_d("Temperature Source") +
  scale_size_discrete("All temp. sources NA") +
  xlab("Date") +
  ylab("Mean temperature, C") +
  facet_wrap(vars(site_id), scales = "free", ncol = 2) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold")) +
  ggtitle("Weekly temperature data completion by provider")

# Note that SCBI seems to have some inaccurate data in early 2019...


# 3. Calculate air VPD ----------------------------------------------------

# Calculate VPD using existing NEON RH values
full_dataset_w_vpd <- tick_neon_riem_clean %>%
  mutate(min_air_vpd = RHtoVPD(RH = min_rh_pct, TdegC = min_temp),
         mean_air_vpd = RHtoVPD(RH = mean_rh_pct, TdegC = mean_temp),
         max_air_vpd = RHtoVPD(RH = max_rh_pct, TdegC = max_temp))

# Well...could be better
full_join(x = full_dataset_w_vpd,
          y = vpd_min_prism %>%
            rename(site_id = field_site_id) %>%
            group_by(site_id, year = year(date), mmwr_week = epiweek(date)) %>%
            summarize(min_vpd = min(vpd_min, na.rm = TRUE)),
          by = c("site_id", "year", "mmwr_week")) %>%
  ggplot() +
  geom_point(aes(x = min_air_vpd, y = min_vpd)) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(vars(site_id)) +
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  xlab("Manual VPD, min.") +
  ylab("PRISM VPD, min")


# What about using GridMet?
gridmet_vpd <- full_join(x = full_gridmet %>%
                           rowwise() %>%
                           mutate(rh_mean = mean(c(rh_min, rh_max))) %>%
                           ungroup() %>%
                           rename(site_id = field_site_id),
                         y = temp_compiled %>%
                           mutate(site_id = toupper(site)) %>%
                           select(site_id, date, mean_temp_c, min_temp_c, max_temp_c),
                         by = c("site_id", "date"))  %>%
  mutate(min_air_vpd = RHtoVPD(RH = rh_min , TdegC = min_temp_c),
         mean_air_vpd = RHtoVPD(RH = rh_mean, TdegC = mean_temp_c),
         max_air_vpd = RHtoVPD(RH = rh_max, TdegC = max_temp_c)) %>%
  group_by(site_id, year = year(date), mmwr_week = epiweek(date)) %>%
  summarize(
    # VPD direct from GridMet
    min_vpd_gridmet = min(vpd, na.rm = TRUE),
    mean_vpd_gridmet = mean(vpd, na.rm = TRUE),
    max_vpd_gridmet = max(vpd, na.rm = TRUE),
    # VPD made from GridMet + NEON
    min_vpd_mix = min(min_air_vpd, na.rm = TRUE),
    mean_vpd_mix = mean(mean_air_vpd, na.rm = TRUE),
    max_vpd_mix = max(max_air_vpd, na.rm = TRUE),
    # RH only
    min_rh_gridmet = min(rh_min, na.rm = TRUE),
    mean_rh_gridmet = mean(rh_mean, na.rm = TRUE),
    max_rh_gridmet = max(rh_max, na.rm = TRUE))

# Mixing GridMet with NEON yields decent results
full_join(x = full_dataset_w_vpd %>%
            select(site_id, year, mmwr_week, contains("_vpd")),
          y = gridmet_vpd,
          by = c("site_id", "year", "mmwr_week")) %>%
  ggplot() +
  geom_point(aes(x = min_air_vpd, y = min_vpd_mix)) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(vars(site_id)) +
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  xlab("NEON VPD, min.") +
  ylab("GridMet VPD, min")

# Using GridMet directly yields less impressive results
full_join(x = full_dataset_w_vpd %>%
            select(site_id, year, mmwr_week, contains("_vpd")),
          y = gridmet_vpd,
          by = c("site_id", "year", "mmwr_week")) %>%
  ggplot() +
  geom_point(aes(x = min_air_vpd, y = min_vpd_gridmet)) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(vars(site_id)) +
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  xlab("NEON VPD, min.") +
  ylab("GridMet VPD, min")

# Merge the GirdMet-NEON VPD data with the existing dataset. But only use the
# new VPD data if there is no NEON-derived available
tick_neon_vpd_fill <- full_join(
  x = full_dataset_w_vpd,
  y = gridmet_vpd %>%
    select(site_id, year, mmwr_week, min_vpd_mix, mean_vpd_mix, max_vpd_mix,
           min_rh_gridmet, mean_rh_gridmet, max_rh_gridmet),
  by = c("site_id", "year", "mmwr_week")) %>%
  # Use the RIEM-sourced data when NEON is NA:
  mutate(
    # Record which source the temp data came from
    rh_vpd_source = case_when(
      
      # Condition: All NEON columns are NA ...AND...
      (!is.finite(mean_air_vpd) & !is.finite(min_air_vpd) & !is.finite(max_air_vpd)) &
        # ...none of the GridMet/mixed columns are NA
        (is.finite(mean_vpd_mix) & is.finite(min_vpd_mix) & is.finite(max_vpd_mix)) ~ "GridMet/NEON",
      
      # Condition: All of the NEON columns are NA ...AND...
      (!is.finite(mean_air_vpd) & !is.finite(min_air_vpd) & !is.finite(max_air_vpd)) &
        # ...all of the RIEM columns are NA as well
        (!is.finite(mean_vpd_mix) & !is.finite(min_vpd_mix) & !is.finite(max_vpd_mix)) ~ "All sources NA",
      
      # Condition: None of the NEON columns are NA
      is.finite(mean_air_vpd) & is.finite(min_air_vpd) & is.finite(max_air_vpd) ~ "NEON",
      
      # If one of the above options doesn't apply...
      TRUE ~ "Unexpected values"),
    
    mean_rh_pct = case_when(
      # If NEON is NA, then use Gridmet
      !is.finite(mean_rh_pct) ~ mean_rh_gridmet,
      # In all other cases use NEON
      TRUE ~ mean_rh_pct),
    min_rh_pct = case_when(
      !is.finite(min_rh_pct) ~ min_rh_gridmet,
      TRUE ~ min_rh_pct),
    max_rh_pct = case_when(
      !is.finite(max_rh_pct) ~ max_rh_gridmet,
      TRUE ~ max_rh_pct),
    
    mean_air_vpd = case_when(
      # If NEON is NA, then use Gridmet
      !is.finite(mean_air_vpd) ~ mean_vpd_mix,
      # In all other cases use NEON
      TRUE ~ mean_air_vpd),
    min_air_vpd = case_when(
      !is.finite(min_air_vpd) ~ min_vpd_mix,
      TRUE ~ min_air_vpd),
    max_air_vpd = case_when(
      !is.finite(max_air_vpd) ~ max_vpd_mix,
      TRUE ~ max_air_vpd),
    
    # Fill in missing dates based on their MMWR week
    time = if_else(condition = is.na(time),
                   true = MMWRweek2Date(MMWRyear = year, MMWRweek = mmwr_week),
                   false = time),
    day_num = if_else(condition = is.na(day_num),
                      true = yday(time),
                      false = day_num)) %>%
  select(site_id, time, year, day_num, mmwr_week, amblyomma_americanum,
         mean_temp, min_temp, max_temp, temp_source, contains("rh_pct"),
         contains("precip"), contains("air_vpd"), rh_vpd_source)

# Export for external use
write_csv(x = tick_neon_vpd_fill,
          file = "data/ticks_with_filled_weather.csv")

# NOTE: We could still further fill this dataset with VPD from GridMet. It would
# be better than PRISM I think, but it would still be less accurate than the
# approach of combining GridMet RH with NEON temps


# 4. Assess coverage of variables -----------------------------------------

# What's RH coverage looking like?
ggplot(data = tick_neon_vpd_fill) +
  geom_point(aes(x = time, y = mean_rh_pct), alpha = 0.45, color = "black",
             fill = "gray20", shape = 21, size = 2) +
  geom_point(data = tick_neon_vpd_fill %>%
               filter(is.na(mean_rh_pct)),
             aes(x = time, y = 0, color = "red", shape = 4), 
             alpha = 0.7) +
  # scale_color_manual("NA data for RH", values = "red") +
  scale_color_identity("", guide = "legend", breaks = "red", labels = "NA data for RH") +
  scale_shape_identity("", guide = "legend", breaks = 4, labels = "NA data for RH") +
  xlab("Date") +
  ylab("RH, %") +
  facet_wrap(vars(site_id), scales = "free", ncol = 2) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold")) +
  ggtitle("Weekly mean RH data completion")

# VPD
ggplot(data = tick_neon_vpd_fill) +
  geom_point(aes(x = time, y = mean_air_vpd), alpha = 0.45, color = "black",
             fill = "gray20", shape = 21, size = 2) +
  geom_point(data = tick_neon_vpd_fill %>%
               filter(is.na(mean_air_vpd)),
             aes(x = time, y = 0, color = "red", shape = 4), 
             alpha = 0.7) +
  # scale_color_manual("NA data for RH", values = "red") +
  scale_color_identity("", guide = "legend", breaks = "red", labels = "NA data for VPD") +
  scale_shape_identity("", guide = "legend", breaks = 4, labels = "NA data for VPD") +
  xlab("Date") +
  ylab("Vapor pressure deficit") +
  facet_wrap(vars(site_id), scales = "free", ncol = 2) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold")) +
  ggtitle("Weekly mean VPD data completion")

ggplot(data = tick_neon_vpd_fill %>%
         filter(rh_vpd_source %in% c("NEON", "GridMet/NEON"))) +
  geom_point(aes(x = time, y = mean_air_vpd, fill = rh_vpd_source),
             color = "black", alpha = 0.45, shape = 21, size = 2) +
  geom_point(data = tick_neon_vpd_fill %>%
               filter(rh_vpd_source == "All sources NA"),
             aes(x = time, y = 0, size = ""), color = "red",
             shape = 4) +
  scale_fill_viridis_d("VPD Source") +
  scale_size_discrete("All VPD sources NA") +
  xlab("Date") +
  ylab("Mean VPD") +
  facet_wrap(vars(site_id), scales = "free", ncol = 2) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold")) +
  ggtitle("Weekly VPD data completion by raw data source")


# Precip
ggplot(data = tick_neon_vpd_fill) +
  geom_point(aes(x = time, y = sum_precip_mm ), alpha = 0.35, color = "black",
             fill = "gray20", shape = 21, size = 2) +
  geom_point(data = tick_neon_vpd_fill %>%
               filter(is.na(sum_precip_mm)),
             aes(x = time, y = 0, color = "red", shape = 4),
             alpha = 0.4) +
  scale_color_manual("NA data for RH", values = "red") +
  scale_color_identity("", guide = "legend", breaks = "red", labels = "NA data for precip") +
  scale_shape_identity("", guide = "legend", breaks = 4, labels = "NA data for precip") +
  xlab("Date") +
  ylab("Precip, mm") +
  facet_wrap(vars(site_id), scales = "free", nrow = 2) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold")) +
  ggtitle("Weekly precip data completion, NEON")

# Make a table for easier assessment
tick_neon_vpd_fill %>%
  summarize(across(everything(),
                   is.finite)) %>%
  summarize(across(everything(),
                   ~ round(sum(.) / nrow(tick_neon_vpd_fill),
                           digits = 2))) %>%
  pivot_longer(names_to = "variable", values_to = "pct_complete", everything()) %>%
  filter(!(variable %in% c("site_id", "time", "year", "day_num", "mmwr_week",
                           "temp_source", "rh_vpd_source", "precip_type"))) %>%
  arrange(pct_complete) %>%
  write_csv(file = "data/na_counts_by_col.csv")


# 5. Look at missing tick counts ------------------------------------------

tick_no_field <- read_csv(file = "data/ticks_missed_sample_dates.csv")

tick_neon_vpd_fill

# Weeks between samples
weeks_between <- tick_neon_vpd_fill %>%
  filter(!is.na(amblyomma_americanum)) %>%
  group_by(site_id) %>%
  arrange(time) %>%
  mutate(t_previous_sample = abs(as.numeric(time %--% lag(x = time, n = 1L), "weeks")))

weeks_between %>%
  ggplot() +
  geom_histogram(aes(t_previous_sample),
                 color = "black", fill = "white", binwidth = 1) +
  geom_histogram(data = filter(weeks_between, t_previous_sample %in% c(3, 6)),
                 aes(t_previous_sample),
                 color = "black", fill = "red", binwidth = 1) +
  facet_wrap(vars(site_id), scales = "free_y") +
  theme_bw() +
  xlab("Number of weeks since last non-NA tick value") +
  ylab("Frequency")

