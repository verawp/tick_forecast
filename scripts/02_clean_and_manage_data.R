library(tidyverse)
library(lubridate)
library(MMWRweek)


# 1. Combine NEON tick and NEON weather data ------------------------------

# The datasets we'll use:
tick_targets <- read_csv(file = "data/ticks_target.csv")
temp_compiled <- read_csv(file = "data/neon_temp_compiled.csv")
rh_compiled <- read_csv(file = "data/neon_rh_compiled.csv")
precip_compiled <- read_csv(file = "data/neon_precip_compiled.csv")

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

# Export for external use
write_csv(x = tick_neon_riem_clean,
          file = "data/ticks_with_filled_weather.csv")



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


# 3. Assess coverage of non-temp vars -------------------------------------

# What's RH coverage looking like?
ggplot(data = tick_neon_riem_clean) +
  geom_point(aes(x = time, y = mean_rh_pct), alpha = 0.45, color = "black",
             fill = "gray20", shape = 21, size = 2) +
  geom_point(data = tick_neon_riem_clean %>%
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
  ggtitle("Weekly RH data completion, NEON")

# What's precip coverage looking like?
ggplot(data = tick_neon_riem_clean) +
  geom_point(aes(x = time, y = sum_precip_mm ), alpha = 0.35, color = "black",
             fill = "gray20", shape = 21, size = 2) +
  geom_point(data = tick_neon_riem_clean %>%
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






