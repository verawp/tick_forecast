# A script to calculate cumulative degree days at the NEON sites using ALL dates
# and not just those where sampling occurred

library(tidyverse)
library(lubridate)
library(zoo)


# 1. Load initial datasets ------------------------------------------------

neon_temps <- read_csv(file = "data/neon_temp_compiled.csv")

head(neon_temps)

# The temperature data from the RIEM package (in script 01)
full_riem <- read_csv(file = "data/temperature_from_riem.csv")

head(full_riem)


# 2. Identify temperature data gaps ---------------------------------------

neon_temps %>% 
  distinct(site, year, date)

# All potential dates covering the years where sampling took place
potential_dates <- neon_temps %>% 
  distinct(site, year) %>%
  # Make each row a list item
  split(sort(as.numeric(rownames(.)))) %>%
  map_df(.x = .,
         .f = ~ {
           
           dates <- seq(from = ymd(paste0(.x$year, "-01-01")),
                        to = ymd(paste0(.x$year, "12-31")),
                        by = 1)
           
           # Output df with site and year included
           tibble(
             site = .x$site,
             year = .x$year,
             date = dates
           )
         })

# Join existing temperature data to the dates above to see what we're missing
neon_with_gaps <- left_join(x = potential_dates,
                            y = neon_temps %>%
                              select(site, year, date,
                                     mean_temp_c, min_temp_c, max_temp_c),
                            by = c("site", "year", "date"))

# Merge the RIEM temperature data with the NEON tick/weather dataset. But only
# use the RIEM temperature data if there is no NEON available...
neon_riem_with_gaps <- full_join(
  x = neon_with_gaps,
  y = full_riem,
  by = c("site" = "site_id", "year", "date"),
  # Columns with identical names will be given suffixes based on which dataset
  # they came from:
  suffix = c("_neon", "_riem")) %>%
  # Use the RIEM-sourced data when NEON is NA:
  mutate(
    mean_temp = case_when(
      # If NEON is NA, then use RIEM
      is.na(mean_temp_c_neon) ~ mean_temp_c_riem,
      # In all other cases use NEON
      TRUE ~ mean_temp_c_neon),
    min_temp = case_when(
      is.na(min_temp_c_neon) ~ min_temp_c_riem,
      TRUE ~ min_temp_c_neon),
    max_temp = case_when(
      is.na(max_temp_c_neon) ~ max_temp_c_riem,
      TRUE ~ max_temp_c_neon),
    
    # Round to two decimal places
    across(.cols = c("mean_temp", "min_temp", "max_temp"),
           .fns = ~round(x = ., digits = 2)),
    # Record which source the temp data came from
    temp_source = case_when(
      
      # Condition: All NEON columns are NA ...AND...
      (is.na(mean_temp_c_neon) & is.na(min_temp_c_neon) & is.na(max_temp_c_neon)) &
        # ...none of the RIEM columns are NA
        (!is.na(mean_temp_c_riem) & !is.na(min_temp_c_riem) & !is.na(max_temp_c_riem)) ~ "RIEM",
      
      # Condition: All of the NEON columns are NA ...AND...
      (is.na(mean_temp_c_neon) & is.na(min_temp_c_neon) & is.na(max_temp_c_neon)) &
        # ...all of the RIEM columns are NA as well
        (is.na(mean_temp_c_riem) & is.na(min_temp_c_riem) & is.na(max_temp_c_riem)) ~ "All temp. sources NA",
      
      # Condition: None of the NEON columns are NA
      !is.na(mean_temp_c_neon) & !is.na(min_temp_c_neon) & !is.na(max_temp_c_neon) ~ "NEON",
      
      # If one of the above options doesn't apply...
      TRUE ~ "Unexpected values")) %>%
  select(site_id = site, year, date, mean_temp, min_temp, max_temp, temp_source) %>%
  mutate(site_id = toupper(site_id))

# Remaining gaps as counts per site and year:
neon_riem_with_gaps %>%
  group_by(site_id, year) %>%
  summarise(across(.cols = everything(), .fns = ~(sum(is.na(.))))) %>%
  data.frame()


# GridMet -----------------------------------------------------------------

# Is GridMet a viable source for gap filling? Compare the correlations between
# the three temperature datasets

gridmet_temp <- read_csv(file = "data/gridmet_temperature_data_compiled.csv")

all_temp_data <- inner_join(x = neon_riem_with_gaps,
                            y = gridmet_temp %>%
                              rename(gridmet_min_temp_c = min_temp_c,
                                     gridmet_max_temp_c = max_temp_c),
                            by = c("site_id" = "field_site_id",
                                   "date"))
# What's the breakdown?
all_temp_data %>%
  count(temp_source)

# What about those unexpected values? They're Inf:
all_temp_data %>%
  filter(temp_source == "Unexpected values")

# To me the NEON relationship looks solid enough for gap filling
all_temp_data %>%
  filter(temp_source %in% c("NEON", "RIEM")) %>%
  ggplot() +
  geom_hex(aes(x = min_temp, y = gridmet_min_temp_c),
           binwidth = 2, color = "gray40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(vars(temp_source)) +
  ylab("GridMet Min. Temp. (C)") +
  xlab("Min. Temp. (C)") +
  xlim(c(-30, 30)) +
  ylim(c(-30, 30)) +
  scale_fill_viridis_c("Count") +
  theme_bw()

all_temp_data %>%
  filter(temp_source %in% c("NEON", "RIEM")) %>%
  ggplot() +
  geom_hex(aes(x = max_temp, y = gridmet_max_temp_c),
           binwidth = 2, color = "gray40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(vars(temp_source)) +
  ylab("GridMet Max. Temp. (C)") +
  xlab("Max. Temp. (C)") +
  xlim(c(-68, 68)) +
  ylim(c(-68, 68)) +
  scale_fill_viridis_c("Count") +
  theme_bw()


# So...if all previous sources were NA...fill them!
# NOTE that I am mostly ignoring mean_temp here, but something could be done
# about it in the future if it's needed
filled_temperature_dataset <- full_join(
  x = neon_riem_with_gaps,
  y = gridmet_temp %>%
    rename(gridmet_min_temp_c = min_temp_c,
           gridmet_max_temp_c = max_temp_c),
  by = c("site_id" = "field_site_id", "date")) %>%
  # Use the RIEM-sourced data when NEON is NA:
  mutate(
    
    new_min_temp = case_when(
      # If min is NA/Inf, then use GridMet
      is.na(min_temp) | is.infinite(min_temp) ~ gridmet_min_temp_c,
      # In all other cases use existing
      TRUE ~ min_temp),
    
    new_max_temp = case_when(
      is.na(max_temp) | is.infinite(max_temp) ~ gridmet_max_temp_c,
      TRUE ~ max_temp),
    
    # Round to two decimal places
    across(.cols = c("mean_temp", "min_temp", "max_temp"),
           .fns = ~round(x = ., digits = 2)),
    
    # Record which source the temp data came from
    temp_source = case_when(
      
      # Condition: Combined columns are NA ...AND...
      (is.na(mean_temp) & is.na(min_temp) & is.na(max_temp)) &
        # ...none of the GridMet columns are NA
        (!is.na(gridmet_min_temp_c) & !is.na(gridmet_max_temp_c)) ~ "GridMet",
      
      # Condition: All of the combined columns are NA ...AND...
      (is.na(mean_temp) & is.na(min_temp) & is.na(max_temp)) &
        # ...all of the GridMet columns are NA as well
        (is.na(gridmet_min_temp_c) & is.na(gridmet_max_temp_c)) ~ "All temp. sources NA",
      
      # Condition: None of the combined columns are NA
      !is.na(mean_temp) & !is.na(min_temp) & !is.na(max_temp) ~ temp_source,
      
      # If one of the above options doesn't apply...
      TRUE ~ "Unexpected values")) %>%
  select(site_id, year, date, mean_temp, min_temp = new_min_temp,
         max_temp = new_max_temp, temp_source)

# Check for duplicated dates
filled_temperature_dataset %>%
  count(site_id, year, date) %>%
  arrange(desc(n))

# None found.

# What's up with the NA years?
filled_temperature_dataset %>%
  filter(is.na(year))

filled_temperature_dataset %>%
  filter(is.na(year)) %>%
  pull(temp_source) %>%
  unique()

# It looks like they occur where GridMet was added

# Any weird temp_sources?
filled_temperature_dataset %>%
  count(temp_source)
# The "All temp. sources NA" are only for 2022

# Are there other NAs remaining?
filled_temperature_dataset %>%
  filter(year < 2022) %>%
  group_by(site_id, year) %>%
  summarise(across(.cols = everything(), .fns = ~(sum(is.na(.))))) %>%
  data.frame()

# NAs only remain in mean_temp! We're good.

filled_temperature_dataset_fix <- filled_temperature_dataset %>%
  mutate(year = year(date)) %>%
  filter(year < 2022)


# 3. Calculate DDs --------------------------------------------------------

# The base temperature to use for DD calculation:
base_temp <- 0

dd_fractions <- filled_temperature_dataset_fix %>%
  arrange(site_id, date) %>%
  group_by(site_id, year) %>%
  rowwise() %>%
  # Using equation 1 from Bouzek et al. 2013
  mutate(dd_mean_temp = mean(c(min_temp, max_temp)),
         dd = if_else(condition = dd_mean_temp < 0,
                      true = 0,
                      false = dd_mean_temp - base_temp)) %>%
  ungroup()

# Now make two new columns:
# thirty_day_dd: The 30-day sum of DDs
# lag_thirty_day_dd: Same as above, but lagged one day
# cume_dd: The cumulative sum of DDs since Jan 01
dd_aggregations <- dd_fractions %>%
  arrange(site_id, date) %>%
  group_by(site_id) %>%
  mutate(thirty_day_dd = rollsum(x = dd,
                                 # 30 days
                                 k = 30,
                                 # Indices with too few data points = NA
                                 fill = NA,
                                 # Move left to right
                                 align = "right"),
         lag_thirty_day_dd = lag(x = thirty_day_dd, n = 1L),
         lag_thirty_day_dd_34wk = lag(x = thirty_day_dd, n = 34L),
         lag_thirty_day_dd_50wk = lag(x = thirty_day_dd, n = 50L),
         lag_thirty_day_dd_42wk = lag(x = thirty_day_dd, n = 42L)) %>%
  group_by(site_id, year) %>%
  mutate(cume_dd = cumsum(dd)) %>%
  ungroup()

# Do the timelines for the 30-day DD and cume DD make sense?
# Red = Cumulative DDs annually
# Blue = Lagged 30-day rolling sum of DDs
dd_aggregations %>%
  filter(site_id == "BLAN") %>%
  ggplot() +
  geom_point(aes(x = date, y = lag_thirty_day_dd, color = "#57106e")) +
  geom_point(aes(x = date, y = cume_dd, color = "#f98e09")) +
  xlab("Date") +
  ylab("Summed DDs") +    
  scale_color_identity(name = "Legend",
                       breaks = c("#f98e09", "#57106e"),
                       labels = c("Cumulative annual", "30-day rolling sum (lag)"),
                       guide = "legend") +
  ggtitle("BLAN degree day aggregation example") +
  theme_bw()

dd_export <- dd_aggregations %>%
  select(-c(mean_temp, min_temp, max_temp, dd_mean_temp))

write_csv(x = dd_export,
          file = "data/daily_dd_calcs.csv")

