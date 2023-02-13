# Generate data for survival analysis
#
# Author: Paolo Inglese, Imperial College London, 2023


library(dplyr)
library(tidyr)
library(readr)
library(magrittr)
library(here)
library(lubridate)


source(here::here("functions.R"))

# Setup ========================================================================

data_root <- here::here("../data")

# Functions ====================================================================

convert_date_death <- function(x) {
  na_values <- is.na(x)
  death_dates <- x[!na_values]
  d_ <- array(NA, length(death_dates))
  d_ <- format(as.POSIXct(death_dates, format = "%a %b %d %Y"),
               format = '%Y/%m/%d')
  d_[is.na(d_)] <- format(as.POSIXct(death_dates[is.na(d_)], format = "%d/%m/%Y"),
                          format = "%Y/%m/%d")
  d <- array(NA, length(x))
  d[!na_values] <- d_
  return(d)
}

convert_date_enrol <- function(x) {
  format(as.POSIXct(x, format = "%d/%m/%Y"), format = '%Y/%m/%d')
}

# Main =========================================================================

# This is the metadata of the 436 subjects
metadata <- load_metadata(data_root) %>%
  mutate(is_plp = factor(as.numeric(genotype == "P/LP")),
         is_deceased = factor(as.numeric(Deceased)),
         sex = factor(sex),
         race = factor(race))
metadata_dict <- read_csv(file.path(data_root, "metadata_dict.csv"))

# The original metadata contains the relevant dates
metadata_orig <- read_csv(file.path(data_root, "rbh_superset_filamentinfo_utf8.csv"))

# Create the event data from the whole cohort
events_full_dataset <- metadata_orig %>%
  filter(!((Isdeceased.x) & (is.na(Date.of.death.x)))) %>%
  select(ID, enrolment.date, Dob.Str.x, Date.of.death.x, Date.mortality.checked,
         Isdeceased.x, genotype, age_at_scan) %>%
  mutate(is_plp = ifelse(genotype == "PLP", 1, 0)) %>%
  select(-genotype) %>%
  mutate(Isdeceased.x = ifelse(Isdeceased.x, 1, 0)) %>%
  dplyr::rename(date_end_study = Date.mortality.checked) %>%
  dplyr::rename(date_death = Date.of.death.x) %>%
  dplyr::rename(date_birth = Dob.Str.x) %>%
  dplyr::rename(date_enrol = enrolment.date) %>%
  dplyr::rename(is_deceased = Isdeceased.x) %>%
  filter(!(is_deceased & is.na(date_death))) %>%
  mutate(date_end_study = as.Date(convert_date_enrol(date_end_study)),
         date_event = convert_date_death(date_death),
         date_birth = as.Date(convert_date_enrol(date_birth)),
         date_enrol = as.Date(convert_date_enrol(date_enrol))) %>%
  mutate_at(vars(date_event), ~replace_na(., "2022/03/02")) %>%
  mutate(date_event = as.Date(date_event)) %>%
  mutate(age_at_enrol = interval(ymd(date_birth), ymd(date_enrol)),
         time_censor = interval(ymd(date_birth), ymd(date_event)))

rm(metadata_orig)

# Event data for the 436 subjects
event_data <- events_full_dataset %>%
  inner_join(
    metadata %>% select(ID, sex, race, genotype)
  ) %>%
  mutate(genotype = factor(genotype))

# Save

write_csv(x = event_data, file = file.path(data_root, "survival_event_data.csv"))
