# Prepare data for unsupervised analysis
#
# Author: Paolo Inglese, Imperial College London, 2022


# Libraries ====================================================================

library(mice)
library(here)
library(magrittr)
library(tidyverse)

# Setup ========================================================================

data_root <- here('../data')

# Load data ====================================================================

metadata <- read_csv(
  file.path(data_root, 'rbh_superset_filamentinfo_utf8.csv')) %>%
  filter(genotype %in% c('NEG', 'VUS', 'PLP')) %>%
  filter(!duplicated(ID))

metadata <- metadata %>%
  as_tibble() %>%
  # select(-starts_with('Hcm.Dateof'), -starts_with('Date'),
  #        -site, -platform, -target_id, -group, -starts_with('pav'),
  #        -starts_with('sarch'), -starts_with('fafh'), -plphet, -NEG,
  #        -starts_with('Biobank'), -starts_with('Consent.'), -starts_with('consent'),
  #        -starts_with('Tests.Text'), -starts_with('SAR'), -PHENOCOPY,
  #        -starts_with('Race.'), -starts_with('Sex.'), -starts_with('Dob'),
  #        -starts_with('Dod'), -starts_with('Age.'), -starts_with('Analysed.'),
  #        -enrolment.date, -days.between.MRI.and.enrolment,
  #        -starts_with('Title.'), -ends_with('Date'),
  #        -ends_with('code'), -ends_with('Comment'), -Ethics) %>%
  select(ID,
         age_at_scan,
         filaments,
         genotype,
         race,
         sex,
         type,
         Isdeceased.x,
         Gen.Acearb,
         Gen.Activityscore,
         Gen.Alcohol,
         Gen.Asaclopi,
         Gen.Betablocker,
         Gen.Cad,
         Gen.Cabg,
         Gen.Ccs,
         Gen.Diastolic,
         Gen.Diuretic,
         Gen.Height,
         Gen.Dm,
         Gen.Ht,
         Gen.Mi,
         Gen.Nyha,
         Gen.Pci,
         Gen.Pulserate,
         Gen.Smoking,
         Gen.Systolic,
         Gen.Weight,
         Hcm.Bsa,
         Hcm.Coincidentinfarction,
         Hcm.Edv,
         Hcm.Ef,
         Hcm.Esv,
         Hcm.Familyhistoryofhcm,
         Hcm.Familyhistoryofscd,
         Hcm.Lvm,
         Hcm.Lvgadolinum,
         Hcm.Lvmostaffectedsegment,
         Hcm.Lvoto,
         Hcm.Lvotpeakvelocity,
         Hcm.Maxlvwallthickness,
         Hcm.Mitralregurgitation,
         Hcm.Mostaffectedlevel,
         Hcm.Perfusiondeficit,
         Hcm.Rvhypertrophy,
         Hcm.Sv) %>%
  mutate(across(  # Set 0 to NA for quantitative measures
    c(Gen.Weight,
      Gen.Height,
      Gen.Ht,
      age_at_scan,
      Hcm.Lvm,
      Hcm.Maxlvwallthickness,
      Gen.Pulserate,
      Gen.Systolic,
      Gen.Diastolic,
      Hcm.Bsa,
      Hcm.Edv,
      Hcm.Esv,
      Hcm.Sv,
      Hcm.Ef,
      Hcm.Mitralregurgitation,
      Hcm.Lvmostaffectedsegment,
      Hcm.Mostaffectedlevel), ~na_if(., 0))) %>%
  mutate(Gen.Smoking = as.factor(na_if(Gen.Smoking, 3))) %>%
  mutate(Gen.Smoking = fct_recode(Gen.Smoking, "0"="1", "1"="2", "2"="4")) %>%
  mutate(Hcm.Lvmostaffectedsegment = as.factor(na_if(Hcm.Lvmostaffectedsegment, 0))) %>%
  mutate(Hcm.Lvmostaffectedsegment = fct_recode(Hcm.Lvmostaffectedsegment,
                                                "Septal" = "1",
                                                "Anterior" = "2",
                                                "Lateral" = "3",
                                                "Inferior" = "4")) %>%
  mutate(Hcm.Mostaffectedlevel = as.factor(na_if(Hcm.Mostaffectedlevel, 0))) %>%
  mutate(Hcm.Mostaffectedlevel = fct_recode(Hcm.Mostaffectedlevel,
                                            "Base" = "1",
                                            "Mid" = "2",
                                            "Apex" = "3")) %>%
  mutate(Isdeceased.x = as.logical(Isdeceased.x)) %>%
  mutate(filaments = as.factor(plyr::mapvalues(filaments, NA, "None"))) %>%
  mutate(type = as.factor(plyr::mapvalues(type, NA, "None"))) %>%
  mutate(sex = as.factor(sex),
         race = as.factor(race),
         genotype = as.factor(genotype)) %>%
  mutate(Gen.Ccs = as.integer(Gen.Ccs),
         Gen.Nyha = as.integer(Gen.Nyha),
         Gen.Activityscore = as.integer(Gen.Activityscore)) %>%
  rename(Deceased = Isdeceased.x)

code_123 <- c(
  'Gen.Ht',
  'Gen.Cad',
  'Gen.Mi',
  'Gen.Dm',
  'Gen.Cabg',
  'Gen.Pci',
  'Gen.Asaclopi',
  'Gen.Diuretic',
  'Gen.Betablocker',
  'Gen.Acearb',
  'Hcm.Rvhypertrophy',
  'Hcm.Familyhistoryofhcm',
  'Hcm.Coincidentinfarction',
  'Hcm.Lvoto',
  'Hcm.Perfusiondeficit',
  'Hcm.Familyhistoryofscd')

metadata <- metadata %>%
  mutate(across(code_123, ~as.factor(na_if(.x, 3))))

code_12345 <- c('Hcm.Lvgadolinum', 'Hcm.Mitralregurgitation')

metadata <- metadata %>%
  mutate(across(code_12345, ~as.integer(na_if(.x, 5))))

mice_metadata <- mice(metadata, method = 'pmm')
metadata_imputed <- metadata

for (i in seq_along(colnames(metadata))) {
  cn <- colnames(metadata)[i]
  if (any(is.na(metadata[[cn]]))) {
    imputed_values <- mice_metadata$imp[[cn]]
    if (is.factor(metadata_imputed[[cn]])) {
      imputed_value <- apply(imputed_values, 1, function(w) {
        names(table(w)[which.max(table(w))])
      })
    } else if (is.integer(metadata_imputed[[cn]])) {
      imputed_value <- apply(imputed_values, 1, function(w) {
        as.integer(names(table(w)[which.max(table(w))]))
      })
    } else {
      imputed_value <- rowMeans(imputed_values)
    }
    metadata_imputed[[cn]][
      as.integer(rownames(imputed_values))] <- imputed_value
  }
}

write_csv(metadata_imputed, file = file.path(data_root, 'metadata_imputed.csv'))

# Write the dictionary for the metadata values =================================

metadata_dict <- array("continuous", ncol(metadata)) %>%
  set_names(colnames(metadata))

metadata_dict[names(metadata_dict) == 'ID'] <- 'character'

metadata_dict[
  names(metadata_dict) %in% (metadata %>%
    select_if(~is.factor(.x)) %>%
    colnames())] <- "factor"

metadata_dict[
  names(metadata_dict) %in% (metadata %>%
    select_if(~is.integer(.x)) %>%
    colnames())] <- "discrete"

metadata_dict[
  names(metadata_dict) %in% (metadata %>%
    select_if(~is.logical(.x)) %>%
    colnames())] <- "logical"

metadata_dict %>%
  as_tibble() %>%
  mutate(variable = names(metadata_dict), .before = value) %>%
  write_csv(file = file.path(data_root, 'metadata_dict.csv'))
