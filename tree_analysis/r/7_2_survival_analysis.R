# Survival analysis HCM
#
# Author: Paolo Inglese, Imperial College London, 2022

library(survival)
library(lubridate)
library(survminer)
library(here)
library(glue)
library(tidyverse)
library(magrittr)
library(gtsummary)


source(here::here("functions.R"))

# Setup ========================================================================

data_root <- here::here("../data")
plots_root <- here::here("../plots")
results_root <- here::here("../results")

# Load =========================================================================

event_data <- read_csv(file.path(data_root, "survival_event_data.csv")) %>%
  mutate(
    sex = factor(sex),
    race = factor(race),
    is_plp = factor(ifelse(is_plp == 1, "Y", "N"))
  )

# Main =========================================================================

# Full model
cox_model_full <- event_data %>%
  mutate(race = factor(ifelse(race == "NFE", "NFE", "non-NFE"))) %>%
  {
    coxph(Surv(time_censor, is_deceased) ~ is_plp + race + sex,
      data = ., robust = TRUE
    )
  }
print(cox_model_full)

# Test validity model
test_ph <- cox.zph(cox_model_full)
print(test_ph)

ggcoxzph(test_ph)

ggcoxdiagnostics(cox_model_full,
  type = "dfbeta",
  linear.predictions = FALSE, ggtheme = theme_bw()
)

# Model only P/LP
cox_model_plp <- coxph(Surv(time_censor, is_deceased) ~ is_plp,
  data = event_data, robust = TRUE
)

print(cox_model_plp)

# Test validity model
test_ph <- cox.zph(cox_model_plp)
print(test_ph)

ggcoxzph(test_ph)

ggcoxdiagnostics(cox_model_plp,
  type = "dfbeta",
  linear.predictions = FALSE, ggtheme = theme_bw()
)

# Print the model params
latex_printout <- lapply(list(cox_model_full, cox_model_plp), tbl_regression,
                         exp = T, digits = 4) %>%
  tbl_merge(
    tab_spanner = c("**Full model**", "**Genotype only**")
  ) %>%
  as_gt() %>%
  gt::as_latex()

print(glue("{latex_printout}"))

# Plot =========================================================================

fit <- survfit(
  Surv(time_censor, is_deceased) ~ is_plp,
  data = event_data
)

ggsurv <- ggsurvplot(
  fit,
  xlim = c(20, 98),
  ylim = c(0, 1),
  xlab = "Life duration (Years)",
  legend = "bottom",
  legend.title = "Genotype", font.legend = c(14),
  palette = c("#BDBDBD", "#A52A2A"),
  risk.table = TRUE,
  risk.table.title = NULL,
  risk.table.height = 0.20,
  linetype = "strata",
  fontsize = 4, font.x = 14, font.y = 14, font.tickslab = 12,
  break.time.by = 10,
  break.y.by = 0.1,
  title = NULL,
  fun = "cumhaz"
)

pdf(file = file.path(plots_root, "survival_model_paolo.pdf"), width = 8,
    height = 8)
print(ggsurv)
dev.off()
