# Plot DDRtree coloured by probability of mortality and P/LP
#
# Author: Paolo Inglese, Imperial College London, 2023

library(tidyverse)
library(magrittr)
library(here)
library(ggrepel)
library(latex2exp)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(paletteer)
library(glue)
library(mgcv)
library(ROCR)
library(survival)
library(lubridate)
library(caret)
library(rms)
library(dsm)
library(pec)
library(tidymv)

source(here::here("functions.R"))

tree_plot_params <- list(
  xlab("Z1"),
  ylab("Z2"),
  coord_fixed(),
  theme_classic(),
  theme(
    legend.position = "top",
    text = element_text(size = unit(6, "pt")),
    axis.title = element_text(size = unit(6, "pt")),
    axis.text = element_text(size = unit(6, "pt")),
    plot.title = element_text(size = unit(6, "pt"), face = "bold"),
    legend.title = element_text(size = unit(6, "pt")),
    legend.key.size = unit(6, "pt"),
    legend.text = element_text(size = unit(6, "pt")),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = unit(c(0, 0, 0, 0), "mm")
  )
)

minmax_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Setup ========================================================================

frames <- c("ED", "ES")
data_root <- here::here("../data")
results_root <- here::here("../results")
plots_root <- here::here("../plots")

# Load data ====================================================================

metadata <- load_metadata(data_root) %>%
  mutate(
    is_plp = factor(as.numeric(genotype == "P/LP")),
    is_deceased = factor(as.numeric(Deceased)),
    sex = factor(sex),
    race = factor(race)
  )
metadata_dict <- read_csv(file.path(data_root, "metadata_dict.csv"))

# Load survival data
event_data <- read_csv(file.path(data_root, "survival_event_data.csv"))

# Load PRS
prs_data <- read_csv(file.path(data_root, "PRScs.HCM.MTAG.exc.MYBPC3.RBH.1KG.LD.csv")) %>%
  mutate(ID = FID)

tree_data <- lapply(c("ED", "ES"), function(frame) {
  # Load DDRtree
  ddrtree <- read_csv(
    file.path(results_root, glue("ddrtree_res_{frame}_decim0.99.csv"))
  ) %>%
    mutate(across(c(cluster, branch, branch_merged), ~ as.factor(.x))) %>%
    mutate(branch = branch_merged)
  ptime <- read_csv(
    file.path(results_root,glue("ddrtree_pseudotime_{frame}_decim_0.99.csv")))
  ddrtree <- ddrtree %>%
    inner_join(ptime)
  if (frame == "ES") {
    # Invert the Z coordinates to match the orientation of ED. This doesn't
    # affect the results, as the tree is parity invariant
    ddrtree$Z1 <- -ddrtree$Z1
    ddrtree$Z2 <- -ddrtree$Z2
    ddrtree$Y1 <- -ddrtree$Y1
    ddrtree$Y2 <- -ddrtree$Y2
    ddrtree <- ddrtree %>%
      mutate(branch = plyr::mapvalues(
        .$branch,
        c("4", "5"),
        c("5", "4")
      )) %>%
      mutate(
        branch = fct_relevel(branch, "1", "2", "3", "4", "5"),
        # Fix ordering from base to tip to calculate survival along branch
        pseudotime = max(pseudotime) - pseudotime
      )
  }
  ddrtree %>%
    select(ID, Z1, Z2, Y1, Y2, branch, pseudotime) %>%
    inner_join(event_data, by = "ID")
}) %>%
  set_names(c("ED", "ES"))

# PRS ==========================================================================

preds_prs <- lapply(frames, function(fr) {
  
  message(glue("Frame = {fr}"))
  
  model_data <- tree_data[[fr]] %>%
    inner_join(prs_data, by = "ID") %>%
    mutate(Y = SCORE) %>%
    as_tibble()
  
  bin_score <- (model_data$SCORE > median(model_data$SCORE)) * 1
  model_data <- model_data %>%
    mutate(score_bin = bin_score)
  
  model_data$score_bin <- factor(model_data$score_bin)
  model_data$score_bin = dplyr::recode_factor(model_data$score_bin, `0` = "lo", `1` = "hi")
  
  model_data <- cbind(model_data, prcomp(model_data %>% select(Z1, Z2) %>% scale())$x)
  
  # Fit GAM to predict if P/LP from the tree coordinates
  set.seed(123)
  
  gam_model <- train(
    x = model_data %>%
      select(Z1, Z2),
    y = model_data$score_bin,
    method = "gam",
    trControl = trainControl(
      method = "repeatedcv",
      number = 10, repeats = 3,
      summaryFunction = twoClassSummary,
      classProbs = TRUE
    ),
    metric = "ROC",
  )
  
  y_pred = predict(gam_model, type = "prob")
  tibble(ID = model_data$ID, prob = y_pred)
}) %>%
  set_names(frames)

# P/LP =========================================================================

# Predict the P/LP using the tree coordinates.

preds_plp <- lapply(frames, function(fr) {
  model_data <- tree_data[[fr]] %>%
    mutate(Y = factor(ifelse(is_plp == 1, "PLP", "non.PLP"))) %>%
    as_tibble()

  # Fit GAM to predict if P/LP from the tree coordinates
  gam_model <- train(
    x = model_data %>%
      select(Z1, Z2),
    y = model_data$Y,
    method = "gam",
    trControl = trainControl(
      method = "repeatedcv",
      number = 10, repeats = 3,
      summaryFunction = twoClassSummary,
      classProbs = TRUE
    ),
    metric = "ROC"
  )
  
  return(list(
    mdl = gam_model,
    preds = predict(gam_model, model_data, type = "prob")[, 2]))
  
}) %>%
  set_names(frames)

saveRDS(preds_plp, file = file.path(results_root, "predictions_tree_plp.rds"))

gam_latex_printout <- gtsummary::tbl_merge(
  list(
    ED = gtsummary::tbl_regression(
      preds_plp$ED$mdl$finalModel,
      intercept = TRUE
    ),
    ES = gtsummary::tbl_regression(
      preds_plp$ES$mdl$finalModel,
      intercept = TRUE
    )
  ),
  tab_spanner = c("**ED**", "**ES**")
) %>%
  gtsummary::as_gt() %>%
  gt::as_latex() %>%
  as.character()

print(glue("{gam_latex_printout}"))

# Odds ratio P/LP branch 4 and 1 ===============================================

OR_obs_b4_b1 <- (tree_data[["ES"]] %>%
  mutate(probs = preds_plp[["ES"]]$preds) %>%
  filter(branch == 4) %>%
  pull(probs) %>%
  median()) /
  (tree_data[["ES"]] %>%
    mutate(probs = preds_plp[["ES"]]$preds) %>%
    filter(branch == 1) %>%
    pull(probs) %>%
    median())

print(OR_obs_b4_b1)

# Permutation test
N_PERM <- 1e4
perm_res <- numeric(N_PERM)

perm_data <- tree_data$ES %>%
  mutate(probs = preds_plp$ES$preds) %>%
  filter(branch %in% c(1, 4))

pb <- progress::progress_bar$new(total = N_PERM, width = 60)

for (i in seq_len(N_PERM)) {
  tmp_data <- perm_data
  tmp_data$branch <- sample(perm_data$branch, replace = FALSE)

  perm_res[i] <- (tmp_data %>%
    filter(branch == 4) %>%
    pull(probs) %>%
    median()) /
    (tmp_data %>%
      filter(branch == 1) %>%
      pull(probs) %>%
      median())

  pb$tick()
}
pb$terminate()

perm_pval <- (sum(perm_res >= OR_obs_b4_b1) + 1) / (N_PERM + 1)
print(glue("P-value = {round(perm_pval, 4)}"))

# Bootstrapping

pb <- progress::progress_bar$new(total = N_PERM, width = 60)
boot_res <- numeric(N_PERM)

for (i in seq_len(N_PERM)) {
  tmp_data <- perm_data[sample(seq_len(nrow(perm_data)), replace = TRUE), ]

  boot_res[i] <- (tmp_data %>%
    filter(branch == 4) %>%
    pull(probs) %>%
    median()) /
    (tmp_data %>%
      filter(branch == 1) %>%
      pull(probs) %>%
      median())

  pb$tick()
}
pb$terminate()

print(
  glue(
    "Bootstrapping 95% CI = {round(quantile(boot_res, c(0.025, 0.975)), 4)}"))

# Odds ratio P/LP branch 4 =====================================================

model_data <- tree_data[["ES"]] %>%
  mutate(probs = preds_plp[["ES"]]$preds) %>%
  filter(branch == 4) %>%
  mutate(pseudotime = minmax_scale(pseudotime)) %>%
  arrange(pseudotime)

(model_data %>%
  filter(pseudotime == max(pseudotime)) %>%
  pull(probs) /
  model_data %>%
    filter(pseudotime == min(pseudotime)) %>%
    pull(probs))

# Mortality ====================================================================

# Use Cox PH with smooth cubic splines of tree coordinates + tree ordering
# (called pseudotime) to estimate survival probability

# Use median observed age as prediction time point
med_surv <- tree_data$ES %>%
  pull(time_censor) %>%
  median() %>%
  round()

preds_surv <- lapply(frames, function(fr) {
  mdl <- cph(Surv(time_censor, is_deceased) ~ pseudotime + rcs(Z1, 3) +
               rcs(Z2, 3),
    data = tree_data[[fr]],
    x = TRUE,
    y = TRUE,
    surv = TRUE,
    iter.max = 100,
    robust = FALSE
  )

  prob <- predictSurvProb(mdl, newdata = tree_data[[fr]], times = med_surv)

  tibble(
    ID = tree_data[[fr]]$ID,
    prob = prob
  )
}) %>%
  set_names(frames)

OR <- lapply(frames, function(fr) {
  median(preds_surv[[fr]]$prob[tree_data[[fr]]$branch == 4]) /
    median(preds_surv[[fr]]$prob[tree_data[[fr]]$branch == 1])
}) %>%
  set_names(frames)

# Plots ========================================================================

plot_tree_hi_prs_probs <- lapply(frames, function(fr) {
  model_data <- tree_data[[fr]] %>%
    inner_join(preds_prs[[fr]])
  
  model_data %>%
    ggplot(aes(x = Z1, y = Z2, color = prob)) +
    geom_point(size = 2, alpha = 0.75) +
    ggtitle(glue("Prob. high PRS - {fr}")) +
    tree_plot_params +
    theme(
      legend.key.size = unit(6, "pt"),
      legend.text = element_text(size = unit(6, "pt")),
      legend.title = element_text(size = unit(8, "pt")),
      legend.key.width = unit(0.75, "cm")
    ) +
    scale_color_viridis_c(option = "magma") +
    coord_fixed()
}) %>%
  set_names(frames)

plot_tree_plp_probs <- lapply(frames, function(fr) {
  model_data <- tree_data[[fr]] %>%
    mutate(prob = preds_plp[[fr]]$preds)

  model_data %>%
    ggplot(aes(x = Z1, y = Z2, color = prob)) +
    geom_point(size = 2, alpha = 0.75) +
    ggtitle(glue("Probability P/LP - {fr}")) +
    tree_plot_params +
    theme(
      legend.key.size = unit(6, "pt"),
      legend.text = element_text(size = unit(6, "pt")),
      legend.title = element_text(size = unit(8, "pt")),
      legend.key.width = unit(0.75, "cm")
    ) +
    coord_fixed()
}) %>%
  set_names(frames)

plot_tree_survival_probs <- lapply(frames, function(fr) {
  model_data <- tree_data[[fr]] %>%
    inner_join(preds_surv[[fr]])

  model_data %>%
    ggplot(aes(x = Z1, y = Z2, color = prob)) +
    geom_point(size = 2, alpha = 0.75) +
    ggtitle(glue("Prob. Survival ({med_surv} years) - {fr}")) +
    tree_plot_params +
    theme(
      legend.key.size = unit(6, "pt"),
      legend.text = element_text(size = unit(6, "pt")),
      legend.title = element_text(size = unit(8, "pt")),
      legend.key.width = unit(0.75, "cm")
    ) +
    coord_fixed()
}) %>%
  set_names(frames)

plot_branch_4 <- lapply(frames, function(fr) {
  model_data <- tree_data[[fr]] %>%
    inner_join(preds_surv[[fr]]) %>%
    filter(branch == 4) %>%
    mutate(pseudotime = minmax_scale(pseudotime)) %>%
    arrange(pseudotime)

  (model_data %>%
    filter(pseudotime == max(pseudotime)) %>%
    pull(prob) /
    model_data %>%
      filter(pseudotime == min(pseudotime)) %>%
      pull(prob))


  gam_model <- gam(prob ~ pseudotime + s(pseudotime, bs = "cs", k = 50),
    data = model_data, method = "REML", select = TRUE
  )

  preds <- predict(gam_model, se.fit = TRUE) %>%
    as_tibble()

  cbind(model_data, preds) %>%
    ggplot(aes(pseudotime, fit)) +
    geom_smooth_ci(color = "red") +
    geom_point(aes(x = pseudotime, y = prob),
      data = model_data, inherit.aes = F,
      size = 1, alpha = 0.5
    ) +
    theme_classic() +
    xlab("Position in branch 4") +
    ylab(glue("Survival probability ({med_surv} yrs)")) +
    ggtitle(glue("Frame {fr}")) +
    theme(text = element_text(family = "sans"))
})

ggsave(
  filename = file.path(plots_root, "surv_prob_branch_4_ES.pdf"),
  width = 105, height = 60, units = "mm", plot = plot_branch_4[[2]]
)

# Apply common scale to P/LP and mortality: only ES

range_prob_plp <- preds_plp[[2]] %>%
  {
    c(min(round(min(.$preds), 4) - 1e-4), max(round(max(.$preds), 4) + 1e-4))
  }

range_prob_survival <- preds_surv[[2]] %>%
  {
    c(min(round(min(.$prob), 4) - 1e-4), max(round(max(.$prob), 4) + 1e-4))
  }

plot_tree_plp_probs <- lapply(plot_tree_plp_probs, function(p) {
  p +
    scale_color_paletteer_c(
      "grDevices::terrain.colors",
      name = NULL,
      limits = range_prob_plp,
      breaks =
        seq(range_prob_plp[1], range_prob_plp[2], length.out = 4),
      labels = round(
        seq(range_prob_plp[1], range_prob_plp[2], length.out = 4), 4
      )
    )
})

plot_tree_survival_probs <- lapply(plot_tree_survival_probs, function(p) {
  p +
    scale_color_paletteer_c(
      "grDevices::ag_GrnYl",
      name = NULL,
      limits = range_prob_survival,
      breaks = seq(
        range_prob_survival[1], range_prob_survival[2],
        length.out = 4
      ),
      labels = round(
        seq(
          range_prob_survival[1], range_prob_survival[2],
          length.out = 4
        ),
        4
      )
    )
})

plot_ES_tree <- plot_grid(
  plotlist = list(plot_tree_plp_probs$ES, plot_tree_survival_probs$ES),
  align = "h", labels = c("A", "B")
)
ggsave(
  filename = file.path(plots_root, "tree_plots_ES_probs.pdf"),
  width = 210, height = 70, dpi = 300, plot = plot_ES_tree, units = "mm"
)
