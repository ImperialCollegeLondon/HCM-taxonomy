# Predict DDRtree Singapore dataset
#
# Author: Paolo Inglese, Imperial College London, 2023

library(magrittr)
library(tidyverse)
library(here)
library(glue)
library(cowplot)
library(ggplot2)
library(reticulate)
library(FNN)
library(caret)

source(here::here("functions.R"))

data_root <- here::here("../data")
results_root <- here::here("../results")
plots_root <- here::here("../plots")

frames <- c("ED", "ES")

# Import Scikit for trustworthiness function
use_condaenv(condaenv = "hcm")
scikit <- import("sklearn")
py_trustworthiness <- scikit['manifold']$trustworthiness

# Functions ====================================================================

# This function plots the min, mean, median and max of WT to compare
# RBH and Singapore cohorts
summary_wt_plot <- function(wt_mat, dataset, is_plp) {
  lapply(c("min", "mean", "median", "max"), function(f) {
    fun <- switch(
      f,
      min = matrixStats::rowMins,
      mean = rowMeans,
      median = matrixStats::rowMedians,
      max = matrixStats::rowMaxs
    )
    tibble(y = fun(wt_mat), is_plp = factor(is_plp), dataset = dataset) %>%
      ggplot(aes(x = dataset, y = y, color = is_plp)) +
      ggbeeswarm::geom_beeswarm() +
      ggtitle(glue("{f} WT from mesh"))
  })
}

strip_values <- function(x) {
  substr(x, 1, 9)
}

# ggplot2 params for tree plotting
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
    plot.title = element_text(size = unit(8, "pt"), face = "bold"),
    legend.title = element_text(size = unit(6, "pt")),
    legend.key.size = unit(6, "pt"),
    legend.text = element_text(size = unit(6, "pt")),
    legend.box.margin = unit(c(0, 0, 0, 0), "mm"),
    plot.margin = unit(c(0, 0, 0, 0), "mm")
  )
)


# Load rbh data ================================================================

metadata_rbh <- load_metadata(data_root) %>%
  mutate(is_plp = as.numeric(genotype == "P/LP"),
         is_deceased = factor(as.numeric(Deceased)),
         sex = factor(sex),
         race = factor(race),
         sbp = as.numeric(Gen.Systolic)) %>%
  dplyr::rename(age = age_at_scan,
                weight = Gen.Weight,
                height = Gen.Height)

tree_data <- lapply(c("ED", "ES"), function(fr) {
  # Load DDRtree
  ddr_path <- file.path(results_root, glue("ddrtree_res_{fr}_decim0.99.csv"))
  if (!file.exists(ddr_path)) return(NULL)
  ddrtree <- read_csv(ddr_path) %>%
    mutate(across(c(cluster, branch, branch_merged), ~ as.factor(.x))) %>%
    mutate(branch = branch_merged)
  if (fr == "ES") {
    ddrtree$Z1 <- -ddrtree$Z1
    ddrtree$Z2 <- -ddrtree$Z2
    ddrtree <- ddrtree %>%
      mutate(branch = plyr::mapvalues(.$branch,
                                      c("4", "5"),
                                      c("5", "4"))) %>%
      mutate(branch = fct_relevel(branch, "1", "2", "3", "4", "5"))
  }
  ddrtree %>%
    select(ID, Z1, Z2, branch)
}) %>%
  set_names(c("ED", "ES"))

# Join RBH with tree coordinates

rbh_data <- lapply(frames, function(fr) {
  load_wallthickness(data_root = data_root, frame = fr,
                     load_curvature = FALSE) %>%
    inner_join(metadata_rbh %>%
                 select(ID, age, race, sex)) %>%
    mutate(race = factor(ifelse(race == "NFE", "NFE", "non-NFE"))) %>%
    inner_join(tree_data[[fr]], by = "ID")
}) %>%
  set_names(frames)

# Load Singapore data ==========================================================

metadata_sing <- read_csv(file.path(data_root, "singapore_metadata.csv")) %>%
  dplyr::rename(age = age_at_time_of_scan,
                race = Ethnicity,
                sex = Gender) %>%
  mutate(height = height * 100) %>%
  mutate(race = "non-NFE") %>%
  select(ID, age, race, sex)

sing_data <- lapply(frames, function(fr) {
  load_wallthickness(data_root = data_root, frame = fr,
                     dataset = "singapore", load_curvature = FALSE) %>%
    mutate_at(vars("ID"), strip_values) %>%
    inner_join(metadata_sing)
}) %>%
  set_names(frames)

# Combine all data =============================================================

data_all <- lapply(frames, function(fr) {
  bind_rows(
  rbh_data[[fr]] %>%
    mutate(dataset = "RBH"),
  sing_data[[fr]] %>%
    mutate(dataset = "Sing")
)  %>%
  select(ID, age, sex, race, starts_with("wt"), dataset)
}) %>%
  set_names(frames)

# Adjust wall thicknesses ======================================================

# WT are adjusted by residualization from a linear model using age at scan
# and sex as covariates. The adjustement is performed within cohort, separately.

# RBH

wt_rbh_adj <- lapply(frames, function(fr) {
  wt_orig <- rbh_data[[fr]] %>%
    select(starts_with(c("wt", "curv"))) %>%
    as.matrix() %>%
    set_colnames(paste0("wt", seq(1, ncol(.)))) %>%
    set_rownames(rbh_data$ID)
  
  n_feats <- ncol(wt_orig)
  
  # Adjust by residualisation
  covariates <- data_all[[fr]] %>%
    filter(dataset == "RBH") %>%
    select(age, sex)
  
  wt_adj <- matrix(NA, sum(data_all[[fr]]$dataset == "RBH"), n_feats,
                   dimnames = dimnames(wt_orig))
  for (i in seq_len(n_feats)) {
    wt_adj[, i] <- data.frame(y = scale(wt_orig[, i]), covariates) %>%
      { resid(
        lm(y ~ -1 + age + sex, data = .)
      ) }
  }
  return(wt_adj)
}) %>%
  set_names(frames)

# Singapore

wt_sing_adj <- lapply(frames, function(fr) {
  wt_orig <- sing_data[[fr]] %>%
    select(starts_with(c("wt", "curv"))) %>%
    as.matrix() %>%
    set_colnames(paste0("wt", seq(1, ncol(.)))) %>%
    set_rownames(rbh_data$ID)
  
  n_feats <- ncol(wt_orig)
  
  # Adjust by residualisation
  covariates <- data_all[[fr]] %>%
    filter(dataset == "Sing") %>%
    select(age, sex)
  
  wt_adj <- matrix(NA, sum(data_all[[fr]]$dataset == "Sing"), n_feats,
                   dimnames = dimnames(wt_orig))
  for (i in seq_len(n_feats)) {
    wt_adj[, i] <- data.frame(y = scale(wt_orig[, i]), covariates) %>%
      { resid(
        lm(y ~ -1 + age + sex, data = .)
      ) }
  }
  
  return(wt_adj)
}) %>%
  set_names(frames)

lapply(frames, function(fr) {
  plot_adj <- summary_wt_plot(rbind(wt_rbh_adj[[fr]], wt_sing_adj[[fr]]),
                              dataset = data_all[[fr]]$dataset,
                              is_plp = 1)
  plot_grid(plotlist = plot_adj)
}) %>%
  {
    plot_grid(plotlist = .)
  }

# Full dataset adjusted

data_adj <- lapply(frames, function(fr) {
  cbind(
    rbind(wt_rbh_adj[[fr]], wt_sing_adj[[fr]]),
    data_all[[fr]] %>% select(-starts_with("wt"))
  )
}) %>%
  set_names(frames)

# Show that the adjustment fixed the batch in PCA space
lapply(frames, function(fr) {
  
  tmp <- data_adj[[fr]] %>%
    cbind(., irlba::prcomp_irlba(
      data_adj[[fr]] %>%
        select(starts_with("wt")), n = 2, scale. = TRUE)$x)
  
  tmp %>%
    ggplot(aes(x = PC1, y = PC2, color = dataset)) +
    geom_point(size = 3) +
    scale_color_brewer(palette = "Dark2") +
    ggtitle(fr)
  
})

# PCA ==========================================================================

# Calculate first 5 PC from RBH and project Singapore onto them. These are
# used as features to predict the tree coordinates.

pca <- lapply(frames, function(fr) {
  data_adj[[fr]] %>%
  filter(dataset == "RBH") %>%
  select(starts_with("wt")) %>%
  {
    irlba::prcomp_irlba(., n = 5, center = TRUE, scale. = TRUE)
  }
}) %>%
  set_names(frames)

scores_rbh <- lapply(frames, function(fr) {
  data_adj[[fr]] %>%
  filter(dataset == "RBH") %>%
  select(starts_with("wt")) %>%
  {
    predict(pca[[fr]], .)
  }
}) %>%
  set_names(frames)

scores_sing <- lapply(frames, function(fr) {
  data_adj[[fr]] %>%
  filter(dataset == "Sing") %>%
  select(starts_with("wt")) %>%
  {
    predict(pca[[fr]], .)
  }
}) %>%
  set_names(frames)

# Predict coords ===============================================================

branch_sing <- plot_projection <- plot_corr <- tw <- vector("list", 2) %>% 
  set_names(frames)

for (fr in frames) {
  
  message(glue("Frame = {fr}"))
  
  coords_rbh <- rbh_data[[fr]] %>%
    inner_join(data_all[[fr]], by = "ID") %>%
    select(Z1, Z2)
  
  message("Fit X model")
  
  model_x <- train(x = scores_rbh[[fr]],
                   y = coords_rbh$Z1,
                   method = "rf",
                   trControl = trainControl(method = "repeatedcv", number = 10,
                                            repeats = 3))
  print(model_x)
  
  message("Fit Y model")
  
  model_y <- train(x = scores_rbh[[fr]],
                   y = coords_rbh$Z2,
                   method = "rf",
                   trControl = trainControl(method = "repeatedcv", number = 10,
                                            repeats = 3))
  print(model_y)
  
  coords_sing_x <- predict(model_x, as.data.frame(scores_sing[[fr]]))
  coords_sing_y <- predict(model_y, as.data.frame(scores_sing[[fr]]))
  
  # Find NN in RBH
  coords_rbh <- rbh_data[[fr]] %>% select(ID, Z1, Z2) %>%
    inner_join(data_adj[[fr]] %>% select(ID), by = "ID") %>%
    select(Z1, Z2)
  
  coords_sing <- cbind(coords_sing_x, coords_sing_y)
  
  nn_index <- get.knnx(coords_rbh,
                       coords_sing,
                       k = 1)$nn.index
  
  # Assing Singapore data to branches
  branch_sing[[fr]] <- (
    rbh_data[[fr]] %>% select(ID, branch) %>%
      inner_join(data_adj[[fr]] %>% select(ID), by = "ID") %>%
      pull(branch)
  )[nn_index]
  
  plot_projection[[fr]] <- rbh_data[[fr]] %>%
    mutate(set = ifelse(ID %in% data_adj[[fr]]$ID, "RBH", "RBH_not_used")) %>%
    select(Z1, Z2, set) %>%
    bind_rows(
      tibble(Z1 = coords_sing[, 1],
             Z2 = coords_sing[, 2],
             set = "Sing.")
    ) %>%
    ggplot(aes(x = Z1, y = Z2, color = set, shape = set)) +
    geom_point(size = 3) +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c("RBH" = 21, "Sing." = 19)) +
    tree_plot_params +
    ggtitle(glue("{fr} - Projected Singaporean data"))
  
  print("Assigned branch composition:")
  print(table(branch_sing[[fr]]))
  
  # Trustworthiness
  tw[[fr]] <- sapply(as.integer(seq(1, 10)), function(k) {
    py_trustworthiness(
      wt_sing_adj[[fr]],
      coords_rbh[nn_index, ],
      n_neighbors = k
    )
  })
  
  # Calculate Spearman correlation between nearest neighbours
  
  # Singapore
  
  data_rbh_matched <- data_adj[[fr]] %>%
    filter(dataset == "RBH") %>%
    filter(ID %in% data_adj[[fr]]$ID)
  
  mat_sing <- (data_adj[[fr]] %>%
                 filter(dataset == "Sing") %>%
                 select(starts_with("wt")) %>%
                 as.matrix() %>%
                 scale())
  mat_rbh <- (data_rbh_matched %>%
                select(starts_with("wt")) %>%
                as.matrix() %>%
                scale())
  
  corr_sing_rbh <- array(NA, nrow(coords_sing))
  for (i in 1:nrow(coords_sing)) {
    x <- mat_sing[i, ]
    y <- mat_rbh[nn_index[i], ]
    corr_sing_rbh[i] <- cor(x, y, method = "spearman")
  }
  
  # RBH
  
  coords_rbh_all <- rbh_data[[fr]] %>% select(Z1, Z2)
  
  # Nearest neighbours of all RBH
  nn_index_rbh <- get.knn(coords_rbh_all, k = 1)$nn.index
  
  mat_rbh <- (data_adj[[fr]] %>% filter(dataset == "RBH") %>%
                select(starts_with("wt")) %>%
                as.matrix() %>%
                scale())
  
  corr_rbh_rbh <- array(NA, nrow(coords_rbh_all))
  for (i in 1:nrow(coords_rbh_all)) {
    x <- mat_rbh[i, ]
    y <- mat_rbh[nn_index_rbh[i], ]
    corr_rbh_rbh[i] <- cor(x, y, method = "spearman")
  }
  
  plot_corr[[fr]] <- bind_rows(
    tibble(rho = corr_rbh_rbh, dataset = "RBH"),
    tibble(rho = corr_sing_rbh, dataset = "Sing.")
  ) %>%
    ggplot(aes(x = dataset, y = rho, color = dataset)) +
    geom_boxplot(width = 0.1, color = "black") +
    ggbeeswarm::geom_beeswarm() +
    xlab("Set") +
    ylab("Spearman's rho") +
    ggtitle(fr) +
    scale_color_brewer(palette = "Dark2") +
    theme_classic()
  
  print("Wilcoxon's test:")
  print(wilcox.test(corr_rbh_rbh, corr_sing_rbh))
  
  print("KS test:")
  print(ks.test(corr_rbh_rbh, corr_sing_rbh, simulate.p.value = TRUE, B = 1e4))

}

plot_projection_combo <- plot_grid(plotlist = plot_projection, labels = c("A", "B"),
                                   align = "h")

ggsave(filename = file.path(plots_root, "plots_for_figures",
                            glue("plot_tree_sing_projection.pdf")),
       width = 210, height = 100, dpi = 300, plot = plot_projection_combo,
       units = "mm")

plot_corr_combo <- plot_grid(
  plotlist = lapply(plot_corr, function(pl) {
    pl$layers[[2]]$aes_params$alpha = 0.3
    pl
  }), labels = c("A", "B"), align = "h")

ggsave(filename = file.path(plots_root, "plots_for_figures",
                            glue("plot_corr_projection.pdf")),
       width = 210, height = 100, dpi = 300,
       plot = plot_corr_combo,
       units = "mm")
