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
    plot.title = element_text(size = unit(8, "pt"), face = "bold"),
    legend.title = element_text(size = unit(6, "pt")),
    legend.key.size = unit(6, "pt"),
    legend.text = element_text(size = unit(6, "pt")),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = unit(c(0, 0, 0, 0), "mm")
  )
)

# Setup ========================================================================

data_root <- here::here("../data")
results_root <- here::here("../results")
plots_root <- here::here("../plots")

# Load data ====================================================================

metadata <- load_metadata(data_root)
metadata_dict <- read_csv(file.path(data_root, "metadata_dict.csv"))

wallthickness <- plot_data <- vector("list", 2) %>% set_names(c("ED", "ES"))

for (frame in c("ED", "ES")) {
  wallthickness[[frame]] <- load_wallthickness(data_root, frame)
  wallthickness[[frame]] <- inner_join(wallthickness[[frame]], metadata)
  
  # Load DDRtree
  ddrtree <- read_csv(
    file.path(results_root, glue("ddrtree_res_{frame}_decim0.99.csv"))
  ) %>%
    mutate(across(c(cluster, branch, branch_merged), ~ as.factor(.x)))
  if (frame == "ES") {
    ddrtree$Z1 <- -ddrtree$Z1
    ddrtree$Z2 <- -ddrtree$Z2
  }
  
  # Change the branch label for ES: 4 <-> 5
  # This is the result of the co-occurrence of clusters run separately
  if (frame == "ES") {
    ddrtree <- ddrtree %>%
      mutate(branch_merged = plyr::mapvalues(.$branch_merged,
                                             c("4", "5"),
                                             c("5", "4"))) %>%
      mutate(branch_merged = fct_relevel(branch_merged, "1", "2", "3", "4", "5"))
  }
  
  plot_data[[frame]] <- wallthickness[[frame]] %>%
    inner_join(ddrtree) %>%
    rename(Genotype = genotype) %>%
    mutate(
      Genotype = as.factor(Genotype),
      Branch = as.factor(branch_merged)
    )
}

# Plot trees ===================================================================

vars_of_interest <- c("is_plp", "Deceased")
name_var <- list("P/LP", "Deceased") %>% set_names(vars_of_interest)
name_plot <- list("Probability of P/LP", "Probability of mortality") %>% 
  set_names(vars_of_interest)

results_cv <- plot_tree_probs <- vector("list", 4) %>%
  set_names(c(paste0("ED.", vars_of_interest),
              paste0("ES.", vars_of_interest)))

for (frame in c("ED", "ES")) {
  model_data <- plot_data[[frame]]
  model_data$is_plp <- model.matrix(~Genotype, model_data)[, 2]
  
  for (i in seq_along(vars_of_interest)) {
    
    res_cv <- cross_validate_gam_tree(
      X = model_data %>% select(Z1, Z2),
      y = (model_data[[vars_of_interest[i]]] == 1) * 1,
      n_cv = 10
    )
    
    results_cv[[glue("{frame}.{vars_of_interest[i]}")]] <- res_cv
    
    print(median(res_cv$auc))
    print(IQR(res_cv$auc))
    
    model_gam <- gam(
      as.formula(glue("{vars_of_interest[i]} ~ -1 + te(Z1, Z2)")),
      data = model_data, family = binomial()
    )
    probs <- predict(model_gam, model_data %>% select(Z1, Z2),
                         type = "response") %>%
      logit_to_probs(.)
    
    lgPredObj <- ROCR::prediction(
      list(res_cv$pred_probs),
      factor(model_data[[vars_of_interest[i]]])
    )
    lgPerfObj <- performance(lgPredObj, "tpr", "fpr")
    # plotting ROC curve
    plot(lgPerfObj, main = "ROC Curve", col = 2, lwd = 2)
    abline(a = 0, b = 1, lwd = 2, lty = 3,col = "black")
    
    # Plot tree overlaid with P/LP probs
    plot_tree_probs[[glue("{frame}.{vars_of_interest[i]}")]] <- model_data %>%
      mutate(prob = probs) %>%
      ggplot(aes(x = Z1, y = Z2, color = prob)) +
      geom_point(size = 2, alpha = 0.75) +
      scale_color_paletteer_c(
        "grDevices::terrain.colors",
        name = NULL,
        limits = c(round(min(probs), 3), round(max(probs), 3)),
        breaks = seq(round(min(probs), 3), round(max(probs), 3), length.out = 4),
        labels = round(seq(min(probs), max(probs), length.out = 4), 3)
      ) +
      ggtitle(name_plot[[i]]) +
      tree_plot_params +
      theme(legend.key.size = unit(6, "pt"),
            legend.text = element_text(size = unit(6, "pt")),
            legend.title = element_text(size = unit(8, "pt")),
            legend.key.width = unit(0.5, "cm")) +
      coord_fixed()
  }
}

plot_tree_probs_combo <- plot_grid(plotlist = plot_tree_probs, align = "hv",
                                   nrow = 1, labels = c("A", NA, "B", NA))
ggsave(filename = file.path(plots_root, "plots_for_figures",
                            "plot_tree_probs.pdf"),
       dpi = 300,
       units = "mm",
       width = 210,
       height = 60,
       plot = plot_tree_probs_combo)

