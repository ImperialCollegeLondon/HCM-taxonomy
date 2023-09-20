# Analysis of correlations between DDRtree and clinical features
#
# Author: Paolo Inglese, Imperial College 2022


# Libraries ====================================================================

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
frame <- "ES"  #  This can be either "ED" or "ES"

if (frame == "ED") {
  sel_res <- "res.0.5"
} else if (frame == "ES") {
  sel_res <- "res.0.1"
} else if (frame == "ED_ES") {
  sel_res <- "res.0.1"
}

# Load data ====================================================================

metadata <- load_metadata(data_root)
metadata_dict <- read_csv(file.path(data_root, "metadata_dict.csv"))
wallthickness <- load_wallthickness(data_root, frame)
wallthickness <- inner_join(wallthickness, metadata)

# Load clusters
sel_partition <- read_csv(
  file.path(
    results_root, glue("cluster_assignments_all_res_{frame}_decim_0.99.csv")
  )
) %>%
  pull(matches(sel_res)) %>%
  as.factor()

# Load DDRtree
ddrtree <- read_csv(
  file.path(results_root, glue("ddrtree_res_{frame}_decim0.99.csv"))
) %>%
  mutate(across(c(cluster, branch, branch_merged), ~ as.factor(.x)))

# Change the branch label for ES: 4 <-> 5
# This is the result of the co-occurrence of clusters run separately
if (frame == "ES") {
  ddrtree <- ddrtree %>%
    mutate(branch_merged = plyr::mapvalues(.$branch_merged,
                                           c("4", "5"),
                                           c("5", "4"))) %>%
    mutate(branch_merged = fct_relevel(branch_merged, "1", "2", "3", "4", "5"))
}

# Test all variables ===========================================================

cont_vars <- metadata_dict %>%
  filter(value %in% "continuous") %>%
  filter(variable != "age_at_scan") %>%
  pull(variable)

cat_vars <- metadata_dict %>%
  filter(value %in% c("logical", "factor", "discrete")) %>%
  pull(variable)

p_vals_cont <- lapply(cont_vars, function(x) {
  print(x)
  df_ <- data.frame(
    branch = factor(ddrtree$branch_merged),
    feature = wallthickness[[x]]
  )
  kruskal.test(feature ~ branch, data = df_)$p.value
}) %>%
  Reduce(c, .) %>%
  set_names(cont_vars)

p_vals_cat <- lapply(cat_vars, function(x) {
  print(x)
  chisq.test(
    x = factor(ddrtree$branch_merged),
    y = factor(wallthickness[[x]]),
    simulate.p.value = TRUE
  )$p.value
}) %>%
  Reduce(c, .) %>%
  set_names(cat_vars)

p_adj_branch <- tibble(
  p_val = p.adjust(c(p_vals_cont, p_vals_cat), "BH"),
  type = c(
    rep("continuous", length(p_vals_cont)),
    rep("discrete", length(p_vals_cat))
  ),
  variable = c(cont_vars, cat_vars)
) %>%
  filter(p_val < 0.05) %>%
  inner_join(metadata_dict)

# Pairwise tests continuous vars ===============================================

sel_cont_vars <- p_adj_branch %>%
  filter(type == "continuous")

# Perform Dunn test for each selected var
p_pairwise <- lapply(sel_cont_vars$variable, function(x) {
  dunn.test::dunn.test(
    wallthickness[[x]], factor(ddrtree$branch_merged),
    method = "bh"
  )
}) %>% set_names(sel_cont_vars$variable)

p_pairwise <- p_pairwise[
  lapply(p_pairwise, function(x) { sum(x$P.adjusted < 0.05) > 0 }) %>% unlist()
]

sel_cont_vars <- sel_cont_vars %>%
  filter(variable %in% names(p_pairwise))

# p_pairwise <- vector("list", nrow(sel_cont_vars)) %>%
#   set_names(sel_cont_vars)

plot_pairwise <- vector("list", nrow(sel_cont_vars)) %>%
  set_names(sel_cont_vars)

# Positions of significance bars
if (frame == "ED") {
  step_0 <- c(1, 0.02, 2, 5, 0.07, 0.5, 1)
  step_size <- c(0.5, 0.02, 2, 5, 0.05, 0.4, 0.5)
} else if (frame == "ES") {
  step_0 <- c(1, 1, 0.02, 1, 0.5, 1, 0.5, 1)
  step_size <- c(0.5, 1, 0.02, 2, 1, 2, 0.2, 1)
} else if (frame == "ED_ES") {
  step_0 <- c(1, 1, 0.02, 1, 0.5, 1, 0.5, 1, 1, 1)
  step_size <- c(0.5, 1, 0.02, 2, 1, 2, 0.2, 1, 1, 1)
}

for (n in seq_len(nrow(sel_cont_vars))) {
  varname <- sel_cont_vars$variable[n]
  
  dat_tmp <- tibble(
    values = wallthickness[[varname]],
    Branch = ddrtree$branch_merged
  ) %>%
    group_by(Branch) %>%
    summarise(
      med = median(values),
      lo = quantile(values, 0.25),
      hi = quantile(values, 0.75)
    )
  
  stat_tmp <- tibble(
    group1 = str_split(
      p_pairwise[[varname]]$comparisons, " -",
      simplify = TRUE
    )[, 1],
    group2 = str_split(
      p_pairwise[[varname]]$comparisons, "- ",
      simplify = TRUE
    )[, 2],
    label = "p.adj",
    p.adj = p_pairwise[[varname]]$P.adjusted,
    y.position = NA
  ) %>%
    filter(p.adj < 0.05)
  
  ordered_groups <- order(
    abs(as.numeric(stat_tmp$group1) - as.numeric(stat_tmp$group2))
  )
  
  stat_tmp <- stat_tmp[ordered_groups, ]
  
  stat_tmp <- stat_tmp %>%
    rowwise() %>%
    mutate(
      max_y_between =
        max(dat_tmp$hi[dat_tmp$Branch %in%
                         c(as.numeric(group1):as.numeric(group2))])
    ) %>%
    arrange(max_y_between) %>%
    group_by(max_y_between) %>%
    mutate(y.position = sapply(
      1:length(max_y_between),
      function(k) {
        (k - 1) * step_0[n] + step_size[n] +
          ifelse(k == 0, max_y_between[k], max(max_y_between[1:k]))
      }
    )) %>%
    ungroup()
  
  y_pos <- c()
  for (k in seq_len(nrow(stat_tmp))) {
    y_ <- stat_tmp$max_y_between[k] + step_0[n]
    if (k > 1) {
      while (any(y_ <= y_pos[1:(k - 1)] + step_size[n])) {
        y_ <- y_ + step_size[n]
      }
    }
    y_pos <- c(y_pos, y_)
  }
  
  # Add y positions
  stat_tmp <- stat_tmp %>%
    mutate(y.position = y_pos,
           x = apply(stat_tmp, 1, function(x) mean(as.numeric(c(x[1], x[2])))),
           signif = gtools::stars.pval(stat_tmp$p.adj))
  
  plot_pairwise[[n]] <- dat_tmp %>%
    ggplot(aes(x = Branch,
               y = med,
               ymin = lo,
               ymax = hi,
               color = Branch)) +
    geom_point(size = 2) +
    geom_errorbar(width = 0.1, linewidth = 1) +
    scale_x_discrete(limits = rev) +
    stat_pvalue_manual(
      data = stat_tmp,
      label = "signif",
      inherit.aes = FALSE, 
      tip.length = 0,
      angle = 270,
      size = unit(4, "pt"),
      coord.flip = TRUE) +
    coord_flip() +
    scale_colour_branded() +
    ggtitle(sel_cont_vars$name_long[n]) +
    ylab(sel_cont_vars$name_short[n]) +
    theme_classic() +
    theme(legend.position = "none",
          text = element_text(size = unit(6, "pt")),
          axis.title = element_text(size = unit(6, "pt")),
          axis.text = element_text(size = unit(6, "pt")),
          plot.title = element_text(size = unit(8, "pt"))) +
    guides(color = guide_legend(title = sel_cont_vars$name_short[n]))
}

plot_cont_vars <- plot_grid(plotlist = plot_pairwise, align = "hv")

# Pairwise tests categorical vars ==============================================

sel_cat_vars <- p_adj_branch %>%
  filter(type == "discrete")

# Calculate p-values for selected variables

calculate_categorical_pvalues <- function(branch, cat_var) {
  ub <- sort(unique(branch))
  ux <- sort(unique(cat_var))
  
  pmat <- matrix(NA, length(ub), length(ux),
                 dimnames = list(paste0("branch_", ub), paste0("pval_", ux))
  )
  
  for (i in seq_along(ub)) {
    for (j in seq_along(ux)) {
      pmat[i, j] <- fisher.test(branch == ub[i],
                                cat_var == ux[j],
                                alternative = "greater"
      )$p.value
    }
  }
  
  padj <- matrix(p.adjust(pmat, "BH"), nrow(pmat), ncol(pmat),
                 dimnames = list(NULL, paste0("padj_", ux))
  )
  
  pmat <- cbind(pmat, padj)
  return(pmat)
}

p_cat_branch <- vector("list", nrow(sel_cat_vars)) %>%
  set_names(sel_cat_vars$variable)

for (n in seq_len(nrow(sel_cat_vars))) {
  p_cat_branch[[sel_cat_vars$variable[n]]] <- calculate_categorical_pvalues(
    ddrtree$branch_merged,
    wallthickness[[sel_cat_vars$variable[n]]]
  )
}

p_cat_branch <- p_cat_branch[
  lapply(p_cat_branch, function(x) {
    sum((x %>%
    as_tibble() %>%
      dplyr::select(starts_with("padj"))) < 0.05) > 0
  }) %>%
    unlist()
]

sel_cat_vars <- sel_cat_vars %>%
  filter(variable %in% names(p_cat_branch))

# Barplots of the selected categorical variables

plot_bar_categorical <- lapply(seq_len(nrow(sel_cat_vars)), function(idx) {
  
  vx <- sel_cat_vars$variable[idx]
  
  p_table <- p_cat_branch[[vx]] %>%
    as_tibble() %>%
    select(starts_with("padj")) %>%
    mutate_all(gtools::stars.pval) %>%
    mutate_all(function(x) ifelse(x == " ", "N.S.", x)) %>%
    rownames_to_column("Branch") %>%
    pivot_longer(cols = colnames(.)[2:ncol(.)]) %>%
    filter(value != "N.S.") %>%
    mutate(Value = sapply(.$name, function(x) strsplit(x, "_")[[1]][2])) %>%
    select(-name)
  
  # Calculate positions of labels
  anno_data <- tibble(
    Branch = factor(ddrtree$branch_merged),
    Value = factor(wallthickness[[vx]])
  ) %>%
    group_by(Branch, Value) %>%
    count() %>%
    filter(Branch %in% p_table$Branch) %>%
    group_by(Branch) %>%
    arrange(desc(Value)) %>%
    mutate(left = c(0, cumsum(n[-length(n)]))) %>%
    mutate(right = cumsum(n)) %>%
    ungroup() %>%
    mutate(centre = rowMeans(cbind(left, right))) %>%
    mutate(Value = factor(Value)) %>%
    inner_join(p_table, by = c("Branch", "Value"))
  
  p_ <- tibble(
    Branch = factor(ddrtree$branch_merged),
    Value = factor(wallthickness[[vx]])) %>%
    group_by(Branch, Value) %>%
    count() %>%
    ggplot(aes(x = Branch, y = n, fill = Value)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(
      data = anno_data,
      mapping = aes(x = Branch, y = centre, label = value),
      inherit.aes = FALSE, color = "black",
      size = unit(4, "pt"), hjust = "center", vjust = "center") +
    scale_x_discrete(limits = rev) +
    ylab("Num. subjects") +
    ggtitle(n) +
    coord_flip() +
    theme_classic() +
    ggtitle(sel_cat_vars$name_long[idx])
  
  if (max(wallthickness[[vx]]) == 2) {
    p_ <- p_ + scale_fill_brewer(
      palette = "Dark2", type = "Seq", labels = c("Yes", "No"))
  } else {
    p_ <- p_ + scale_fill_brewer(palette = "Dark2", type = "Seq")
  }
  p_ + 
    theme(
      text = element_text(size = unit(6, "pt")),
      axis.title = element_text(size = unit(6, "pt")),
      axis.text = element_text(size = unit(6, "pt")),
      plot.title = element_text(size = unit(8, "pt")),
      legend.title = element_text(size = unit(6, "pt")),
      legend.key.size = unit(6, "pt"),
      legend.text = element_text(size = unit(6, "pt")),
      legend.box.margin = margin(0, 0, 0, 0))

})

plot_cat_vars <- plot_grid(plotlist = plot_bar_categorical, align = "hv")

# Plot tree ====================================================================

plot_data <- wallthickness %>%
  inner_join(ddrtree) %>%
  rename(Genotype = genotype) %>%
  mutate(
    Cluster = sel_partition,
    Genotype = as.factor(Genotype),
    Branch = as.factor(branch_merged)
  )

p_tree_branch <- plot_data %>%
  ggplot(aes(x = Z1, y = Z2, fill = Branch, color = Branch)) +
  geom_point(
    x = plot_data$Y1, y = plot_data$Y2,
    shape = 20,
    inherit.aes = F,
    size = 1,
    color = "lightgrey"
  ) +
  scale_x_continuous(limits = c(-25, 20)) +
  scale_y_continuous(limits = c(-12, 10)) +
  geom_point(size = 2, shape = 21, alpha = 0.7) +
  scale_colour_branded() +
  scale_fill_branded() +
  ggtitle(paste0("DDRTree branches - ", frame)) +
  theme_classic() +
  theme(
    legend.key.size = unit(0.1, "cm"),
    legend.position = "top"
  ) + coord_fixed()


plot_branch_combo <- plot_grid(
  
  plotlist = list(
   p_tree_branch,
   plot_grid(plotlist = list(plot_cont_vars, plot_cat_vars), nrow = 1)
  ),
  nrow = 2
)

ggsave(
  filename = file.path(
    plots_root,
    paste0("ddrtree_test_branches_", frame, "_decim_0.99.pdf")
  ),
  dpi = 300, units = "mm", width = 230, height = 220, plot = plot_branch_combo
)

# Plot the tree coloured by the selected variables =============================

plot_tree <- function(plot_data, var_data, is_discrete) {
  plot_out <- plot_data %>%
    ggplot(
      aes_string(
        x = "Z1_fixed",
        y = "Z2_fixed",
        color = var_data$variable)
    ) +
    # geom_point(
    #   x = plot_data$Y1_fixed,
    #   y = plot_data$Y2_fixed,
    #   shape = 20,
    #   size = 1,
    #   color = "lightgrey",
    #   inherit.aes = FALSE
    # ) +
    geom_point(size = 1, alpha = 0.5) +
    xlab("Z1") +
    ylab("Z2")
  
  if (is_discrete) {
    unique_labels <- sort(unique(plot_data[[var_data$variable]]))
    if (length(unique_labels) == 2) {
      if (all(unique_labels == c(1, 2))) {
        legend_labels <- c("Yes", "No")
        plot_out <- plot_out +
          scale_color_brewer(
            palette = "Dark2",
            type = "Seq",
            labels = legend_labels,
            name = var_data$name_short
          )
      } else {
        print(var_data$variable)
        legend_labels <- c("No", "Yes")  # this is for `is_plp`
        plot_out <- plot_out +
          scale_color_manual(
            values = brewer.pal(3, "Dark2")[2:1],
            labels = legend_labels,
            name = var_data$name_short
          ) +
          guides(color = guide_legend(reverse = TRUE))
      }
    } else {
      plot_out <- plot_out +
        scale_color_brewer(
          palette = "Dark2",
          type = "Seq",
          name = var_data$name_short
        )
    }
  } else {
    plot_out <- plot_out +
      scale_color_paletteer_c(
        "ggthemes::Red-Blue Diverging",
        name = var_data$name_short,
        trans = "reverse"
      ) +
      guides(color = guide_legend(reverse = TRUE))
  }
  
  plot_out <- plot_out +
    ggtitle(var_data$name_long) +
    tree_plot_params +
    theme(plot.title = element_text(size = unit(6, "pt"), face = "bold"),
          legend.position = "right")
  
  return(plot_out)
}

# sel_vars <- tibble(
#   var_name = c(sel_cont_vars, sel_cat_vars),
#   type = c(
#     rep("continuous", length(sel_cont_vars)),
#     rep("discrete", length(sel_cat_vars))
#   )
# ) %>%
#   mutate(
#     name = c(
#       "Weight", "BSA", "EDV", "Left Ventricle Mass",
#       "LVOTO peak vel.",
#       "Max LV Wall Thickness", "SV", "ACE/ARB", "ASA/Clopi",
#       "Beta Blocker",
#       "Hypertension", "Coinc. Infarction", "LVOTO",
#       "Most affected level"
#     ),
#     short_name = c(
#       "weight", "BSA", "EDV", "LVM", "vel",
#       "WT", "SV", "ACE/ARB", "ASA/Clopi", "Beta Blocker",
#       "HT", "Coinc. Infarction", "LVOTO", "Most. Aff. Level"
#     )
#   )

sel_vars <- rbind(sel_cont_vars, sel_cat_vars)

plot_data <- wallthickness %>%
  # select(ID,
  #        genotype,
  #        matches(sel_vars$variable)) %>%
  inner_join(ddrtree) %>%
  rename(Genotype = genotype) %>%
  mutate(
    Cluster = sel_partition,
    Genotype = as.factor(Genotype),
    Branch = as.factor(branch_merged)
  ) %>%
  mutate(
    across(
      which(colnames(.) %in%
              sel_vars$variable[sel_vars$type == "discrete"]),
      ~ as.factor(.x)
    )
  )

if (frame == "ES") {
  plot_data$Z1_fixed <- -plot_data$Z1
  plot_data$Z2_fixed <- -plot_data$Z2
  plot_data$Y1_fixed <- -plot_data$Y1
  plot_data$Y2_fixed <- -plot_data$Y2
} else {
  plot_data$Z1_fixed <- plot_data$Z1
  plot_data$Z2_fixed <- plot_data$Z2
  plot_data$Y1_fixed <- plot_data$Y1
  plot_data$Y2_fixed <- plot_data$Y2
}

# Plot tree continuous variables

plot_tree_features_cont <- lapply(1:nrow(sel_cont_vars), function(i) {
  plot_tree(plot_data, sel_cont_vars[i, ], FALSE)
})

plot_tree_features_cat <- lapply(seq_len(nrow(sel_cat_vars)), function(i) {
  plot_tree(plot_data, sel_cat_vars[i, ], TRUE)
})

plot_tree_features_combo <- plot_grid(
  plotlist = c(plot_tree_features_cont,
               plot_tree_features_cat),
          ncol = 2, align = "hv")

ggsave(
  filename = file.path(plots_root,
    paste0("plot_trees_", frame, "_decim_0.99_branch_vars.pdf")
  ),
  dpi = 300,
  units = "mm",
  width = 105,
  height = ifelse(frame == "ED", 240, 240),
  plot = ggplot_gtable(ggplot_build(plot_tree_features_combo))
)
