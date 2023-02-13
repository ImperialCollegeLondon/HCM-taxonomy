# Run unsupervised analysis on wall thickness from 3D meshes
#
# Author: Paolo Inglese, Imperial College, 2022


# Libraries ====================================================================

# Necessary to run umap-learn and calculate the trustworthiness
library(reticulate)
use_condaenv(
  condaenv = "hcm",
  conda = "C:/Users/pi514/AppData/Local/mambaforge/condabin/mamba.bat"
)
sklearn_py <- import("sklearn")

library(magrittr)
library(tidyverse)
library(stringr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(umap)
library(matrixStats)
library(ggbeeswarm)
library(clustree)
library(fpc)
library(purrr)
library(cowplot)
library(doParallel)
library(ggsignif)
library(here)

source(here::here("functions.R"))

# Setup ========================================================================

data_root <- here("../data")
results_root <- here("../results")
frame <- "ED"  # it can be either "ED" or "ES"

# Load data ====================================================================

metadata <- load_metadata(data_root)
wallthickness <- load_wallthickness(data_root, frame)
wallthickness <- inner_join(wallthickness, metadata)

# Adjust for covariates ========================================================

wt_values_orig <- wallthickness %>%
  select(starts_with("wt")) %>%
  as.matrix()

wt_values_adj <- adjust_wt_covariates(wt_values_orig, wallthickness)

# Calculate partitions =========================================================

test_res <- seq(0.1, 1, 0.1)

results_cluster <- cluster_wallthickness(
  wt_values_adj,
  test_res,
  dim_red = "PCA")

# Save clusters
results_cluster$clusters %>%
  mutate(ID = wallthickness$ID) %>%
  write.csv(., file = file.path(
    results_root, paste0(
      "cluster_assignments_all_res_", frame,
      "_decim_0.99.csv"
    )))

# Cluster Stability ============================================================

boot_res <- vector("list", length(test_res)) %>%
  set_names(paste0("res.", test_res))

# Stability test configuration
if (frame == 'ED') {
  stab_test_res <- c(0.4, 0.5)
} else if (frame == 'ES') {
  stab_test_res <- seq(0.1, 0.6, 0.1)
} else if (frame == "ED_ES") {
  stab_test_res <- seq(0.1, 0.6, 0.1)
}
stability_test_conf <- list(n_rep = 1e3, res = stab_test_res)

boot_res <- lapply(stability_test_conf$res, function(r) {
  clusterboot(wt_values_adj,
              B = stability_test_conf$n_rep,
              bootmethod = "subset", 
              clustermethod = partial(get_partitions_fpc, res = r),
              subtuning = round(nrow(wt_values_adj) * 0.8))
})

saveRDS(
  boot_res,
  file = file.path(
    results_root,
    paste0("subset_cluster_stability_", frame, "_decim_0.99.rds")
  )
)

# Print latex results

boot_res <- readRDS(file.path(
  results_root,
  paste0("subset_cluster_stability_", frame, "_decim_0.99.rds")
))

out_xtable <- lapply(boot_res, function(x) x$subsetmean)
n_cols <- max(lapply(out_xtable, length) %>% unlist())
mat_vals <- matrix(NA, length(out_xtable), n_cols,
  dimnames = list(names(boot_res), paste0("cluster.", seq_len(n_cols)))
)
for (i in seq_len(length(out_xtable))) {
  mat_vals[i, 1:length(out_xtable[[i]])] <- out_xtable[[i]]
}

mat_vals <- as.data.frame(mat_vals) %>%
  set_rownames(stability_test_conf$res)

# Print for LaTeX
print(xtable::xtable(mat_vals, digits = 4))

# Select optimal partitioning ==================================================

if (frame == "ES") {
  sel_res <- "res.0.1"
} else if (frame == "ED") {
  sel_res <- "res.0.5"
} else if (frame == "ED_ES") {
  sel_res <- "res.0.1"
}
sel_partition <- results_cluster$clusters[[sel_res]]

# UMAP =========================================================================

pca <- prcomp_irlba(wt_values_adj, n = 50)

# Calculate trustworthiness to choose the optimal metric
metrics <- c("euclidean", "cosine", "manhattan")

tw <- sapply(metrics, function(m) {
  umap_metric <- m
  umap_config <- umap.defaults
  umap_config$metric <- umap_metric
  
  mapped <- pca$x %>%
    set_rownames(paste0("obs.", seq_len(nrow(wallthickness)))) %>%
    umap(., config = umap_config, method = "umap-learn")
  
  sklearn_py$manifold$trustworthiness(
    wt_values_adj, mapped$layout,
    metric = umap_metric, n_neighbors = as.integer(umap_config$n_neighbors)
  )
}) %>%
  set_names(metrics)

# Print the trustworthiness for latex
data.frame(trustworthiness = tw) %>%
  rownames_to_column() %>%
  rename("metric" = "rowname") %>%
  xtable::xtable(., digits = 4)

# Calculate UMAP using the metric with the highest trustworthiness
umap_metric <- "cosine"
umap_config <- umap.defaults
umap_config$metric <- umap_metric
umap_config$random_state <- 123

mapped <- pca$x %>%
  set_rownames(make.names(paste0("obs.", seq_len(nrow(wallthickness))))) %>%
  umap(., config = umap_config, method = "umap-learn")

mapped_layout <- mapped$layout %>%
  set_colnames(c("UMAP1", "UMAP2")) %>%
  as_tibble() %>%
  mutate(
    ID = wallthickness$ID,
    Cluster = sel_partition
  )

plot_umap_clusters <- wallthickness %>%
  inner_join(mapped_layout) %>%
  rename(Genotype = genotype) %>%
  ggpubr::ggscatter(
    data = .,
    x = "UMAP1",
    y = "UMAP2",
    color = "Cluster",
    alpha = 0.5
  ) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed()

plot_umap_geno <- wallthickness %>%
  inner_join(mapped_layout) %>%
  rename(Genotype = genotype) %>%
  ggpubr::ggscatter(
    data = .,
    x = "UMAP1",
    y = "UMAP2",
    color = "Genotype",
    alpha = 0.5
  ) +
  scale_color_brewer(palette = "Set1") +
  coord_fixed()

plot_umap_combo <- plot_grid(
  plotlist = list(plot_umap_clusters, plot_umap_geno), nrow = 2,
  labels = c("A", "B")
)

rm(plot_umap_clusters, plot_umap_geno)

# ggsave(
#   filename = here(file.path(
#     "../plots", paste0("plot_umap_", frame, "_decim_", decim, ".pdf")
#   )),
#   dpi = 300, width = 12, height = 6, plot = p_cluster_combo
# )

# Enrichment of genotype =======================================================

pvals_list <- vector("list", ncol(results_cluster$clusters)) %>%
  set_names(colnames(results_cluster$clusters))
padj_list <- pvals_list

enrich_results <- expand.grid(
  sort(unique(sel_partition)),
  sort(unique(wallthickness$genotype))) %>%
  set_colnames(c("Cluster", "Genotype")) %>%
  mutate(p_val = NA,
         p_adj = NA)
for (i in seq_len(nrow(enrich_results))) {
  enrich_results$p_val[i] <- fisher.test(
    sel_partition == enrich_results$Cluster[i],
    wallthickness$genotype == enrich_results$Genotype[i],
    alternative = "greater"
  )$p.value
}
enrich_results$p_adj <- p.adjust(enrich_results$p_val, "BH")

# Print latex table
enrich_results %>%
  select(Cluster, Genotype, p_adj) %>%
  spread(., key = "Genotype", value = "p_adj") %>%
  xtable::xtable(., digits = 4)

plot_composition <- tibble(
  cluster = sel_partition,
  genotype = wallthickness$genotype
) %>%
  group_by(cluster, genotype) %>%
  tally() %>%
  rename(Genotype = genotype) %>%
  ggplot(aes(x = cluster, y = n, fill = Genotype)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_brewer(palette = "Set1") +
  xlab("Cluster") +
  ylab("Number of subjects") +
  theme_minimal_hgrid() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

plot_scatter <- plot_grid(
  plotlist = list(
    plot_umap_combo,
    plot_grid(results_cluster$plot_clustree,
              plot_composition, ncol = 1,
              labels = c("C", "D"))),
  nrow = 1, align = 'hv'
)

ggsave(
  file = paste0(
    "../plots/plot_umap_", frame, "_decim_0.99.pdf"
  ),
  plot = plot_scatter, width = 8, height = 8, dpi = 300
)

# WT vs genotype ===============================================================

vals <- wt_values_adj %>%
  as_tibble() %>%
  mutate(genotype = wallthickness$genotype) %>%
  rowwise() %>%
  mutate(
    wt_med = median(c_across(starts_with("wt"))),
    wt_sd = sd(c_across(starts_with("wt")))
  ) %>%
  ungroup() %>%
  mutate(Cluster = sel_partition)

# Median WT vs genotype

# Calculate p-values for pairwise comparisons
pvals_by_geno <- c()
pvals_by_geno <- c(
  pvals_by_geno,
  kruskal.test(
    vals$wt_med[vals$genotype %in% c("NEG", "PLP")],
    vals$genotype[vals$genotype %in% c("NEG", "PLP")]
  )$p.value
)
pvals_by_geno <- c(
  pvals_by_geno,
  kruskal.test(
    vals$wt_med[vals$genotype %in% c("VUS", "PLP")],
    vals$genotype[vals$genotype %in% c("VUS", "PLP")]
  )$p.value
)
pvals_by_geno <- c(
  pvals_by_geno,
  kruskal.test(
    vals$wt_med[vals$genotype %in% c("NEG", "VUS")],
    vals$genotype[vals$genotype %in% c("NEG", "VUS")]
  )$p.value
)
p_adj_by_geno <- p.adjust(pvals_by_geno, "BH")

# Manually set the coordinates for significance symbols
if (frame == "ED") {
  y_pos <- c(3, 2.3, 3.2)
  y_range <- c(-1.7, 3.5)
} else if (frame == "ES") {
  y_pos <- c(3.15, 2.95, 3.4)
  y_range <- c(-1.7, 3.5)
} else if (frame == "ED_ES") {
  y_pos <- c(3.15, 2.95, 3.4)
  y_range <- c(-1.7, 3.5)
}

p_by_geno_df <- data.frame(
  group1 = c("NEG", "VUS", "NEG"),
  group2 = c("P/LP", "P/LP", "VUS"),
  p_adj = p_adj_by_geno,
  y.position = y_pos,
  p.signif = gtools::stars.pval(p_adj_by_geno)
) %>%
  mutate(p.signif = replace(p.signif, p.signif == " ", "N.S."))

p_1 <- vals %>%
  select(wt_med, genotype) %>%
  rename(Genotype = genotype) %>%
  ungroup() %>%
  mutate(genotype = replace(Genotype, Genotype == "PLP", "P/LP")) %>%
  mutate(Cluster = sel_partition) %>%
  ggplot(aes(x = genotype, y = wt_med, color = Cluster)) +
  geom_violin(color = "black") +
  geom_boxplot(color = "black", width = 0.1) +
  geom_beeswarm() +
  stat_pvalue_manual(
    data = p_by_geno_df, label = "p.signif", map_signif_level = TRUE,
    size = 5
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_y_continuous(limits = y_range) +
  ylab("Median adjusted WT") +
  xlab("Genotype") +
  ggtitle("Median WT by genotype") +
  theme_minimal_hgrid()

# Median WT vs cluster

# Calculate p-values for pairwise comparisons
if (frame == "ED") {
  pvals_by_cluster <- c()
  pvals_by_cluster <- c(
    pvals_by_cluster,
    kruskal.test(
      vals$wt_med[vals$Cluster %in% c("0", "1")],
      vals$Cluster[vals$Cluster %in% c("0", "1")]
    )$p.value
  )
  pvals_by_cluster <- c(
    pvals_by_cluster,
    kruskal.test(
      vals$wt_med[vals$Cluster %in% c("1", "2")],
      vals$Cluster[vals$Cluster %in% c("1", "2")]
    )$p.value
  )
  pvals_by_cluster <- c(
    pvals_by_cluster,
    kruskal.test(
      vals$wt_med[vals$Cluster %in% c("0", "2")],
      vals$Cluster[vals$Cluster %in% c("0", "2")]
    )$p.value
  )
  p_adj_by_cluster <- p.adjust(pvals_by_cluster, "BH")
  y_pos <- c(3.05, 3.25, 3.45)
  p_df_by_cluster <- data.frame(
    group1 = c("0", "1", "0"),
    group2 = c("1", "2", "2"),
    p_adj = p_adj_by_cluster,
    y.position = y_pos,
    p.signif = gtools::stars.pval(p_adj_by_cluster)
  ) %>%
    mutate(p.signif = replace(p.signif, p.signif == " ", "N.S."))
} else if (frame %in% c("ES", "ED_ES")) {
  pvals_by_cluster <- c()
  pvals_by_cluster <- c(
    pvals_by_cluster,
    kruskal.test(
      vals$wt_med[vals$Cluster %in% c("0", "1")],
      vals$Cluster[vals$Cluster %in% c("0", "1")]
    )$p.value
  )
  p_adj_by_cluster <- p.adjust(pvals_by_cluster, "BH")
  y_pos <- c(2.95)
  p_df_by_cluster <- data.frame(
    group1 = "0",
    group2 = "1",
    p_adj = p_adj_by_cluster,
    y.position = y_pos,
    p.signif = gtools::stars.pval(p_adj_by_cluster)
  ) %>%
    mutate(p.signif = replace(p.signif, p.signif == " ", "N.S."))
}

p_2 <- vals %>%
  select(wt_med, genotype) %>%
  ungroup() %>%
  mutate(genotype = replace(genotype, genotype == "PLP", "P/LP")) %>%
  mutate(Cluster = sel_partition) %>%
  rename(Genotype = genotype) %>%
  ggplot(aes(x = Cluster, y = wt_med, color = Genotype)) +
  geom_violin(color = "black", trim = TRUE) +
  geom_boxplot(color = "black", width = 0.1) +
  geom_beeswarm() +
  stat_pvalue_manual(
    data = p_df_by_cluster, label = "p.signif",
    map_signif_level = TRUE, size = 5
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = y_range) +
  ylab("Median adjusted WT") +
  xlab("Cluster") +
  ggtitle("Median WT by cluster") +
  theme_minimal_hgrid()

p_wt <- cowplot::plot_grid(
  plotlist = list(p_1, p_2), labels = c("A", "B"),
  nrow = 1
)
ggsave(
  filename = paste0("../plots/median_wt_", frame, "_decim_0.99.pdf"),
  plot = p_wt, dpi = 300, width = 12, height = 6
)
