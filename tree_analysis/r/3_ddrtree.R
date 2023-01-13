# Perform DDTree analysis
#
# Author: Paolo Inglese, Imperial College London, 2022


# Libraries ====================================================================

library(monocle)
library(glue)
library(tidyverse)
library(Seurat)
library(igraph)
library(DDRTree)
library(irlba)
library(magrittr)
library(cowplot)
library(here)
library(plyr)

# Setup ========================================================================

data_root <- here::here("../data")
results_root <- here::here("../results")
decim <- 0.99
frame <- "ED_ES"

if (frame == "ED") {
  sel_res <- "res.0.5"
} else if (frame == "ES") {
  sel_res <- "res.0.1"
} else if (frame == "ED_ES") {
  sel_res <- "res.0.1"
}

# Load data ====================================================================

# Load imputed metadata
metadata <- read_csv(
  file.path(data_root, "metadata_imputed.csv")
) %>%
  filter(genotype %in% c("NEG", "VUS", "PLP")) %>%
  filter(!duplicated(ID))

if (frame == "ED_ES") {
  wt_single_frame <- vector("list", 2) %>% set_names(c("ED", "ES"))
  for (fr in c("ED", "ES")) {
    wt_single_frame[[fr]] <- read.csv(
      file.path(
        data_root,
        paste0("wallthickness_", fr, "_decimated_", decim, ".csv")
      )
    ) %>%
      select(-X) %>%
      rename_with(~ gsub("X", paste0("wt-", fr, "-"), .x), starts_with("X")) %>%
      mutate_at("ID", list(~ str_split(., "_") %>%
                             map_chr(., 1))) %>%
      filter(!duplicated(ID))
  }
  wallthickness <- inner_join(wt_single_frame[["ED"]],
                              wt_single_frame[["ES"]])
  rm(wt_single_frame)
} else {
  wallthickness <- read.csv(
    file.path(
      data_root,
      paste0("wallthickness_", frame, "_decimated_", decim, ".csv")
    )
  ) %>%
    select(-X) %>%
    rename_with(~ gsub("X", "wt", .x), starts_with("X")) %>%
    mutate_at("ID", list(~ str_split(., "_") %>%
                           map_chr(., 1))) %>%
    filter(!duplicated(ID))
}

wallthickness <- inner_join(wallthickness, metadata)

partitions <- read.csv(
  file.path(here::here(
    "../results",
    paste0(
      "cluster_assignments_all_res_", frame, "_decim_",
      decim, ".csv"
    )
  ))
)
sel_partition <- factor(partitions[[sel_res]])

# Adjust for covariates ========================================================

wt_values_orig <- wallthickness %>%
  select(starts_with("wt")) %>%
  as.matrix()

rownames(wt_values_orig) <- wallthickness$ID
seu_obj <- CreateSeuratObject(counts = t(wt_values_orig), assay = "WT")
seu_obj <- SetAssayData(seu_obj, "data", t(wt_values_orig))
meta_data <- data.frame(
  sex = factor(wallthickness$sex),
  age = as.numeric(wallthickness$age_at_scan),
  race = factor(wallthickness$race)
) %>%
  set_rownames(wallthickness$ID)
seu_obj <- AddMetaData(seu_obj, meta_data)

seu_obj <- ScaleData(
  seu_obj,
  vars.to.regress = c("sex", "age", "race"), assay = "WT"
)
wt_values_adj <- t(seu_obj@assays$WT@scale.data)

rm(seu_obj)

# DDRTree analysis =============================================================

pca <- prcomp_irlba(wt_values_adj, n = 50)

obs_metadata <- data.frame(
  ID = metadata$ID,
  genotype = metadata$genotype
) %>%
  set_rownames(rownames(wt_values_adj))

feat_metadata <- data.frame(vert = seq(1, ncol(pca$x))) %>%
  set_rownames(colnames(pca$x))

cds <- monocle::newCellDataSet(t(pca$x),
  phenoData = AnnotatedDataFrame(obs_metadata),
  featureData = AnnotatedDataFrame(feat_metadata),
  expressionFamily = uninormal()
)

cds <- reduceDimension(cds,
  reduction_method = "DDRTree", norm_method = "none",
  pseudo_expr = 0, scaling = FALSE, verbose = TRUE,
  relative_expr = FALSE, ncenter = NULL, maxIter = 100,
  tol = 1e-6
)

cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "State")

Z <- cds@reducedDimS
Y <- cds@reducedDimK
branch <- as.numeric(cds$State)

if (frame == "ED") {
  merged_branch <- branch %>%
    plyr::mapvalues(., c(2, 3, 14, 15), c(1, 1, 1, 1)) %>%
    plyr::mapvalues(., c(4, 5, 6, 7, 8), c(2, 2, 2, 2, 2)) %>%
    plyr::mapvalues(., c(9, 10), c(3, 3)) %>%
    plyr::mapvalues(., c(11, 12, 13), c(4, 4, 4))
} else if (frame == "ES") {
  merged_branch <- branch %>%
    plyr::mapvalues(., c(1, 2, 3, 10, 11), c(1, 1, 1, 1, 1)) %>%
    plyr::mapvalues(., c(7, 8, 9), c(2, 2, 2)) %>%
    plyr::mapvalues(., c(5), c(3)) %>%
    plyr::mapvalues(., c(6), c(5))
} else if (frame == "ED_ES") {
  merged_branch <- branch %>%
    plyr::mapvalues(., c(3, 16, 17, 13, 14, 11, 10, 18, 5, 6, 2, 9),
                    c(1, 15, 15, 12, 12, 12, 12, 8, 8, 8, 1, 8)) %>%
    plyr::mapvalues(., c(4, 7, 19, 8, 15, 12), c(2, 3, 4, 5, 6, 7))
}

data.frame(Z1 = Z[1, ], Z2 = Z[2, ], branch = factor(merged_branch)) %>%
  ggplot(aes(x = Z1, y = Z2, color = branch)) +
  geom_point()

data.frame(
  ID = wallthickness$ID,
  Z1 = Z[1, ], Z2 = Z[2, ],
  Y1 = Y[1, ], Y2 = Y[2, ],
  cluster = sel_partition,
  branch = branch,
  branch_merged = merged_branch
) %>%
  write.csv(file = file.path(
    here::here("../results/"),
    paste0("ddrtree_res_", frame, "_decim", decim, ".csv")
  ))
