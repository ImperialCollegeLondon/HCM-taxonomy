# Functions used in the scripts
#
# Author: Paolo Inglese, Imperial College, 2022

library(magrittr)
library(dplyr)


# Estimate clusters using SNN and multilevel Louvain
find_clusters <- function(X, res) {
  require("Seurat")
  require("magrittr")
  
  if (is.null(rownames(X))) {
    X <- X %>%
      set_rownames(make.names(paste0("obs.", seq_len(nrow(X)))))
  }
  
  snn_obj <- X %>%
    FindNeighbors(., annoy.metric = "cosine", verbose = FALSE)
  
  partitions <- FindClusters(
    snn_obj$snn,
    algorithm = 2,
    resolution = res,
    verbose = FALSE
  )
  
  return(partitions)
}


# Helper function used to estimate cluster stability
get_partitions_fpc <- function(X, res) {
  parts <- find_clusters(X, res)
  parts <- parts[[paste0("res.", res)]]
  parts <- as.numeric(parts)
  up <- unique(parts)
  
  clusterlist <- vector("list", length(up)) %>%
    set_names(up)
  for (k in seq_along(up)) {
    clusterlist[[k]] <- (parts == up[k])
  }
  
  list(
    result = NULL,
    nc = length(up),
    clusterlist = clusterlist,
    partition = parts
  )
}

# Load imputed metadata
load_metadata <- function(data_root) {
  metadata <- read_csv(
    file.path(data_root, "metadata_imputed.csv")
  ) %>%
    filter(genotype %in% c("NEG", "VUS", "PLP")) %>%
    filter(!duplicated(ID)) %>%
    mutate(genotype = replace(genotype, genotype == "PLP", "P/LP")) %>%
    mutate(genotype = factor(genotype),
           is_plp = genotype == "P/LP")
  
  return(metadata)
}


# Load wall thickness data
load_wallthickness <- function(data_root, frame, decim = 0.99, dataset = "rbh",
                               load_curvature = FALSE) {
  
  require(glue)
  
  if (!frame %in% c("ED", "ES", "ED_ES")) {
    stop("frame must be one of 'ED', 'ES', 'ED_ES'")
  }
  
  if (frame == "ED_ES") {
    wt_single_frame <- vector("list", 2) %>% set_names(c("ED", "ES"))
    for (fr in c("ED", "ES")) {
      
      if (dataset == "rbh") {
        fname <- glue("wallthickness_{fr}_decimated_{decim}.csv")
      } else if (dataset == "singapore") {
        fname <- glue("wallthickness_{fr}_decimated_{decim}_singapore.csv")
      }
      
      wt_single_frame[[fr]] <- read.csv(
        file.path(
          data_root,
          fname
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
    
    if (dataset == "rbh") {
      wt_fname <- glue("wallthickness_{frame}_decimated_{decim}.csv")
      curv_fname <- glue("curvature_{frame}_decimated_{decim}.csv")
    } else if (dataset == "singapore") {
      wt_fname <- glue("wallthickness_{frame}_decimated_{decim}_singapore.csv")
      curv_fname <- glue("curvature_{frame}_decimated_{decim}_singapore.csv")
    }
    
    wallthickness <- read.csv(
      file.path(
        data_root,
        wt_fname
      )
    ) %>%
      select(-X) %>%
      rename_with(~ gsub("X", "wt", .x), starts_with("X"))
    
    if (load_curvature) {
      curvature <- read.csv(
        file.path(data_root, curv_fname)
      ) %>%
        select(-X) %>%
        rename_with(~ gsub("X", "curv", .x), starts_with("X"))
      
      wallthickness <- inner_join(wallthickness, curvature)
    }
    
    wallthickness <- wallthickness %>%
      mutate_at("ID", list(~ str_split(., "_") %>%
                             map_chr(., 1))) %>%
      filter(!duplicated(ID))
  }
  
  return(wallthickness)
  
}


# Adjust wall thickness for covariates
adjust_wt_covariates <- function(wt_values, wallthickness) {
  
  require(Seurat)
  
  # wallthickness is the data frame containing the covariates
  # wt_values is the matrix with the only wall thickness values
  
  if (is.null(rownames(wt_values))) {
    rownames(wt_values) <- wallthickness$ID
  }
  
  seu_obj <- CreateSeuratObject(counts = t(wt_values), assay = "WT")
  seu_obj <- SetAssayData(seu_obj, "data", t(wt_values))
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
  
  return(wt_values_adj)
  
}


# Perform clustering
cluster_wallthickness <- function(wt_values, test_resolution, dim_red) {
  
  require("irlba")
  
  if (dim_red == "PCA") {
    pca <- irlba::prcomp_irlba(wt_values, n = 50)
    input_X <- pca$x
  } else if (dim_red == "none") {
    input_X <- wt_values
  } else {
    stop("'dim_red' must be either 'PCA' or 'none'")
  }
  
  partitions <- find_clusters(input_X, test_resolution)
  
  # Clustree
  p_clustree <- clustree(
    partitions,
    prefix = "res.",
    prop_filt = 0.05,
    node_size = 5
  ) +
    guides(
      edge_colour = "none",
      edge_alpha = "none",
      point_size = "none",
      node_size = "none"
    ) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    theme(
      legend.key.size = unit(0.1, "cm"),
      legend.position = "top",
      aspect.ratio = 1
    )

  # orig_partitions <- vector("list", ncol(partitions)) %>%
  #   set_names(colnames(partitions))
  # for (i in 1:ncol(partitions)) {
  #   orig_partitions[[i]] <- data.frame(
  #     obs = rownames(partitions),
  #     partition_orig = partitions[, i]
  #   )
  # }
  
  return(list(plot_clustree = p_clustree, clusters = partitions))
  
}


# Set the colour map
branded_colors <- list(
  "blue" = "#00798c", 
         "red" = "#d1495b",
         "yellow" = "#edae49",
         "green" = "#66a182",
         "navy" = "#2e4057",
         "grey" = "#8d96a3",
         "darkgrey" = "#f37735")
         
branded_pal <- function(
    primary = "blue",
    other = "grey",
    direction = 1
) {
  stopifnot(primary %in% names(branded_colors))
  
  function(n) {
    if (n > 7) warning("Branded Color Palette only has 7 colors.")
    
    if (n == 2) {
      other <- if (!other %in% names(branded_colors)) {
        other
      } else {
        branded_colors[other]
      }
      color_list <- c(other, branded_colors[primary])
    } else {
      color_list <- branded_colors[1:n]
    }
    
    color_list <- unname(unlist(color_list))
    if (direction >= 0) color_list else rev(color_list)
  }
}

scale_colour_branded <- function(
    primary = "blue",
    other = "grey",
    direction = 1,
    ...
) {
  ggplot2::discrete_scale(
    "colour", "branded",
    branded_pal(primary, other, direction),
    ...
  )
}

scale_fill_branded <- function(
    primary = "blue",
    other = "grey",
    direction = 1,
    ...
) {
  ggplot2::discrete_scale(
    "fill", "branded",
    branded_pal(primary, other, direction),
    ...
  )
}


# Logits to probs
logit_to_probs <- function(logit) {
  odds <- exp(logit)
  probs <- odds / (1 + odds)
  return(probs)
}


# Perform cross validation GAM
cross_validate_gam_tree <- function(form, X, n_cv, type, method) {
  require(mgcv)
  require(caret)
  require(pROC)
  require(ROCR)
  
  set.seed(123)
  flds <- createFolds(X$Y, k = n_cv, list = TRUE, returnTrain = FALSE)
  
  all_preds <- numeric(nrow(X))
  
  pb <- progress::progress_bar$new(total = n_cv, width = 60,
                                   format = "CV [:bar]", clear = F)
  
  if (type == "classification") {
    family <- binomial()
  } else if (type == "regression") {
    family <- gaussian()
  }
  
  for (cv in seq_len(n_cv)) {
    pb$tick()
    dat_train <- X[-flds[[cv]], ]
    mdl <- gam(form, data = dat_train, family = family,
               method = method, select = TRUE)
    preds <- predict(mdl, newdata = X[flds[[cv]], ], type = "response")
    all_preds[flds[[cv]]] <- preds
  }
  pb$terminate()
  
  if (type == "classification") {
    perf_all <- auc(X$Y, all_preds)
  } else if (type == "regression") {
    perf_all <- cor(X$Y, all_preds) ** 2
  }
  
  return(list(perf_all = perf_all, preds = all_preds))
}