# Analysis of correlations between DDRtree and clinical features
#
# Author: Paolo Inglese, Imperial College 2022


# Libraries ====================================================================

library(igraph)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(glue)
library(here)
library(cowplot)
library(monocle)


source(here::here("functions.R"))

results_root = "../results/"

frame = "ES"

# Load the monocle object - this should already have the states calculated by
# the function orderCells
cds = readRDS(file.path(results_root, paste0("cds_", frame, ".rds")))

tree_backbone = cds@auxOrderingData$DDRTree$pr_graph_cell_proj_dist
dp_mst = cds@auxOrderingData$DDRTree$pr_graph_cell_proj_tree
reduced_dims = reducedDimS(cds)

# Fix the orientation of the tree for ES
if (frame == "ES") {
  cds@reducedDimS[2, ] = -reducedDimS(cds)[2, ]
  cds@reducedDimK[2, ] = -reducedDimK(cds)[2, ]
}

# Find the root state ==========================================================

# The root state corresponds to the centre of the tree

find_tree_centre_node = function(tree_backbone) {
  centre_pos = apply(tree_backbone, 1, median)
  dist_from_centre = dist(rbind(centre_pos, t(tree_backbone))) %>% as.matrix()
  return(names(which.min(dist_from_centre[-1, 1])))
}

tree_centre_node = find_tree_centre_node(tree_backbone)
root_state = cds$State[colnames(cds) == tree_centre_node]
message("Root state: ", root_state)

unique_states = sort(unique(cds$State))

# Find leaf branches to merge ==================================================

# Find branch and leaf nodes

adj_mat = get.adjacency(dp_mst)
is_leaf = which(apply(adj_mat, 1, function(x) {
  sum(x == 1) == 1
}))
is_branch = which(apply(adj_mat, 1, function(x) {
  sum(x == 1) == 3
}))

leaf_states = sort(cds$State[is_leaf])

states_connected_to_branch_nodes = vector("list", length(is_branch))
for (i in seq_along(is_branch)) {
  states_connected_to_branch_nodes[[i]] = sort(
    unique(
      cds$State[adj_mat[is_branch[i], ] != 0] %>%
        as.character() %>%
        as.numeric()))
}
names(states_connected_to_branch_nodes) = names(is_branch)

# Find the states connected to each branch

conn_branch_states = vector("list", length(is_branch)) %>%
  set_names(names(is_branch))
for (i in seq_along(is_branch)) {
  nn_idx = which(adj_mat[is_branch[i], ] != 0)
  conn_branch_states[[i]] = cds$State[nn_idx]
}

# Determine the short leaf states - these correspond to small variability.

cds$State = as.numeric(as.character(cds$State))
unique_states = sort(unique(cds$State))

get_geodesic_length_branches = function(dp_mst, cds) {
  
  us = sort(unique(cds$State))
  tb = cds@auxOrderingData$DDRTree$pr_graph_cell_proj_dist  # tree backbone
  lb = vector("list", length(us)) %>%
    set_names(us)
  for (i in seq_along(us)) {
    
    # Extract sub graph corresponding to each state
    subtree = igraph::subgraph(dp_mst, vids = colnames(cds)[cds$State == us[i]])
    edges_df = igraph::as_data_frame(subtree, "edges")
    
    # Calculate the length as the sum of the euclidean distance between the connected
    # nodes. Find them by looking at the edges
    ll = 0
    for (j in 1:nrow(edges_df)) {
      idx_v1 = which(colnames(cds) == edges_df[j, 1])
      idx_v2 = which(colnames(cds) == edges_df[j, 2])
      ll = ll + sum((tb[, idx_v1] - tb[, idx_v2]) ** 2)
    }
    
    lb[[as.character(us[i])]] = ll
  }
  
  return(lb)
}

length_branches = get_geodesic_length_branches(dp_mst, cds)

length_leafs = lapply(leaf_states, function(x) length_branches[[as.character(x)]]) %>%
  set_names(leaf_states)
graph_length = max(apply(tree_backbone, 1, function(x) max(x) - min(x)))
rel_state_lengths = lapply(length_leafs, function(x) x / graph_length)
small_states = names(rel_state_lengths)[unlist(rel_state_lengths) < 0.05]

message("Short leaf states: ", paste0(small_states, collapse = ", "))

# Merge them with the closest state.

conn_states_to_small = vector("list", length(small_states)) %>%
  set_names(small_states)
for (i in seq_along(small_states)) {
  tmp_list = list()
  for (j in seq_along(conn_branch_states)) {
    if (small_states[i] %in% conn_branch_states[[j]]) {
      tmp_list = append(tmp_list, conn_branch_states[[j]])
    }
  }
  conn_states_to_small[[i]] = unique(tmp_list)
  rm(tmp_list)
}

# 1. Merge short ===============================================================

p1 = plot_cell_trajectory(cds, color_by = factor(cds$State)) + coord_fixed()

old_state = as.numeric(small_states)

new_state = array(NA, length(old_state))
for (i in seq_along(old_state)) {
  
  # Calculate distance between state points and others connected
  idx_state = which(cds$State == old_state[i])
  nn_states = conn_states_to_small[[as.character(old_state[i])]]
  nn_states = as.numeric(as.character(nn_states[!nn_states %in% old_state[i]]))
  
  # Prefer merging with the closest connected leaf is present, otherwise to
  # the longest connected branch
  if (sum(nn_states %in% leaf_states) != 0) {
    nn_states = nn_states[nn_states %in% leaf_states]
    idx_conn = which(cds$State %in% nn_states)

    nns = data.frame(distance=numeric(), name=character())
    for (j in seq_along(idx_state)) {
      d = as.matrix(dist(rbind(reduced_dims[, idx_state[j]],
                               reduced_dims[, idx_conn] %>% t())))
      d = d[1, -1]
      nns = rbind(nns,
                  data.frame(distance=min(d), name=names(which.min(d))))
    }
    closest_node = nns$name[which.min(nns$distance)]
    new_state[i] = cds$State[colnames(cds) == closest_node]
  } else {
    conn_lengths = lapply(nn_states, function(x) length_branches[[as.character(x)]]) %>% unlist()
    new_state[i] = nn_states[which.max(conn_lengths)]
  }
}

changed_states = cbind(old_state, new_state)
tmp = lapply(seq_len(nrow(changed_states)), function(i) sort(changed_states[i, ])) %>%
  Reduce(rbind, .)
changed_states = changed_states[!duplicated(tmp), ]

cds$State_merged = cds$State
cds$State_merged = plyr::mapvalues(cds$State_merged, changed_states[, 1], changed_states[, 2])

p2 = plot_cell_trajectory(cds, color_by = factor(cds$State_merged)) + coord_fixed()

# 2. Merge branches that do not have bifurcation ===============================

merge_branches_nobi = function(cds, adj_mat, is_branch, old_state, new_state) {
  
  old_state_2 = c()
  new_state_2 = c()
  
  unique_state_2 = sort(unique(cds$State_merged))
  
  conn_branch_states_new = lapply(is_branch, function(x) {
    sort(unique(cds$State_merged[adj_mat[x, ] != 0]))
  })
  
  for (i in seq_len(length(new_state))) {
    
    # Find the connected states to the old and new
    conn_ = list()
    for (j in seq_len(length(conn_branch_states_new))) {
      x = conn_branch_states_new[[j]]
      if (new_state[i] %in% x) {
        conn_ = c(conn_, list(x))
      }
    }
    
    # Merge only with states that are neighbours and have a single connection
    candidates = lapply(conn_, function(x) {
      if (length(x) == 2) {
        return(x[x != unique_state_2[i]])
      }
    }) %>% unlist()
    candidates = unique(candidates)
    if (length(candidates) == 0) {
      next()
    }
    # check that the candidates are not present in a 3 connection
    
    for (j in seq_along(candidates)) {
      
      max_nn = max(lapply(conn_, function(x) {
        length(x)
      }) %>% unlist())
      
      if (max_nn == 2) {
        good_conn_ = lapply(conn_, function(x) {
          if (candidates[j] %in% x) {
            return(x)
          } else {
            return(NA)
          }
        })
        good_conn_ = good_conn_[!is.na(good_conn_)]
        for (j in 1:length(good_conn_)) {
          old_state_2 = c(old_state_2, good_conn_[[j]][1])
          new_state_2 = c(new_state_2, good_conn_[[j]][2])
        }
      }
    }
  }
  
  changed_states = cbind(old_state_2, new_state_2)
  tmp = lapply(seq_len(nrow(changed_states)), function(i) sort(changed_states[i, ])) %>%
    Reduce(rbind, .)
  changed_states = changed_states[!duplicated(tmp), ]

  unique_mapped = unique(c(changed_states))
  occurr_mapped = sapply(unique_mapped, function(x) sum(changed_states == x))
  mapped_states = cds$State_merged
  
  mapped_ = changed_states
  mapped_list = list()

  while(nrow(mapped_) > 0) {
    
    vals = mapped_[1, ]
    mapped_ = mapped_[-1, , drop = F]

    found = FALSE
    
    if (length(mapped_list) > 0) {
      for (k in 1:length(mapped_list)) {
        if (any(vals %in% mapped_list[[k]])) {
          mapped_list[[k]] = c(mapped_list[[k]], vals)
          found = TRUE
        }
      }
    }

    if (!found) {
      mapped_list = c(mapped_list, list(vals))
    }
  }
  
  for (i in seq_along(mapped_list)) {
    mapped_states[mapped_states %in% mapped_list[[i]]] = max(mapped_states) + i
  }
  
  return(mapped_states)
  
}

new_states = merge_branches_nobi(cds, adj_mat, is_branch, old_state, new_state)
table(new_states)

cds$State_merged = new_states
cds$State_merged = as.factor(cds$State_merged)

# 3. Merge branches with less than 5 nodes =====================================

cds$State_merged = as.numeric(as.character(cds$State_merged))

unique_values_branches = sort(unique(cds$State_merged))
branch_size = sapply(unique_values_branches, function(x) sum(cds$State_merged == x))
small_branches = as.numeric(unique_values_branches)[branch_size < 5]

mean_branch_pos = lapply(unique_values_branches, function(x) {
  mean_coords = apply(cds@reducedDimS[, cds$State_merged == x, drop = FALSE], 1, mean)
}) %>%
  Reduce(rbind, .) %>%
  set_rownames(unique_values_branches)

dist_between_branches = dist(mean_branch_pos) %>% as.matrix()
# Find the closest branch to the small ones
merged_branches = c()
for (small_branch in small_branches) {
  ord = order(dist_between_branches[unique_values_branches == small_branch, ])
  for (i in 2:length(ord)) {
    if (!ord[i] %in% small_branches) {
      merged_branches = rbind(merged_branches,
                              c(small_branch, unique_values_branches[ord[i]]))
      break()
    }
  }
}

cds$State_merged = plyr::mapvalues(cds$State_merged, merged_branches[, 1], merged_branches[, 2])
cds$State_merged = factor(cds$State_merged)

p3 = plot_cell_trajectory(cds, color_by = factor(cds$State_merged)) + coord_fixed()

# 4. Hierarchy of the tree =====================================================

# Find the hierarchy and show the merged states as the second level starting
# from the center node

tree_centre_node = find_tree_centre_node(tree_backbone)
root_state = cds$State_merged[colnames(cds) == tree_centre_node]
message("Root state: ", root_state)

p4 = plot_complex_cell_trajectory(cds, root_states = root_state, color_by = "State_merged")

p_combo = cowplot::plot_grid(plotlist = list(p1, p2, p3, p4))

# Save =========================================================================

ggsave(filename = glue::glue("../plots/hierarchy_merged_branches_{frame}_new.pdf"),
       dpi = 300, width = 18, height = 8, plot = p_combo)

saveRDS(cds, file = glue::glue("{results_root}/cds_{frame}.rds"))