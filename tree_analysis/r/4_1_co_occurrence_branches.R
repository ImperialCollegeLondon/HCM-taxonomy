# Calculate the co-occurrence of the DDRTree branches in the ED and ES phases
#
# Author: Paolo Inglese, Imperial College London, 2022


# Libraries ====================================================================

library(tidyverse)
library(magrittr)
library(here)

# Setup ========================================================================

results_root <- here("../results")
plots_root <- here("../plots")
data_root <- here("../data")

# Functions ====================================================================

jaccard_index <- function(x, y) {
  i_x <- which(x == 1)
  i_y <- which(y == 1)
  return(length(intersect(i_x, i_y)) / length(union(i_x, i_y)))
}

# Start ========================================================================

branches_ED <- read_csv(file.path(results_root, "ddrtree_res_ED_decim0.99.csv")) %>%
  select(ID, branch_merged) %>%
  rename("branch_ED" = "branch_merged")
branches_ES <- read_csv(file.path(results_root, "ddrtree_res_ES_decim0.99.csv")) %>%
  select(ID, branch_merged) %>%
  rename("branch_ES" = "branch_merged") %>%
  mutate(branch_ES = plyr::mapvalues(.$branch_ES, c("4", "5"), c("5", "4"))) %>%
  mutate(branch_ES = fct_relevel(branch_ES, "1", "2", "3", "4", "5"))

branches <- inner_join(branches_ED, branches_ES, by = "ID")

ub_ed <- sort(unique(branches$branch_ED))
ub_es <- sort(unique(branches$branch_ES))

co_occ_mat <- matrix(NA, length(ub_ed), length(ub_es),
                     dimnames = list(paste0("ED_", ub_ed), paste0("ES_", ub_es)))
for (i in seq_along(ub_ed)) {
  for (j in seq_along(ub_es)) {
    co_occ_mat[i, j] <- jaccard_index(branches$branch_ED == ub_ed[i],
                                      branches$branch_ES == ub_es[j])
  }
}

pdf(file = file.path(plots_root, "co_occurrence_branches.pdf"))
corrplot::corrplot(co_occ_mat, order = "original", addCoef.col = 1)
dev.off()