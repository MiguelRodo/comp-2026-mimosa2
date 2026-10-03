# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# 22_sim_heterogeneity.R: Study 3 - Effect-size heterogeneity
# Replaces Heterogeneous_effect_sim.R 
# ==============================================================================
# Summary of what changed relative to Heterogeneous_effect_sim.R:
# - BUG: the effect sizes of the two subgroups were swapped, so the evaluated
#   "small" group actually had the LARGER effect (scenarios.R, Study 3).
# - Version 2 called MIMOSA2() with no error handling or time limit, so one
#   failing fit lost the whole task. Both fits now have FIT_TIMEOUT and their
#   status is recorded.
# - All subjects (focal and other group) are stored with a Group column, so
#   the effect of pooling on the strong responders can also be looked at.
# - DiD calls with Benjamini-Hochberg are adjusted within the focal group.
# - nsim = 100 in both profiles (each task has two MIMOSA2 fits).
# =============================================================================
source("sims/_header.R")
STUDY <- "heterogeneity"

scen   <- scenarios_heterogeneity()
design <- build_design(STUDY, scen)
message(sprintf("Heterogeneity: %d scenarios; this profile runs %d datasets",
                nrow(scen), nrow(select_profile(design, STUDY))))

task_fun <- function(row) {
  sim <- simulate_heterogeneity(row)
  P <- row$P_focal + row$P_other
  focal <- seq_len(row$P_focal)
  analyse_dataset(sim, task_id = row$Task_ID,
                  eval_idx = focal,
                  mimosa_subsets = list(MIMOSA2_pooled = seq_len(P),     # fitted to all subjects
                                        MIMOSA2_separate = focal),       # fitted to the focal group only
                  group = rep(c("focal", "other"), c(row$P_focal, row$P_other)))
}

run_study(STUDY, design, task_fun, n_fits = 2)
if (N_CHUNKS == 1) combine_study(STUDY)
