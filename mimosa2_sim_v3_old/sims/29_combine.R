# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# 29_combine.R: Combines task files int o one reuslts file/study
# ==============================================================================
# Only needed after a SLURM job ARRAY (several chunks)
# Single-chunk runs combine automatically at the end. 
# Usage:
#   SIM_PROFILE=standard 
#   Rscript sims/29_combine.R baseline prior heterogeneity imbalance maxit_check
# ==============================================================================
for (f in c("R/config.R", "R/runner.R")) source(f)
studies <- commandArgs(trailingOnly = TRUE)
if (length(studies) == 0) studies <- c("baseline", "prior", "heterogeneity", "imbalance", "maxit_check")
for (s in studies) {
  if (!file.exists(file.path(OUT_DIR, s, "design_full.rds"))) { message("skip ", s, " (not run)"); next }
  combine_study(s)
}
