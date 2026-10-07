# =============================================================================
# 29_combine.R : combine task files into one results file per study.
# Only needed after a SLURM job ARRAY (several chunks); single-chunk runs
# combine automatically at the end. Usage:
#   SIM_PROFILE=standard Rscript sims/29_combine.R baseline prior heterogeneity imbalance maxit_check
# =============================================================================
for (f in c("R/config.R", "R/runner.R")) source(f)
studies <- commandArgs(trailingOnly = TRUE)
if (length(studies) == 0 && nzchar(Sys.getenv("STUDIES"))) studies <- strsplit(Sys.getenv("STUDIES"), "[ ,]+")[[1]]   # set by slurm/combine.sh
if (length(studies) == 0) studies <- c("baseline", "prior", "heterogeneity", "imbalance", "maxit_check")
for (s in studies) {
  if (!file.exists(file.path(OUT_DIR, s, "design_full.rds"))) { message("skip ", s, " (not run)"); next }
  tryCatch(combine_study(s), error = function(e) message("NOT combined: ", conditionMessage(e)))   # [v3.1]
}
