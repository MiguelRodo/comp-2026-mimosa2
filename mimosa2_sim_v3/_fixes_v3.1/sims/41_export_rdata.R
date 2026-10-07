# =============================================================================
# 41_export_rdata.R : save a copy of every results file as .Rdata  (NEW v3.1)
# =============================================================================
# The pipeline itself keeps using .rds (every script reads <study>_results_<profile>.rds),
# so this only ADDS a .Rdata copy next to each .rds file; nothing is changed or removed.
#
# Each .Rdata file contains:
#   <study>_results          the whole results list (design, datasets, fits, estimates,
#                            subjects, task_time, run_info, ...), identical to the .rds
#   <study>_design, <study>_datasets, <study>_fits, <study>_estimates, <study>_subjects
#                            the main tables as separate data frames, for convenience
# The names include the study, so several files can be load()-ed into one session.
#
# Usage (from the mimosa2_sim_v3 folder):
#   SIM_PROFILE=standard Rscript sims/41_export_rdata.R            # all studies
#   SIM_PROFILE=standard Rscript sims/41_export_rdata.R baseline   # one study
# In R / RStudio:  load("_simulations/baseline/baseline_results_standard.Rdata"); ls()
# =============================================================================
source("R/config.R")
studies <- commandArgs(trailingOnly = TRUE)
if (length(studies) == 0 && nzchar(Sys.getenv("STUDIES"))) studies <- strsplit(Sys.getenv("STUDIES"), "[ ,]+")[[1]]
if (length(studies) == 0) studies <- c("baseline", "prior", "heterogeneity", "imbalance", "maxit_check")

for (st in studies) {
  f_rds <- file.path(OUT_DIR, st, sprintf("%s_results_%s.rds", st, PROFILE))
  if (!file.exists(f_rds)) { message("skip ", st, ": ", f_rds, " not found"); next }
  res <- readRDS(f_rds)
  env <- new.env()
  assign(paste0(st, "_results"), res, envir = env)
  for (tb in c("design", "datasets", "fits", "estimates", "subjects"))
    if (!is.null(res[[tb]])) assign(paste0(st, "_", tb), res[[tb]], envir = env)
  f_rdata <- sub("\\.rds$", ".Rdata", f_rds)
  tmp <- paste0(f_rdata, ".tmp")
  save(list = ls(env), envir = env, file = tmp, compress = TRUE)
  file.rename(tmp, f_rdata)
  message(sprintf("[%s] saved %s (objects: %s)", st, f_rdata, paste(sort(ls(env)), collapse = ", ")))
}
