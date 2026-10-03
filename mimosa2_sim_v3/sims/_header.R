# Shared start-up for every simulation script (run from the pipeline root folder):
#   cd mimosa2_sim_v3 ; SIM_PROFILE=standard Rscript sims/20_sim_baseline.R
# MIMOSA2 is fitted inside forked processes: stop BLAS/OpenMP from starting
# extra threads in every worker (that would oversubscribe the node).
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
for (f in c("R/config.R", "R/dgm.R", "R/methods.R", "R/runner.R", "R/scenarios.R")) {
  if (!file.exists(f)) stop("Run this script from the mimosa2_sim_v3 folder (cannot find ", f, ")")
  source(f)
}
suppressPackageStartupMessages(library(MIMOSA2))   # loaded once; forked workers inherit it
