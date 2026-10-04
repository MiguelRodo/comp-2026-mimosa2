# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# config.R
# - paths
# - seeds 
# - run profiles
# - fitting settings 
# ==============================================================================
# Settings have been updated to ensure reproducibility, Morris et al. (2019, Table 1 and Section 4)
# Settings and results saved (runner.R -> run_info_chunkNN.rds, also copied into the combined results file as $run_info)
# Datasets are identical across profiles
#    Standard run is exact subset of extended run (same task IDs, same random-number streams).
#    If you run standard first then extended, extended only computes missidng tasks
# ==============================================================================
# Profile Specification: 
#   SIM_PROFILE=smoked:   test pipeline (short run)
#   SIM_PROFILE=standard: size for ~20 cores = ~375 core-hours (~19h on 20 cores)
#   SIM_PROFILE=extended: larger nsim, ~1650 core-hours for 64+ cores/several nodes  (SLURM job array)
# ==============================================================================

# ------------------------ Library path on the cluster -------------------------
MY_LIBS <- Sys.getenv("MIMOSA2_R_LIBS", "/scratch/abrmoe030/R_libs")
if (dir.exists(MY_LIBS)) .libPaths(c(MY_LIBS, .libPaths()))

# -------------------------------- Run profile ---------------------------------
PROFILE <- tolower(Sys.getenv("SIM_PROFILE", "standard"))
stopifnot(PROFILE %in% c("smoke", "standard", "extended"))

# ----------------------------- Output locations -------------------------------
# Relative to folder from which Rscript was run from 
OUT_DIR   <- Sys.getenv("MIMOSA2_OUT_DIR", "_simulations")
FIG_DIR   <- Sys.getenv("MIMOSA2_FIG_DIR", "_fig")
TAB_DIR   <- Sys.getenv("MIMOSA2_TAB_DIR", "_tables")
for (d in c(OUT_DIR, FIG_DIR, TAB_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ------------------------------ Parallel workers ------------------------------
# Reads N_WORKERS first (set by slurm/ scripts), then SLURM_CPUS_PER_TASK, then SLURM_NTASKS. 
# All workers must be on ONE node because scheduler forks processes (parallel::mcparallel). 
# To use several nodes submit a SLURM job array: each array element runs its own share of tasks.
.first_num <- function(...) {
  for (v in c(...)) { x <- suppressWarnings(as.numeric(Sys.getenv(v, ""))); if (!is.na(x) && x > 0) return(x) }
  NA_real_
}
N_WORKERS <- .first_num("N_WORKERS", "SLURM_CPUS_PER_TASK", "SLURM_NTASKS")
if (is.na(N_WORKERS)) N_WORKERS <- 2

# SLURM job array:  chunk k of K (k = 0..K-1). 
# Without an array: one chunk.
CHUNK_ID <- .first_num("SLURM_ARRAY_TASK_ID"); if (is.na(CHUNK_ID)) CHUNK_ID <- 0
N_CHUNKS <- .first_num("N_CHUNKS", "SLURM_ARRAY_TASK_COUNT"); if (is.na(N_CHUNKS)) N_CHUNKS <- 1

# If  array does not start at 0 (e.g. --array=1-4), subtract minimum:
ARRAY_MIN <- .first_num("SLURM_ARRAY_TASK_MIN"); if (!is.na(ARRAY_MIN)) CHUNK_ID <- CHUNK_ID - ARRAY_MIN

# --------------------------- MIMOSA2 fitting settings -------------------------
# Identical across every study 
MAXIT       <- 30    # (package default is 100); see sims/24_sim_maxit_check.R
FIT_TIMEOUT <- 600   # ONE time limit (seconds) per MIMOSA2 fit in every study.
    
ALPHAS      <- c(0.01, 0.05)  # nominal levels for responder calls

# -------------------- Seeds: ONE master seed per study ------------------------
# L'Ecuyer-CMRG generator, a master seed per study,and an independent random-number STREAM for every dataset (Morris 4.1.1).
# The seeds of the old scripts are kept where they existed.
MASTER_SEED <- c(
  baseline      = 20260626,  # new (Simulation.R had no seed)
  prior         = 9854,      # as Prior_simulations.R
  heterogeneity = 4471,      # as Heterogeneous_effect_sim.R
  imbalance     = 3857,      # as Count_Imbalanced_Sim.R
  maxit_check   = 20260626   # re-uses the baseline streams on purpose (same datasets)
)

# ----------- Number of repetitions per scenario (nsim) by profile -------------
# nsim is now justified by the Monte Carlo SE it buys (Morris 5.3):
#   MCSE(mean of a per-dataset proportion) = SD / sqrt(nsim) <= 0.5 / sqrt(nsim)
#   nsim = 100 -> worst case 0.050, typical (SD ~ 0.2) 0.020
#   nsim = 200 -> worst case 0.035, typical (SD ~ 0.2) 0.014
#   per-dataset AUC (SD ~ 0.05-0.15): MCSE 0.005-0.015 at nsim = 100
# The analysis scripts report the achieved MCSE and the nsim needed for a
# target MCSE, so the choice can be checked after the run (Morris 5.3).
NSIM <- list(
  smoke    = c(baseline = 2,   prior = 2,   heterogeneity = 2,   imbalance = 2,   maxit_check = 2),
  standard = c(baseline = 100, prior = 100, heterogeneity = 100, imbalance = 100, maxit_check = 50),
  extended = c(baseline = 200, prior = 200, heterogeneity = 100, imbalance = 200, maxit_check = 100)
)

# Largest nsim that any profile uses. 
# Seeds are generated for this many repetitions so that every profile draws the same datasets.
NSIM_MAX <- c(baseline = 200, prior = 200, heterogeneity = 100, imbalance = 200, maxit_check = 100)
for (p in names(NSIM)) stopifnot(all(NSIM[[p]] <= NSIM_MAX[names(NSIM[[p]])]))

nsim_for <- function(study) unname(NSIM[[PROFILE]][study])

# --------------------------- Shared DGM constants -----------------------------
BASELINE_BACKGROUND  <- 1e-4    # mu_0
BASELINE_STIM_EFFECT <- 2.5e-4  # delta_0
BG_EFFECT            <- 0       # gamma

message(sprintf("[config] profile=%s workers=%d chunk=%d/%d maxit=%d timeout=%ds",
                PROFILE, N_WORKERS, CHUNK_ID + 1, N_CHUNKS, MAXIT, FIT_TIMEOUT))
