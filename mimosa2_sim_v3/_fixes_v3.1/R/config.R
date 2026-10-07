# =============================================================================
# config.R : one place for paths, seeds, run profile and fitting settings
# MIMOSA2 simulation study, Isabella Lethbridge & Tayyeb Abrahams
# Version 3 (re-run version, October 2026)
# =============================================================================
#
# [CHANGE] NEW FILE. In version 2, library paths, seeds, time limits and the
# number of workers were hard-coded separately in every script, and they
# differed between scripts:
#   * the time limit was 600 s (baseline), 120 s (prior, earlier version),
#     180 s (prior, imbalance) and none (heterogeneity);
#   * the baseline script set no seed at all;
#   * the cell-count ranges, phi and P were copied by hand.
# Morris et al. (2019, Table 1 and Section 4) ask for the code to reproduce the
# study exactly, so every setting that affects the results now lives here and
# is saved alongside the results (runner.R -> run_info_chunkNN.rds, also copied into
# the combined results file as $run_info).
#
# HOW TO CHOOSE A PROFILE (set the environment variable SIM_PROFILE):
#   SIM_PROFILE=smoke     tiny run (minutes) to test the whole pipeline
#   SIM_PROFILE=standard  sized for ~20 cores: ~375 core-hours (~19 h on 20 cores)
#   SIM_PROFILE=extended  full grids, larger nsim: ~1650 core-hours, for 64+ cores or
#                         several nodes (SLURM job array)
# Datasets are identical across profiles: the standard run is an exact subset
# of the extended run (same task IDs, same random-number streams), so if you
# run "standard" first and later get more compute, running "extended" only
# computes the tasks that are still missing.
# =============================================================================

# ---- Library path on the cluster (unchanged from your scripts) --------------
# [CHANGE 4 Oct] Exactly your two lines, applied unconditionally and FIRST.
# The path is also exported as R_LIBS so that every child R process started
# from here (the smoke test runs each script with Rscript) uses it too.
my_libs <- Sys.getenv("MIMOSA2_R_LIBS", "/scratch/abrmoe030/R_libs")
if (!dir.exists(my_libs)) warning("Library folder not found on this machine: ", my_libs)
.libPaths(c(my_libs, .libPaths()))
Sys.setenv(R_LIBS = paste(unique(c(my_libs, strsplit(Sys.getenv("R_LIBS"), ":")[[1]])), collapse = ":"))
MY_LIBS <- my_libs
message("[config] R library paths: ", paste(.libPaths(), collapse = " | "))

# ---- Run profile -------------------------------------------------------------
PROFILE <- tolower(Sys.getenv("SIM_PROFILE", "standard"))
stopifnot(PROFILE %in% c("smoke", "standard", "extended"))

# ---- Output locations (relative to the folder you run Rscript from) ---------
OUT_DIR   <- Sys.getenv("MIMOSA2_OUT_DIR", "_simulations")
FIG_DIR   <- Sys.getenv("MIMOSA2_FIG_DIR", "_fig")
TAB_DIR   <- Sys.getenv("MIMOSA2_TAB_DIR", "_tables")
for (d in c(OUT_DIR, FIG_DIR, TAB_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- Parallel workers --------------------------------------------------------
# [CHANGE] Reads N_WORKERS first (set by the slurm/ scripts), then
# SLURM_CPUS_PER_TASK, then SLURM_NTASKS (what you used before). All workers must be on ONE node because the
# scheduler forks processes (parallel::mcparallel). To use several nodes,
# submit a SLURM job array: each array element runs its own share of tasks.
.first_num <- function(...) {
  for (v in c(...)) { x <- suppressWarnings(as.numeric(Sys.getenv(v, ""))); if (!is.na(x) && x > 0) return(x) }
  NA_real_
}
N_WORKERS <- .first_num("N_WORKERS", "SLURM_CPUS_PER_TASK", "SLURM_NTASKS")
if (is.na(N_WORKERS)) N_WORKERS <- 2

# SLURM job array: chunk k of K (k = 0..K-1). Without an array: one chunk.
CHUNK_ID <- .first_num("SLURM_ARRAY_TASK_ID"); if (is.na(CHUNK_ID)) CHUNK_ID <- 0
N_CHUNKS <- .first_num("N_CHUNKS", "SLURM_ARRAY_TASK_COUNT"); if (is.na(N_CHUNKS)) N_CHUNKS <- 1
# If your array does not start at 0 (e.g. --array=1-4), subtract the minimum:
ARRAY_MIN <- .first_num("SLURM_ARRAY_TASK_MIN"); if (!is.na(ARRAY_MIN)) CHUNK_ID <- CHUNK_ID - ARRAY_MIN

# ---- MIMOSA2 fitting settings (identical in every study) ---------------------
MAXIT       <- 30    # as in version 2 (package default is 100); see sims/24_sim_maxit_check.R
FIT_TIMEOUT <- 600   # [CHANGE] ONE time limit (seconds) per MIMOSA2 fit in every study.
                     # Version 2 used 600/180/120/none. A fit that runs longer is
                     # killed and recorded as a "timeout" failure (Morris 5.1).
ALPHAS      <- c(0.01, 0.05)  # nominal levels for responder calls
# [CHANGE v3.1, 6 Oct] Continuity correction used in the VARIANCE of the DiD
# Wald test: p_tilde = (x + DID_CC) / (N + 2 * DID_CC). With DID_CC = 0 (v3.0)
# an assay with a zero count contributed zero variance, which made the test
# badly anti-conservative when s0 or u1 had few cells (imbalance study: 26% of
# non-responders called at alpha = 0.05). The estimate itself is unchanged.
DID_CC      <- 0.5

# ---- Seeds: ONE master seed per study (Morris 4.1: set the seed once) --------
# [CHANGE] Version 2: baseline had no seed; the others called set.seed() but
# used the default Mersenne-Twister generator with forked workers, so a run
# was not reproducible. Now: L'Ecuyer-CMRG generator, a master seed per study,
# and an independent random-number STREAM for every dataset (Morris 4.1.1).
# The seeds of the old scripts are kept where they existed.
MASTER_SEED <- c(
  baseline      = 20260626,  # new (Simulation.R had no seed)
  prior         = 9854,      # as Prior_simulations.R
  heterogeneity = 4471,      # as Heterogeneous_effect_sim.R
  imbalance     = 3857,      # as Count_Imbalanced_Sim.R
  maxit_check   = 20260626   # re-uses the baseline streams on purpose (same datasets)
)

# ---- Number of repetitions per scenario (nsim) by profile --------------------
# [CHANGE] nsim is now justified by the Monte Carlo SE it buys (Morris 5.3):
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
# Largest nsim that any profile uses. Seeds are generated for this many
# repetitions so that every profile draws the same datasets.
NSIM_MAX <- c(baseline = 200, prior = 200, heterogeneity = 100, imbalance = 200, maxit_check = 100)
for (p in names(NSIM)) stopifnot(all(NSIM[[p]] <= NSIM_MAX[names(NSIM[[p]])]))

nsim_for <- function(study) unname(NSIM[[PROFILE]][study])

# ---- Shared DGM constants (as in version 2) ----------------------------------
BASELINE_BACKGROUND  <- 1e-4    # mu_0
BASELINE_STIM_EFFECT <- 2.5e-4  # delta_0
BG_EFFECT            <- 0       # gamma

message(sprintf("[config] profile=%s workers=%d chunk=%d/%d maxit=%d timeout=%ds",
                PROFILE, N_WORKERS, CHUNK_ID + 1, N_CHUNKS, MAXIT, FIT_TIMEOUT))
