# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# 20_sim_baseline.R: Study 1 - Baseline sensitivity analysis 
# Replaces Simulation.R
# ==============================================================================
# Summary of what changed relative to Simulation.R:
# - Simulation.R could not run as written: it put `did_glm_shrink_prob` (never
#   defined) into the results data frame, and used withTimeout() without
#   loading R.utils. It therefore cannot be the script that produced
#   Simulation_2.0.Rdata (audit A8). This script reproduces the same design.
# - Seed set once (L'Ecuyer-CMRG) + one stream per dataset, stored (runner.R).
#   Responder allocation fixed (dgm.R CHANGE 1): realised rho now equals the
#   nominal rho up to rounding of rho * P (stored per dataset).
# - Null scenarios (rho = 0) added as a check with known properties.
#   nsim = 100 (standard) / 200 (extended) instead of 30, chosen for MCSE.
#   Per-dataset TP/FP/TN/FN for MIMOSA2 (Bayesian FDR via getResponse) and
#   DiD (unadjusted and Benjamini-Hochberg) at alpha = 0.01 and 0.05, and the
#   per-dataset AUC of both methods, are stored (audit A3, A5).
# - Iterations / time / failure status of every MIMOSA2 fit are stored.
# - "Sparse" (15,000 cells) is called "Low", as in the other studies.
# ==============================================================================
source("sims/_header.R")
STUDY <- "baseline"

scen   <- scenarios_baseline()
design <- build_design(STUDY, scen)
message(sprintf("Baseline: %d scenarios; this profile runs %d datasets",
                nrow(scen), nrow(select_profile(design, STUDY))))

task_fun <- function(row) analyse_dataset(simulate_baseline(row), task_id = row$Task_ID)

run_study(STUDY, design, task_fun, n_fits = 1)
if (N_CHUNKS == 1) combine_study(STUDY)
