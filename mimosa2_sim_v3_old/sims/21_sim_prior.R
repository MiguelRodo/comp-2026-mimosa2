# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
#21_sim_prio.R: Study 2 - Robustness to prior distribution 
# Replaces Prior_simulations.R
# ==============================================================================
# Summary of what changed relative to Prior_simulations.R:
# - Every non-Beta family is calibrated to the Beta(phi = 10000) mean AND
#   variance at EACH condition mean (dgm.R CHANGE 2-4). Version 2 matched the
#   variance only at mu_0 = 1e-4, so at the stimulated responder mean the EG
#   and LN families had 2-3x the Beta SD, and the LN mean was 1.28 mu.
#   Results from version 2 for EG and LN are therefore not comparable with
#   Beta and should not be reported as a pure "shape" effect.
# - The bimodal family is applied to all stimulated draws (CHANGE 4).
# - Version 2 grid had 200 replicates x 300 scenarios = 60,000 tasks (about
#   33 h on 20 cores) while the saved results had only 10 replicates
#   (500 rows per cell); now nsim = 100 (standard, P = 50 only) or 200
#   (extended, P = 20, 50, 100).
# - "Extremely Low" (3,000 cells) is now kept in the plots (it was dropped as
#   NA by Prior_simulation_plots.R, whose factor levels omitted it).
# - Seeds/streams, failure recording and per-dataset measures as baseline.
# =============================================================================
source("sims/_header.R")
STUDY <- "prior"

scen   <- scenarios_prior()
design <- build_design(STUDY, scen)
# Solve the EG / LN calibration ONCE in the parent process; workers inherit it.
n_cal <- prewarm_calibration(scen$Distribution, scen$Effect, scen$Phi)
message(sprintf("Prior family: %d scenarios; %d calibrated (family, mean) pairs; this profile runs %d datasets",
                nrow(scen), n_cal, nrow(select_profile(design, STUDY))))

task_fun <- function(row) analyse_dataset(simulate_prior(row), task_id = row$Task_ID)

run_study(STUDY, design, task_fun, n_fits = 1)
if (N_CHUNKS == 1) combine_study(STUDY)
