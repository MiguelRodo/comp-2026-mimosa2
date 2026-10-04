# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# 24_sim_maxit_check.R: Study 5 - justification of maxit=30
# ==============================================================================
# Uses the BASELINE design and its stored random-number streams, so the
# datasets are exactly the baseline datasets with the same Task_ID. 
# Each is fitted with maxit = 30 (as in every study) and maxit = 100 (package
# default). 
# The analysis (analysis/35_maxit_check.R) checks
#   (1) MIMOSA2_maxit30 results equal the baseline results (reproducibility),
#   (2) how much AUC / calls change with maxit = 100.
# ==============================================================================
source("sims/_header.R")
STUDY <- "maxit_check"

scen_b   <- scenarios_baseline()
design_b <- build_design("baseline", scen_b)            # same streams as the baseline study
design   <- design_b[maxit_check_subset(design_b) & design_b$Rep <= NSIM_MAX[["maxit_check"]], ]
design$In_standard <- TRUE
design$In_smoke    <- design$P == 20 & design$Effect == 2.5e-4 & design$Cell_range == "Medium"
message(sprintf("maxit check: %d scenarios; this profile runs %d datasets",
                length(unique(design$Scenario_ID)), nrow(select_profile(design, STUDY))))

task_fun <- function(row) {
  sim <- simulate_baseline(row)
  P <- length(sim$ns1)
  analyse_dataset(sim, task_id = row$Task_ID,
                  mimosa_subsets = list(MIMOSA2_maxit30 = seq_len(P), MIMOSA2_maxit100 = seq_len(P)),
                  maxit = c(MIMOSA2_maxit30 = 30, MIMOSA2_maxit100 = 100))
}

run_study(STUDY, design, task_fun, n_fits = 2)
if (N_CHUNKS == 1) combine_study(STUDY)
