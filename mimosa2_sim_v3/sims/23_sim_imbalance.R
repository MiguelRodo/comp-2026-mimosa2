# =============================================================================
# 23_sim_imbalance.R : Study 4, cell-count imbalance
# Replaces Count_Imbalanced_Sim.R (version 2).
# See STUDY_DESCRIPTIONS.md (Study 4).
# =============================================================================
# [CHANGE] Summary of what changed relative to Count_Imbalanced_Sim.R
#  * Scenarios are generated from (depleted assays, fraction) rather than a
#    hand-typed list of 120 vectors (scenarios.R); the depleted assays and the
#    matched balanced control (same total cells) are stored as columns, so the
#    analysis no longer has to parse names. (Count_Imbalanced_Plot.R looked
#    for names ending in "_Imb"/"_Bal", but the names end in "_Imb_s1",
#    "_Bal_1", ..., so every scenario was labelled "Balanced" in the ROC plot
#    and every non-reference scenario was dropped from the AUC odds-ratio plot.)
#  * set.seed() with forked workers was not reproducible; now streams.
#  * sim_progress.log shows 450 tasks (15 scenarios x 30), so the saved
#    Count_Imbalance.Rdata came from an older 15-scenario version of the
#    script, not from Count_Imbalanced_Sim.R (audit A8).
# =============================================================================
source("sims/_header.R")
STUDY <- "imbalance"

scen   <- scenarios_imbalance()
design <- build_design(STUDY, scen)
message(sprintf("Imbalance: %d scenarios; this profile runs %d datasets",
                nrow(scen), nrow(select_profile(design, STUDY))))

task_fun <- function(row) analyse_dataset(simulate_imbalance(row), task_id = row$Task_ID)

run_study(STUDY, design, task_fun, n_fits = 1)
if (N_CHUNKS == 1) combine_study(STUDY)
