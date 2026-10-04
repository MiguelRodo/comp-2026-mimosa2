# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# 30_performance_tables.R: Performance measures + MSE for every study
# ==============================================================================
# SIM_PROFILE=standard Rscript analysis/30_performance_tables.R
# Writes to _tables/:
#   <study>_performance.csv   one row per scenario x method x rule x alpha
#                             (TPR, TNR, FPR, FDR, P(any FP), AUC; each with MCSE;
#                             failure rate)
#   <study>_failures.csv      MIMOSA2 failures, timeouts, maxit, fit times
#   <study>_dAUC_*.csv        paired AUC differences with MCSE
#   <study>_realised_dgm.csv  realised rho and realised Delta per scenario
#   precision_report.csv      max MCSE per study/method and the nsim that a
#                             target MCSE would need (Morris 5.3)
# =============================================================================
source("analysis/analysis_functions.R")

studies <- c("baseline", "prior", "heterogeneity", "imbalance")
prec <- list()
for (st in studies) {
  res <- tryCatch(load_results(st), error = function(e) { message("skip ", st, ": ", conditionMessage(e)); NULL })
  if (is.null(res)) next
  message("\n==== ", st, " ====")
  perf <- performance_table(res)
  save_tab(perf, paste0(st, "_performance"))
  fails <- failure_table(res)
  save_tab(fails, paste0(st, "_failures"))
  fm <- fails %>% filter(grepl("^MIMOSA2", Method))
  message(sprintf("MIMOSA2 fits: %d planned datasets, failure rate %.2f%% (errors %d, timeouts %d); whole tasks failed: %d",
                  nrow(res$design), 100 * (1 - sum(fm$n_ok) / sum(fm$n_planned)),
                  sum(fm$n_error), sum(fm$n_timeout), sum(res$fits$Status == "task_failed")))
  # NB: MIMOSA2 does not return its EM iteration count, so Iter / Hit_maxit are NA;
  # whether maxit = 30 is enough is answered by Study 5 (35_maxit_check.R).
  # realised DGM quantities per scenario (rho, true Delta among responders)
  real <- res$datasets %>% left_join(res$design %>% select(Task_ID, Scenario_ID), by = "Task_ID") %>%
    group_by(Scenario_ID) %>%
    summarise(rho_realised = mean(rho_realised), n_resp_eval = mean(n_resp_eval),
              mean_true_delta_resp = mean(mean_true_delta_resp, na.rm = TRUE),
              share_zero_ns1 = mean(n_zero_ns1 / P_total), .groups = "drop") %>%
    left_join(distinct(res$design[, scenario_cols(res$design), drop = FALSE]), by = "Scenario_ID")
  save_tab(real, paste0(st, "_realised_dgm"))
  # paired AUC differences
  mims <- grep("^MIMOSA2", unique(res$fits$Method), value = TRUE)   # not the "ALL" rows of failed tasks
  for (m in mims) save_tab(paired_auc_diff(res, m, "DiD"), sprintf("%s_dAUC_%s_vs_DiD", st, m))
  if (all(c("MIMOSA2_pooled", "MIMOSA2_separate") %in% mims))
    save_tab(paired_auc_diff(res, "MIMOSA2_pooled", "MIMOSA2_separate"), sprintf("%s_dAUC_pooled_vs_separate", st))
  prec[[st]] <- precision_report(perf, st)
}
if (length(prec)) {
  pr <- bind_rows(prec)
  save_tab(pr, "precision_report")
  message("\nMonte Carlo precision achieved (max over scenarios) and nsim needed for")
  message("target MCSE of 0.02 (TPR), 0.01 (FDR), 0.01 (AUC):")
  print(as.data.frame(pr), digits = 3)
}
