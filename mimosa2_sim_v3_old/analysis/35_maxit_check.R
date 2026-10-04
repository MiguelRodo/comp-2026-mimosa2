# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# 35_maxit_check.R: Study 5 - justification of maxit=30 + reproducibility 
# ==============================================================================
# (1) REPRODUCIBILITY  
#     The maxit-check study regenerated the baseline datasets from their stored 
#     streams. Its maxit = 30 fit must give the same AUC as the baseline study for 
#     the same Task_ID. Any difference means the datasets were not reproduced 
#     exactly.
# (2) MAXIT. 
#     Paired differences (maxit 100 - maxit 30) in per-dataset AUC and
#     TPR, agreement of the responder calls, and fit time. If the differences
#     are within Monte Carlo error, maxit = 30 is justified (Chapter 3 TODO).
# ==============================================================================
source("analysis/analysis_functions.R")
res <- load_results("maxit_check")

# --------------------------- (1) reproducibility ------------------------------
base <- tryCatch(load_results("baseline"), error = function(e) NULL)
if (!is.null(base)) {
  a <- res$estimates %>% filter(Method == "MIMOSA2_maxit30", Alpha == ALPHAS[1], Status == "ok") %>% select(Task_ID, AUC30 = AUC, TP30 = TP)
  b <- base$estimates %>% filter(Method == "MIMOSA2", Alpha == ALPHAS[1], Status == "ok") %>% select(Task_ID, AUCb = AUC, TPb = TP)
  ab <- inner_join(a, b, by = "Task_ID")
  ok <- nrow(ab) > 0 && isTRUE(all.equal(ab$AUC30, ab$AUCb, tolerance = 1e-8)) && all(ab$TP30 == ab$TPb)
  message(sprintf("Reproducibility: %d datasets in both studies; identical AUC and TP: %s", nrow(ab), ok))
  if (!ok && nrow(ab) > 0) warning("maxit=30 refits differ from the baseline results: datasets were NOT reproduced exactly.")
} else message("Baseline results not found: reproducibility check skipped.")

# --------------------------- (2) maxit 100 vs 30 ------------------------------
dm <- res$estimates %>% add_dataset_measures() %>% filter(Rule == "BFDR") %>%
  select(Task_ID, Method, Alpha, AUC, TPR, FDP, Status) %>%
  pivot_wider(names_from = Method, values_from = c(AUC, TPR, FDP, Status)) %>%
  filter(Status_MIMOSA2_maxit30 == "ok", Status_MIMOSA2_maxit100 == "ok") %>%
  mutate(dAUC = AUC_MIMOSA2_maxit100 - AUC_MIMOSA2_maxit30,
         dTPR = TPR_MIMOSA2_maxit100 - TPR_MIMOSA2_maxit30,
         dFDP = FDP_MIMOSA2_maxit100 - FDP_MIMOSA2_maxit30) %>%
  left_join(res$design %>% select(Task_ID, Scenario_ID, P, Effect, Cell_range), by = "Task_ID")
calls <- res$subjects %>% filter(Eval) %>%
  mutate(agree01 = (MIMOSA2_maxit30_q < 0.01) == (MIMOSA2_maxit100_q < 0.01),
         dprob = abs(MIMOSA2_maxit100_prob - MIMOSA2_maxit30_prob)) %>%
  group_by(Task_ID) %>% summarise(call_agreement = mean(agree01, na.rm = TRUE), max_dprob = suppressWarnings(max(dprob, na.rm = TRUE)), .groups = "drop")
tm <- res$fits %>% select(Task_ID, Method, Fit_time, Hit_maxit) %>%
  pivot_wider(names_from = Method, values_from = c(Fit_time, Hit_maxit))
summ <- dm %>% left_join(calls, by = "Task_ID") %>% left_join(tm, by = "Task_ID") %>%
  group_by(P, Effect, Cell_range, Alpha) %>%
  summarise(n = n(),
            dAUC_mcse = mc_se(dAUC), dAUC = mean(dAUC, na.rm = TRUE),   # MCSE first (same name re-used)
            dTPR_mcse = mc_se(dTPR), dTPR = mean(dTPR, na.rm = TRUE),
            dFDP_mcse = mc_se(dFDP), dFDP = mean(dFDP, na.rm = TRUE),
            call_agreement = mean(call_agreement, na.rm = TRUE),
            max_dprob = max(max_dprob, na.rm = TRUE),
            hit_maxit_30 = mean(Hit_maxit_MIMOSA2_maxit30, na.rm = TRUE),
            time_30 = median(Fit_time_MIMOSA2_maxit30, na.rm = TRUE),
            time_100 = median(Fit_time_MIMOSA2_maxit100, na.rm = TRUE), .groups = "drop")
save_tab(summ, "maxit_check_summary")
print(knitr::kable(summ, digits = 4, caption = "maxit = 100 minus maxit = 30 on the same datasets (mean, MCSE)"))
message("Interpretation: if |dAUC| and |dTPR| are within ~2 MCSE of 0 and call agreement is ~1,")
message("then stopping at 30 EM iterations does not change the conclusions.")
