# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# analysis_functions.R: Performance measures with Monte Carlo Standard Errors
# ==============================================================================
# Used by analysis/30-35. 
# Run the analysis scripts from the mimosa2_sim_v3 folder, on the cluster or 
# locally after copying _simulations/ across:
#   SIM_PROFILE=standard Rscript analysis/30_performance_tables.R
# (In RStudio: setwd() to the mimosa2_sim_v3 folder and 
# Sys.setenv(SIM_PROFILE = "standard") before sourcing.)
#
# Key definitions:
# - All measures are computed PER DATASET and then averaged over the nsim
#   datasets of a scenario; the Monte Carlo SE of a mean is SD / sqrt(n). 
#   Scenarios are never pooled.
# - TPR (sensitivity) = TP / (TP + FN) in datasets with >= 1 responder.
# - TNR (specificity) = TN / (TN + FP); FPR = 1 - TNR.
# - FDP = FP / (TP + FP), defined as 0 when nothing is called; the FDR is the
#   mean FDP over datasets (this is the quantity the Bayesian FDR rule and
#   Benjamini-Hochberg aim to control).
# - P(any FP) = proportion of datasets with at least one false positive
#   (equals the FDR in the null scenarios, where every call is false).
# - AUC = per-dataset Mann-Whitney AUC (probability that a random responder
#   in the dataset outscores a random non-responder). Version 2
#   pooled subjects ACROSS datasets before computing one AUC, which mixes
#   posterior probabilities from different fitted models.
# - Delta AUC = paired difference AUC_MIMOSA2 - AUC_DiD in the same dataset;
#   its MCSE is SD(differences) / sqrt(n). Replaces the Hanley-McNeil intervals, 
#   which treated all pooled subjects as independent, fixed 50% positives and 
#   assumed r = 0.5.
# - Fit failures (error / timeout / task failed) are the first performance
#   measure and are reported per scenario. Other measures are computed on 
#   successful fits only.
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(ggplot2)
})
for (f in c("R/config.R", "R/scenarios.R")) source(f)

# ------------------------------------------------------------------------------
# Loading
# ------------------------------------------------------------------------------
load_results <- function(study, profile = PROFILE) {
  f <- file.path(OUT_DIR, study, sprintf("%s_results_%s.rds", study, profile))
  if (!file.exists(f)) stop("No results for ", study, " (", f, "). Run the simulation / sims/29_combine.R first.")
  readRDS(f)
}
scenario_cols <- function(design) setdiff(names(design), c("Rep", "Task_ID", "Seed", "In_standard", "In_smoke"))

# ------------------------------------------------------------------------------
# Monte Carlo SE helpers 
# ------------------------------------------------------------------------------
mc_n    <- function(x) sum(!is.na(x))
mc_mean <- function(x) if (mc_n(x) == 0) NA_real_ else mean(x, na.rm = TRUE)
mc_se   <- function(x) if (mc_n(x) < 2) NA_real_ else sd(x, na.rm = TRUE) / sqrt(mc_n(x))
nsim_needed <- function(sd, target) ceiling((sd / target)^2)       

# ------------------------------------------------------------------------------
# Per-dataset measures
# ------------------------------------------------------------------------------
add_dataset_measures <- function(est) {
  est %>% mutate(
    TPR    = ifelse(Status == "ok" & n_resp > 0, TP / (TP + FN), NA_real_),
    TNR    = ifelse(Status == "ok" & n_nonresp > 0, TN / (TN + FP), NA_real_),
    FPR    = 1 - TNR,
    Calls  = ifelse(Status == "ok", TP + FP, NA_real_),
    FDP    = ifelse(Status == "ok", ifelse(TP + FP > 0, FP / (TP + FP), 0), NA_real_),
    Any_FP = ifelse(Status == "ok", as.numeric(FP > 0), NA_real_),
    AUC    = ifelse(Status == "ok", AUC, NA_real_))
}

# ------------------------------------------------------------------------------
# Performance table: one row per scenario x method x rule x alpha
# ------------------------------------------------------------------------------
performance_table <- function(res) {
  sc <- scenario_cols(res$design)
  planned <- res$design %>% count(Scenario_ID, name = "n_planned")
  est <- res$estimates %>% add_dataset_measures() %>%
    left_join(res$design %>% select(Task_ID, Scenario_ID), by = "Task_ID")
  tab <- est %>%
    group_by(Scenario_ID, Method, Rule, Alpha) %>%
    # NB: inside summarise() a new column immediately replaces the old one, so
    # every SD / MCSE / n is computed BEFORE the mean that re-uses its name.
    summarise(n_ok = sum(Status == "ok"),
              TPR_mcse = mc_se(TPR), TPR_n = mc_n(TPR), TPR_sd = sd(TPR, na.rm = TRUE), TPR = mc_mean(TPR),
              TNR_mcse = mc_se(TNR), TNR = mc_mean(TNR),
              FPR_mcse = mc_se(FPR), FPR = mc_mean(FPR),
              FDR_mcse = mc_se(FDP), FDP_sd = sd(FDP, na.rm = TRUE), FDR = mc_mean(FDP),
              Any_FP_mcse = mc_se(Any_FP), Any_FP = mc_mean(Any_FP),
              Calls = mc_mean(Calls),
              AUC_mcse = mc_se(AUC), AUC_n = mc_n(AUC), AUC_sd = sd(AUC, na.rm = TRUE), AUC = mc_mean(AUC),
              .groups = "drop") %>%
    left_join(planned, by = "Scenario_ID") %>%
    mutate(Fail_rate = 1 - n_ok / n_planned,
           Fail_mcse = sqrt(Fail_rate * (1 - Fail_rate) / n_planned)) %>%
    left_join(distinct(res$design[, sc, drop = FALSE]), by = "Scenario_ID")
  # The mean of the per-dataset FDP over successful fits ignores datasets
  # whose fit failed; flag scenarios where that is non-trivial (Morris 5.1).
  tab %>% mutate(Fail_flag = Fail_rate > 0.01)
}

# Paired difference of per-dataset AUC between two methods (same datasets)
paired_auc_diff <- function(res, m1, m2) {
  a <- res$estimates %>% filter(Alpha == ALPHAS[1]) %>%
    filter(Status == "ok") %>%
    distinct(Task_ID, Method, AUC)
  w <- a %>% filter(Method %in% c(m1, m2)) %>%
    pivot_wider(names_from = Method, values_from = AUC) %>%
    { if (!all(c(m1, m2) %in% names(.))) .[0, ] %>% mutate(!!m1 := numeric(0), !!m2 := numeric(0)) else . } %>%
    filter(!is.na(.data[[m1]]), !is.na(.data[[m2]])) %>%
    mutate(dAUC = .data[[m1]] - .data[[m2]]) %>%
    left_join(res$design %>% select(Task_ID, Scenario_ID), by = "Task_ID")
  out <- w %>% group_by(Scenario_ID) %>%
    summarise(n_pairs = n(), dAUC_mean = mean(dAUC), dAUC_sd = sd(dAUC),
              dAUC_mcse = sd(dAUC) / sqrt(n()),
              corr_AUC = if (length(.data[[m1]]) < 3) NA_real_ else suppressWarnings(cor(.data[[m1]], .data[[m2]])), .groups = "drop") %>%
    mutate(Comparison = paste(m1, "-", m2),
           lower = dAUC_mean - 1.96 * dAUC_mcse, upper = dAUC_mean + 1.96 * dAUC_mcse) %>%
    left_join(distinct(res$design[, scenario_cols(res$design), drop = FALSE]), by = "Scenario_ID")
  out
}

# ------------------------------------------------------------------------------
# FDR calibration curve (per dataset, then averaged) 
# Version 2 thresholded each subject's LOCAL posterior (1 - P(R) <= a) and
# pooled subjects across datasets and P. The rule actually used is the
# cumulative Bayesian FDR (q-value) of each dataset, so here a subject is
# called at level a if q < a (exactly as MIMOSA2::getResponse), the FDP is computed per dataset, and FDR(a) is the mean FDP over datasets with its MCSE.
# ------------------------------------------------------------------------------
fdr_curve <- function(res, alphas = seq(0.01, 0.20, by = 0.01), prob_cols = NULL) {
  s <- res$subjects %>% filter(Eval)
  if (is.null(prob_cols)) prob_cols <- grep("^MIMOSA2.*_q$", names(s), value = TRUE)
  split_s <- split(s, s$Task_ID)
  one <- function(d) {
    out <- list()
    for (qc in prob_cols) {
      q <- d[[qc]]; if (all(is.na(q))) next
      calls <- outer(q, alphas, "<")      # getResponse() rule: q < alpha
      nc <- colSums(calls); fp <- colSums(calls & d$Truth == 0)
      out[[qc]] <- data.frame(Method = sub("_q$", "", qc), Alpha = alphas, FDP = ifelse(nc > 0, fp / nc, 0), Calls = nc)
    }
    p_bh <- p.adjust(d$DiD_p, "BH")
    calls <- outer(p_bh, alphas, "<=")
    nc <- colSums(calls); fp <- colSums(calls & d$Truth == 0)
    out[["DiD"]] <- data.frame(Method = "DiD (BH)", Alpha = alphas, FDP = ifelse(nc > 0, fp / nc, 0), Calls = nc)
    bind_rows(out) %>% mutate(Task_ID = d$Task_ID[1])
  }
  per_ds <- bind_rows(lapply(split_s, one))
  per_ds %>% left_join(res$design %>% select(Task_ID, Scenario_ID), by = "Task_ID") %>%
    group_by(Scenario_ID, Method, Alpha) %>%
    summarise(n = n(), FDR = mean(FDP), FDR_mcse = sd(FDP) / sqrt(n()), Calls = mean(Calls), .groups = "drop") %>%
    left_join(distinct(res$design[, scenario_cols(res$design), drop = FALSE]), by = "Scenario_ID")
}

# ------------------------------------------------------------------------------
# Vertically averaged ROC curves: TPR of each dataset at a grid of FPR values, 
# averaged over datasets.
# Version 2 drew ONE ROC curve per scenario from subjects pooled across datasets 
# (and across P).
# ------------------------------------------------------------------------------
roc_at_grid <- function(score, truth, grid) {
  ok <- !is.na(score); score <- score[ok]; truth <- truth[ok]
  n1 <- sum(truth == 1); n0 <- sum(truth == 0)
  if (n1 == 0 || n0 == 0) return(rep(NA_real_, length(grid)))
  o <- order(score, decreasing = TRUE); s <- score[o]; t <- truth[o]
  last <- c(s[-1] != s[-length(s)], TRUE)           # end of each block of tied scores
  tp <- c(0, cumsum(t == 1)[last] / n1)
  fp <- c(0, cumsum(t == 0)[last] / n0)
  tp[findInterval(grid, fp)]
}
roc_curves <- function(res, score_cols, grid = seq(0, 1, by = 0.02), keep_scenarios = NULL) {
  s <- res$subjects %>% filter(Eval)
  if (!is.null(keep_scenarios)) {
    ids <- res$design$Task_ID[res$design$Scenario_ID %in% keep_scenarios]
    s <- s %>% filter(Task_ID %in% ids)
  }
  per <- lapply(split(s, s$Task_ID), function(d) {
    bind_rows(lapply(names(score_cols), function(m) {
      sc <- d[[score_cols[[m]]]]
      if (all(is.na(sc))) return(NULL)
      data.frame(Method = m, FPR = grid, TPR = roc_at_grid(sc, d$Truth, grid))
    })) %>% mutate(Task_ID = d$Task_ID[1])
  })
  bind_rows(per) %>%
    left_join(res$design %>% select(Task_ID, Scenario_ID), by = "Task_ID") %>%
    group_by(Scenario_ID, Method, FPR) %>%
    summarise(TPR_mean = mean(TPR, na.rm = TRUE), TPR_mcse = sd(TPR, na.rm = TRUE) / sqrt(sum(!is.na(TPR))),
              n = sum(!is.na(TPR)), .groups = "drop") %>%
    left_join(distinct(res$design[, scenario_cols(res$design), drop = FALSE]), by = "Scenario_ID")
}

# ------------------------------------------------------------------------------
# Failures and convergence
# ------------------------------------------------------------------------------
failure_table <- function(res) {
  sc <- scenario_cols(res$design)
  planned <- res$design %>% count(Scenario_ID, name = "n_planned")
  res$fits %>%
    left_join(res$design %>% select(Task_ID, Scenario_ID), by = "Task_ID") %>%
    group_by(Scenario_ID, Method) %>%
    summarise(n_fits = n(), n_ok = sum(Status == "ok"),
              n_error = sum(Status == "error"), n_timeout = sum(Status == "timeout"),
              n_task_failed = sum(Status == "task_failed"),
              hit_maxit = mean(Hit_maxit, na.rm = TRUE),
              median_fit_sec = median(Fit_time, na.rm = TRUE),
              max_fit_sec = suppressWarnings(max(Fit_time, na.rm = TRUE)),
              getresponse_agree = mean(GetResponse_agree, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(planned, by = "Scenario_ID") %>%
    mutate(Fail_rate = 1 - n_ok / n_planned) %>%
    left_join(distinct(res$design[, sc, drop = FALSE]), by = "Scenario_ID")
}

# ------------------------------------------------------------------------------
# Precision report: was nsim large enough? 
# ------------------------------------------------------------------------------
precision_report <- function(perf, study, targets = c(TPR = 0.02, FDR = 0.01, AUC = 0.01)) {
  perf %>% filter(Alpha == ALPHAS[1]) %>%
    group_by(Method, Rule) %>%
    summarise(nsim_max = max(n_planned),
              max_TPR_mcse = max(TPR_mcse, na.rm = TRUE), max_FDR_mcse = max(FDR_mcse, na.rm = TRUE),
              max_AUC_mcse = max(AUC_mcse, na.rm = TRUE),
              nsim_for_TPR_target = nsim_needed(max(TPR_sd, na.rm = TRUE), targets[["TPR"]]),
              nsim_for_FDR_target = nsim_needed(max(FDP_sd, na.rm = TRUE), targets[["FDR"]]),
              nsim_for_AUC_target = nsim_needed(max(AUC_sd, na.rm = TRUE), targets[["AUC"]]),
              .groups = "drop") %>%
    mutate(Study = study, .before = 1)
}

# ------------------------------------------------------------------------------
# Plot style 
# ------------------------------------------------------------------------------
CELL_COLOURS <- c("High" = "deeppink", "Medium" = "steelblue3", "Low" = "orange",
                  "Very Low" = "purple", "Extremely Low" = "darkgreen")
METHOD_LINETYPES <- c("MIMOSA2" = "solid", "DiD" = "dotted", "DiD (BH)" = "dotted",
                      "DiD (unadjusted)" = "dashed",
                      "MIMOSA2_pooled" = "solid", "MIMOSA2_separate" = "dashed")
theme_mimosa <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey30"),
          plot.caption = element_text(size = 8, colour = "grey30", hjust = 0),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text = element_text(face = "bold", size = 9),
          strip.background = element_rect(fill = "grey92", colour = NA),
          panel.grid.minor = element_blank(),
          legend.position = "bottom", legend.box = "horizontal")
}
effect_lab <- function(x) formatC(x, format = "e", digits = 1)
res_lab <- function(x) paste0(round(100 * res_rho(x)), "% responders")
method_label <- function(method, rule) {
  case_when(method == "DiD" & rule == "BH" ~ "DiD (BH)",
            method == "DiD" & rule == "unadjusted" ~ "DiD (unadjusted)",
            TRUE ~ method)
}
save_fig <- function(p, name, width, height) {
  f <- file.path(FIG_DIR, paste0(name, ".pdf"))
  ggsave(f, plot = p, width = width, height = height)
  message("saved ", f)
}
save_tab <- function(x, name) {
  f <- file.path(TAB_DIR, paste0(name, ".csv"))
  write.csv(x, f, row.names = FALSE)
  message("saved ", f)
}
# "0.912 (0.004)" formatting for tables: estimate (MCSE)
fmt_mcse <- function(est, se, digits = 3) ifelse(is.na(est), "-", sprintf(paste0("%.", digits, "f (%.", digits, "f)"), est, se))
