# =============================================================================
# analysis_functions.R : performance measures WITH Monte Carlo SEs
# (Morris et al. 2019, Sections 5.1-5.4 and 6.2).  NEW FILE.
# =============================================================================
# Used by analysis/30-35. Run the analysis scripts from the mimosa2_sim_v3
# folder, on the cluster or locally after copying _simulations/ across:
#   SIM_PROFILE=standard Rscript analysis/30_performance_tables.R
# (In RStudio: setwd() to the mimosa2_sim_v3 folder and
#  Sys.setenv(SIM_PROFILE = "standard") before sourcing.)
#
# Key definitions (state these in Chapter 3, "Performance measures"):
#  * All measures are computed PER DATASET and then averaged over the nsim
#    datasets of a scenario; the Monte Carlo SE of a mean is SD / sqrt(n)
#    (Morris Table 6). Scenarios are never pooled (Morris 5.2).
#  * TPR (sensitivity) = TP / (TP + FN) in datasets with >= 1 responder.
#  * TNR (specificity) = TN / (TN + FP); FPR = 1 - TNR.
#  * FDP = FP / (TP + FP), defined as 0 when nothing is called; the FDR is the
#    mean FDP over datasets (this is the quantity the Bayesian FDR rule and
#    Benjamini-Hochberg aim to control).
#  * P(any FP) = proportion of datasets with at least one false positive
#    (equals the FDR in the null scenarios, where every call is false).
#  * AUC = per-dataset Mann-Whitney AUC (probability that a random responder
#    in the dataset outscores a random non-responder). [CHANGE] Version 2
#    pooled subjects ACROSS datasets before computing one AUC, which mixes
#    posterior probabilities from different fitted models (audit A5).
#  * Delta AUC = paired difference AUC_MIMOSA2 - AUC_DiD in the same dataset;
#    its MCSE is SD(differences) / sqrt(n) (Morris 4.3: same datasets for both
#    methods). [CHANGE] Replaces the Hanley-McNeil intervals, which treated all
#    pooled subjects as independent, fixed 50% positives and assumed r = 0.5.
#  * Fit failures (error / timeout / task failed) are the first performance
#    measure (Morris 5.1) and are reported per scenario. Other measures are
#    computed on successful fits only.
# =============================================================================
# [CHANGE 4 Oct] config.R first, so .libPaths() includes your R_libs folder
# BEFORE any package is loaded (previously dplyr etc. were loaded first).
for (f in c("R/config.R", "R/scenarios.R")) source(f)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(ggplot2)
})

# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
# Profile used for each study in the write-up (12 Oct 2026)
STUDY_PROFILE <- c(baseline = "standard", maxit_check = "standard",
                   prior = "extended", heterogeneity = "extended", imbalance = "extended")
load_results <- function(study, profile = NULL) {
  if (is.null(profile)) profile <- if (study %in% names(STUDY_PROFILE)) STUDY_PROFILE[[study]] else PROFILE
  f <- file.path(OUT_DIR, study, sprintf("%s_results_%s.rds", study, profile))
  if (!file.exists(f)) stop("No results for ", study, " (", f, "). Run the simulation / sims/29_combine.R first.")
  message("loading ", f)
  readRDS(f)
}
scenario_cols <- function(design) setdiff(names(design), c("Rep", "Task_ID", "Seed", "In_standard", "In_smoke"))

# ---------------------------------------------------------------------------
# Monte Carlo SE helpers (Morris Table 6)
# ---------------------------------------------------------------------------
mc_n    <- function(x) sum(!is.na(x))
mc_mean <- function(x) if (mc_n(x) == 0) NA_real_ else mean(x, na.rm = TRUE)
mc_se   <- function(x) if (mc_n(x) < 2) NA_real_ else sd(x, na.rm = TRUE) / sqrt(mc_n(x))
nsim_needed <- function(sd, target) ceiling((sd / target)^2)       # Morris 5.3, eq. (1) generalised

# ---------------------------------------------------------------------------
# Per-dataset measures
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Performance table: one row per scenario x method x rule x alpha
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# FDR calibration curve (per dataset, then averaged)  [CHANGE, audit A2]
# Version 2 thresholded each subject's LOCAL posterior (1 - P(R) <= a) and
# pooled subjects across datasets and P. The rule actually used is the
# cumulative Bayesian FDR (q-value) of each dataset, so here a subject is
# called at level a if q < a (exactly as MIMOSA2::getResponse), the FDP is computed per dataset, and FDR(a) is
# the mean FDP over datasets with its MCSE.
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Vertically averaged ROC curves (Fawcett 2006): TPR of each dataset at a grid
# of FPR values, averaged over datasets.  [CHANGE] Version 2 drew ONE ROC
# curve per scenario from subjects pooled across datasets (and across P).
# ---------------------------------------------------------------------------
roc_at_grid <- function(score, truth, grid) {
  ok <- !is.na(score); score <- score[ok]; truth <- truth[ok]
  n1 <- sum(truth == 1); n0 <- sum(truth == 0)
  if (n1 == 0 || n0 == 0) return(rep(NA_real_, length(grid)))
  o <- order(score, decreasing = TRUE); s <- score[o]; t <- truth[o]
  last <- c(s[-1] != s[-length(s)], TRUE)           # end of each block of tied scores
  tp <- c(0, cumsum(t == 1)[last] / n1)
  fp <- c(0, cumsum(t == 0)[last] / n0)
  # [CHANGE v3.1, 6 Oct] A block of TIED scores (e.g. the ~10% of subjects whose
  # MIMOSA2 posterior is exactly 0, or DiD z = 0) is a straight diagonal segment
  # of the ROC curve (Fawcett 2006, the same convention as the Mann-Whitney AUC
  # with ties counted 1/2). v3.0 read the curve as a step function
  # (findInterval), which drew such a block as a long flat line followed by a
  # jump at FPR = 1, and made those ROC curves look worse than their AUC.
  # Linear interpolation; at vertical segments the top value is used.
  stats::approx(fp, tp, xout = grid, method = "linear", ties = max, rule = 2)$y
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

# ---------------------------------------------------------------------------
# [NEW v3.1] Detection by response profile. The eight DGM profiles behave very
# differently: e.g. R3 (p(s,1) = p(s,0)) has a true Delta that does not grow
# with delta, so it caps the TPR of every method; NSR / NR2 drive the false
# calls. Call rate = share of subjects of that profile called at level alpha,
# pooled over the datasets of the scenario (descriptive; no MCSE).
# ---------------------------------------------------------------------------
profile_table <- function(res, alpha = 0.01) {
  s <- res$subjects %>% filter(Eval)
  qcols <- grep("^MIMOSA2.*_q$", names(s), value = TRUE)
  s <- s %>% group_by(Task_ID) %>% mutate(DiD_BH = p.adjust(DiD_p, "BH")) %>% ungroup()
  calls <- lapply(qcols, function(qc) data.frame(Task_ID = s$Task_ID, Profile = s$Profile, True_delta = s$True_delta,
                                                 Method = sub("_q$", "", qc), Call = s[[qc]] < alpha))
  calls[[length(calls) + 1]] <- data.frame(Task_ID = s$Task_ID, Profile = s$Profile, True_delta = s$True_delta,
                                           Method = "DiD (BH)", Call = s$DiD_BH <= alpha)
  bind_rows(calls) %>%
    left_join(res$design %>% select(Task_ID, Scenario_ID), by = "Task_ID") %>%
    group_by(Scenario_ID, Method, Profile) %>%
    summarise(n_subjects = n(), call_rate = mean(Call, na.rm = TRUE),
              median_true_delta = median(True_delta), .groups = "drop") %>%
    mutate(Alpha = alpha) %>%
    left_join(distinct(res$design[, scenario_cols(res$design), drop = FALSE]), by = "Scenario_ID")
}

# ---------------------------------------------------------------------------
# Failures and convergence (Morris 5.1)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Precision report: was nsim large enough? (Morris 5.3, 6.2)
# ---------------------------------------------------------------------------
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

# =============================================================================
# Thesis plot style  [v3.2, 14 Oct 2026; replaces the v3.1 style section]
# Decisions are logged in the project doc claude/figure_conventions.md.
#  * Every figure is drawn FIG_W = 8 in wide and included at \textwidth
#    (about 6.3 in), i.e. scaled by about 0.79. The 14 pt axis titles and
#    10 pt legend text then print at about 11 pt and 8 pt.
#  * Numbers on effect-size axes and in strips: 2 decimals x 10^k
#    (6.25 x 10^-5, 1.00 x 10^-3), never "6.2e-05".
#  * Strips: rho = 0.1, P = 50, delta = ... (plotmath, label_parsed).
#  * Legends list MIMOSA2 first. DiD calls (TPR, FDR) are "DiD-BH"; for AUC,
#    ROC and Delta AUC the label is LAB_DID_AUC (BH does not change the ranking,
#    so the AUC is the same with or without it).
#  * PDFs are written with cairo_pdf, which embeds the fonts. The plain pdf()
#    device does not, which is why the Greek Delta disappeared from the v3.1
#    figures. If cairo is unavailable on the cluster, run the 3x scripts locally.
# =============================================================================
FIG_W <- 8

# ---- Colours --------------------------------------------------------------------
PAL_GREEN  <- "#006400"   # accent 1 (two-way comparisons)
PAL_ORCHID <- "#DA70D6"   # accent 2 (two-way comparisons)
PAL_REF    <- "black"     # reference scenario / reference family
PAL_DID    <- "#4D4D4D"   # DiD comparator
PAL_OCHRE  <- "#B07D2B"   # extra muted tone (only when 5 categories are needed)
PAL_SLATE  <- "#6B7A8F"   # extra muted tone
CELL_COLOURS <- c("High" = "#C2477F", "Medium" = "#8A52B8", "Low" = "#3F7CC0",
                  "Very Low" = "#2E9488", "Extremely Low" = "#4E9A45")
CELL_LABELS  <- c("High" = "High", "Medium" = "Medium", "Low" = "Low",
                  "Very Low" = "Very low", "Extremely Low" = "Extremely low")
CELL_TITLE   <- "Parent cell count"
scale_cell <- function(name = CELL_TITLE, aesthetics = c("colour", "fill"))
  scale_colour_manual(values = CELL_COLOURS, labels = CELL_LABELS, name = name, aesthetics = aesthetics)

# ---- Methods ----------------------------------------------------------------------
LAB_DID_BH  <- "DiD-BH"
LAB_DID_AUC <- "DiD"          # set to "DiD-BH" to use that label on the AUC figures too
METHOD_ORDER     <- c("MIMOSA2", LAB_DID_BH, LAB_DID_AUC, "DiD (unadjusted)")
METHOD_LINETYPES <- setNames(c("solid", "dotted", "dotted", "dashed"), METHOD_ORDER)
METHOD_TITLE     <- "Model"
# Factor with MIMOSA2 first: the legend follows the factor levels
method_factor <- function(x) factor(x, levels = intersect(unique(METHOD_ORDER), unique(x)))
plot_method <- function(method, rule, auc = FALSE) {
  lab <- case_when(method == "DiD" & auc ~ LAB_DID_AUC,
                   method == "DiD" & rule == "BH" ~ LAB_DID_BH,
                   method == "DiD" & rule == "unadjusted" ~ "DiD (unadjusted)",
                   TRUE ~ method)
  method_factor(lab)
}
scale_method <- function(name = METHOD_TITLE) scale_linetype_manual(values = METHOD_LINETYPES, name = name)

# ---- Theme --------------------------------------------------------------------------
theme_mimosa <- function() {
  theme_bw(base_size = 11) +
    theme(axis.title        = element_text(size = 14),
          axis.title.x      = element_text(margin = margin(t = 8)),
          axis.title.y      = element_text(margin = margin(r = 8)),
          axis.text         = element_text(size = 10),
          strip.text        = element_text(face = "bold", size = 11),
          strip.background  = element_blank(),
          panel.border      = element_rect(fill = NA, colour = "grey30", linewidth = 0.3),
          panel.spacing     = unit(0.6, "lines"),
          panel.grid.minor  = element_blank(),
          legend.position   = "bottom", legend.box = "horizontal",
          legend.title      = element_text(size = 10),
          legend.text       = element_text(size = 10),
          legend.key.width  = unit(2, "lines"),        # long enough to tell dotted from solid
          legend.background = element_blank(),
          legend.key        = element_blank())
}
# 45-degree x labels: only for effect-size axes in multi-panel figures, where
# 2-decimal scientific labels do not fit horizontally
# (extra horizontal panel spacing so the tilted labels of neighbouring panels do not touch)
tilt_x <- function() theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                           panel.spacing.x = unit(1.4, "lines"))

save_fig <- function(p, name, height, width = FIG_W) {
  # No title, subtitle or caption inside the figure (captions go in LaTeX)
  p <- p + theme(plot.title = element_blank(), plot.subtitle = element_blank(),
                 plot.caption = element_blank())
  f <- file.path(FIG_DIR, paste0(name, ".pdf"))
  dev <- if (isTRUE(capabilities("cairo"))) cairo_pdf else grDevices::pdf
  ggsave(f, plot = p, width = width, height = height, device = dev)
  message("saved ", f)
}

# ---- Labels -------------------------------------------------------------------------
# Plotmath strings: sci_math(6.25e-5) = '"6.25" %*% 10^-5'; 0 stays "0"
sci_math <- function(x, digits = 2) {
  out <- rep("''", length(x))
  ok <- !is.na(x)
  out[ok & x == 0] <- "0"
  nz <- which(ok & x != 0)
  if (length(nz)) {
    e <- floor(log10(abs(x[nz])) + 1e-9)
    m <- x[nz] / 10^e
    up <- abs(round(m, digits)) >= 10            # 9.999 -> 10.00: move to the next power
    m[up] <- m[up] / 10; e[up] <- e[up] + 1
    # mantissa quoted, so plotmath prints it as text and keeps trailing zeros (2.50, 1.00)
    out[nz] <- sprintf(paste0('"%.', digits, 'f" %%*%% 10^%d'), m, e)
  }
  out
}
sci_axis    <- function(digits = 2) function(b) parse(text = sci_math(b, digits))
effect_axis <- sci_axis(2)                                    # use as scale_x_*(labels = effect_axis)
math_factor <- function(x, fmt) { u <- sort(unique(x)); factor(fmt(x), levels = unique(fmt(u))) }
rho_math    <- function(x) paste0("rho == ", res_rho(x))      # x = "Prop_0.10" -> rho = 0.1
P_math      <- function(x) paste0("italic(P) == ", x)
delta_math  <- function(x) paste0("delta == ", sci_math(x))
X_EFFECT    <- expression("Effect size " * delta)
DAUC_LAB    <- function(what) bquote(Delta * "AUC (" * .(what) * ")")

# ---- Table helpers (unchanged: the table code reads these column names) ------------
effect_lab <- function(x) formatC(x, format = "e", digits = 1)
res_lab <- function(x) paste0(round(100 * res_rho(x)), "% responders")
method_label <- function(method, rule) {
  case_when(method == "DiD" & rule == "BH" ~ "DiD (BH)",
            method == "DiD" & rule == "unadjusted" ~ "DiD (unadjusted)",
            TRUE ~ method)
}
save_tab <- function(x, name) {
  f <- file.path(TAB_DIR, paste0(name, ".csv"))
  write.csv(x, f, row.names = FALSE)
  message("saved ", f)
}
# "0.912 (0.004)" formatting for tables: estimate (MCSE)
fmt_mcse <- function(est, se, digits = 3) ifelse(is.na(est), "-", sprintf(paste0("%.", digits, "f (%.", digits, "f)"), est, se))
