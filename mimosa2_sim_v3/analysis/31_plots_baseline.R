# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# 31_plots_baseline.R: Figures and Tables for Study 1
# Replaces Baseline_simulation_plots.R and Simulation_2.0.R
# =============================================================================
# What changed relative to Baseline_simulation_plots.R:
# - Bands are mean +/- 1.96 x Monte Carlo SE of the per-dataset measure
#   (Morris 5.2, 6.2). Version 2 used bootstrap CIs (mean_cl_boot), which
#   estimate the same thing; MCSE is now computed explicitly and stored.
# - The DiD comparator is drawn next to MIMOSA2 in every figure (Morris 6.2:
#   put the methods side by side). DiD calls use Benjamini-Hochberg so that
#   both methods target the FDR at the same nominal level.
# - NEW figures: specificity / observed FDR (audit A3), per-dataset AUC and
#   paired Delta AUC with MCSE (A4, A5), FDR calibration using the correct
#   cumulative q-value rule (A2), and the null scenarios (B4).
# - ROC curves are vertically averaged per-dataset curves for ONE value of P
#   at a time (version 2 pooled subjects across datasets and across P; A6).
# - Fixed: x-axis label "1-Sensitivity" (should be 1 - specificity), y label
#   typo "Observed FDF", AUROC_matrix used an object (AUROC) whose code was
#   commented out, and the first AUROC table had its two method columns
#   swapped.
# - The Hanley-McNeil Delta-AUC intervals are removed (not valid here; see
#   analysis_functions.R).
# ==============================================================================
source("analysis/analysis_functions.R")
res  <- load_results("baseline")
perf <- performance_table(res) %>%
  mutate(Label = method_label(Method, Rule),
         Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")))
main <- perf %>% filter(Res_prop != "Prop_0.00",
                        (Method == "MIMOSA2" & Rule == "BFDR") | (Method == "DiD" & Rule == "BH")) %>%
  mutate(Res_lab = factor(res_lab(Res_prop), levels = res_lab(sort(unique(Res_prop)))),
         P_lab = factor(paste0("P = ", P), levels = paste0("P = ", sort(unique(P)))))
effects <- sort(unique(main$Effect))
P_SHOW <- if (50 %in% main$P) 50 else max(main$P)          # single P for ROC / calibration figures

band_plot <- function(d, y, se, ylab, title, subtitle, hline = NULL, ylim = c(0, 1)) {
  p <- ggplot(d, aes(x = Effect, y = .data[[y]], colour = Cell_range, fill = Cell_range,
                     linetype = Label, group = interaction(Cell_range, Label))) +
    geom_ribbon(aes(ymin = pmax(.data[[y]] - 1.96 * .data[[se]], ylim[1]),
                    ymax = pmin(.data[[y]] + 1.96 * .data[[se]], ylim[2])), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
    scale_x_log10(breaks = effects, labels = effect_lab(effects)) +
    scale_colour_manual(values = CELL_COLOURS, name = "Parent cell count") +
    scale_fill_manual(values = CELL_COLOURS, name = "Parent cell count") +
    scale_linetype_manual(values = c("MIMOSA2" = "solid", "DiD (BH)" = "dotted"), name = "Model") +
    coord_cartesian(ylim = ylim) +
    labs(title = title, subtitle = subtitle, x = "Effect size", y = ylab,
         caption = sprintf("Lines: mean over datasets; bands: +/- 1.96 Monte Carlo SE (nsim = %d per scenario).",
                           max(d$n_planned))) +
    theme_mimosa()
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", colour = "grey20")
  p
}

# ----------- 1. Sensitivity at 1% nominal FDR (all scenarios) -----------------
d <- main %>% filter(Alpha == 0.01)
p <- band_plot(d, "TPR", "TPR_mcse", "True positive rate", "Sensitivity analysis of MIMOSA2.",
               "TPR at 1% nominal FDR (MIMOSA2: Bayesian FDR; DiD: Benjamini-Hochberg).") +
  facet_grid(Res_lab ~ P_lab)
save_fig(p, "base_plot", 11, 11)

# ---------- 2. Clean version: P = 10, 50, 100; 10/50/90% responders -----------
d2 <- d %>% filter(P %in% c(10, 50, 100), Res_prop %in% c("Prop_0.10", "Prop_0.50", "Prop_0.90"))
if (nrow(d2)) save_fig(band_plot(d2, "TPR", "TPR_mcse", "True positive rate", "Sensitivity analysis of MIMOSA2.",
                                 "TPR at 1% nominal FDR.") + facet_grid(Res_lab ~ P_lab), "base_plot_clean", 8, 7)

# ------------- 3. Specificity and observed FDR (NEW, audit A3) ----------------
for (a in ALPHAS) {
  da <- main %>% filter(Alpha == a)
  save_fig(band_plot(da, "TNR", "TNR_mcse", "Specificity (TNR)", "Specificity of responder calls.",
                     sprintf("Nominal FDR %g%%.", 100 * a), ylim = c(0.5, 1)) + facet_grid(Res_lab ~ P_lab),
           sprintf("base_specificity_%03d", round(1000 * a)), 11, 11)
  save_fig(band_plot(da, "FDR", "FDR_mcse", "Observed FDR (mean FDP)", "False discovery rate of responder calls.",
                     sprintf("Nominal level %g%% (dashed). Above the line = FDR not controlled.", 100 * a),
                     hline = a, ylim = c(0, max(0.3, 2.5 * a))) + facet_grid(Res_lab ~ P_lab),
           sprintf("base_fdr_%03d", round(1000 * a)), 11, 11)
}

# ---------------- 4. Per-dataset AUC and paired Delta AUC ---------------------
dA <- perf %>% filter(Res_prop != "Prop_0.00", Alpha == 0.01,
                      (Method == "MIMOSA2" & Rule == "BFDR") | (Method == "DiD" & Rule == "unadjusted")) %>%
  mutate(Label = ifelse(Method == "DiD", "DiD (BH)", "MIMOSA2"),   # same AUC for any DiD rule
         Res_lab = factor(res_lab(Res_prop), levels = res_lab(sort(unique(Res_prop)))),
         P_lab = factor(paste0("P = ", P), levels = paste0("P = ", sort(unique(P)))))
save_fig(band_plot(dA, "AUC", "AUC_mcse", "Mean per-dataset AUC", "Discrimination of MIMOSA2 and DiD.",
                   "AUC computed within each dataset, then averaged.", ylim = c(0.5, 1)) +
           scale_linetype_manual(values = c("MIMOSA2" = "solid", "DiD (BH)" = "dotted"),
                                 labels = c("MIMOSA2" = "MIMOSA2", "DiD (BH)" = "DiD"), name = "Model") +
           facet_grid(Res_lab ~ P_lab), "base_auc", 11, 11)

dd <- paired_auc_diff(res, "MIMOSA2", "DiD") %>% filter(Res_prop != "Prop_0.00") %>%
  mutate(Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")),
         Res_lab = factor(res_lab(Res_prop), levels = res_lab(sort(unique(Res_prop)))),
         P_lab = factor(paste0("P = ", P), levels = paste0("P = ", sort(unique(P)))))
save_tab(dd, "baseline_dAUC_table")
p <- ggplot(dd, aes(Effect, dAUC_mean, colour = Cell_range)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.15), size = 0.25) +
  scale_x_log10(breaks = effects, labels = effect_lab(effects)) +
  scale_colour_manual(values = CELL_COLOURS, name = "Parent cell count") +
  facet_grid(Res_lab ~ P_lab) +
  labs(title = "Paired difference in AUC (MIMOSA2 - DiD).", subtitle = "Same datasets for both methods.",
       x = "Effect size", y = expression(Delta * " AUC"),
       caption = "Points: mean paired difference; bars: +/- 1.96 Monte Carlo SE. Above 0 = MIMOSA2 better.") +
  theme_mimosa()
save_fig(p, "base_dAUC", 11, 11)

# -------------- 5. Vertically averaged ROC curves, P = P_SHOW -----------------
sel_eff <- effects[unique(c(1, ceiling(length(effects) / 2), length(effects)))]
keep <- res$design %>% filter(P == P_SHOW, Effect %in% sel_eff,
                              Res_prop %in% c("Prop_0.10", "Prop_0.50", "Prop_0.90")) %>%
  pull(Scenario_ID) %>% unique()
roc <- roc_curves(res, c(MIMOSA2 = "MIMOSA2_prob", DiD = "DiD_z"), keep_scenarios = keep) %>%
  mutate(Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")),
         Res_lab = factor(res_lab(Res_prop), levels = res_lab(c("Prop_0.10", "Prop_0.50", "Prop_0.90"))),
         Eff_lab = factor(paste0("Effect: ", effect_lab(Effect)), levels = paste0("Effect: ", effect_lab(sel_eff))))
p <- ggplot(roc, aes(FPR, TPR_mean, colour = Cell_range, linetype = Method, group = interaction(Method, Cell_range))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey10", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = CELL_COLOURS, name = "Parent cell count") +
  scale_linetype_manual(values = c(MIMOSA2 = "solid", DiD = "dotted"), name = "Model") +
  facet_grid(Res_lab ~ Eff_lab) + coord_equal() +
  labs(title = "Simulation performance of MIMOSA2.", subtitle = sprintf("Vertically averaged ROC curves, P = %d.", P_SHOW),
       x = "1 - Specificity", y = "Sensitivity",
       caption = "Each curve averages the per-dataset ROC curves (TPR at fixed FPR) over the nsim datasets.") +
  theme_mimosa()
save_fig(p, "ROC_plot_clean", 9, 9)

# --------------------- 6. FDR calibration (correct rule) ----------------------
fc <- fdr_curve(res) %>% filter(P == P_SHOW, Res_prop %in% c("Prop_0.10", "Prop_0.50", "Prop_0.90"), Effect %in% sel_eff) %>%
  mutate(Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")),
         Res_lab = factor(res_lab(Res_prop), levels = res_lab(c("Prop_0.10", "Prop_0.50", "Prop_0.90"))),
         Eff_lab = factor(paste0("Effect: ", effect_lab(Effect)), levels = paste0("Effect: ", effect_lab(sel_eff))))
save_tab(fc, "baseline_fdr_calibration")
p <- ggplot(fc, aes(Alpha, FDR, colour = Cell_range, fill = Cell_range, linetype = Method, group = interaction(Method, Cell_range))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey10", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmax(FDR - 1.96 * FDR_mcse, 0), ymax = FDR + 1.96 * FDR_mcse), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = CELL_COLOURS, name = "Parent cell count") +
  scale_fill_manual(values = CELL_COLOURS, name = "Parent cell count") +
  scale_linetype_manual(values = c("MIMOSA2" = "solid", "DiD (BH)" = "dotted"), name = "Model") +
  coord_cartesian(ylim = c(0, 0.3)) + facet_grid(Res_lab ~ Eff_lab) +
  labs(title = "False discovery rate.", subtitle = sprintf("Observed FDR vs nominal level, P = %d. On or below the dashed line = controlled.", P_SHOW),
       x = "Nominal FDR", y = "Observed FDR (mean FDP)",
       caption = "Per dataset: call subjects with q-value <= nominal level; FDP = FP / calls (0 if none). Bands: +/- 1.96 MCSE.") +
  theme_mimosa()
save_fig(p, "fdr_plot_clean", 9, 9)

# -------------- 7. Null scenarios (rho = 0): known properties -----------------
null <- perf %>% filter(Res_prop == "Prop_0.00") %>%
  transmute(P, Cell_range, Method = method_label(Method, Rule), Alpha,
            `P(any false call)` = fmt_mcse(Any_FP, Any_FP_mcse),
            `Per-subject FPR` = fmt_mcse(FPR, FPR_mcse), `Mean calls` = round(Calls, 2), Fail_rate) %>%
  arrange(Alpha, Method, Cell_range, P)
save_tab(null, "baseline_null_scenarios")
message("\nNull scenarios (no responders). FDR-controlling rules (MIMOSA2 BFDR, DiD BH) should give")
message("P(any false call) <= alpha; the unadjusted DiD test should give a per-subject FPR <= alpha.")
print(knitr::kable(null))

# ---------- 8. Main summary table at P = P_SHOW (estimate (MCSE)) -------------
tab <- perf %>% filter(P == P_SHOW, Res_prop != "Prop_0.00", Alpha == 0.01,
                       (Method == "MIMOSA2" & Rule == "BFDR") | (Method == "DiD" & Rule == "BH")) %>%
  arrange(Res_prop, desc(Effect), Cell_range) %>%
  transmute(Res_prop, Effect = effect_lab(Effect), Cell_range, Method = method_label(Method, Rule),
            TPR = fmt_mcse(TPR, TPR_mcse), FDR = fmt_mcse(FDR, FDR_mcse), AUC = fmt_mcse(AUC, AUC_mcse),
            Fail = sprintf("%.1f%%", 100 * Fail_rate)) %>%
  pivot_wider(names_from = Method, values_from = c(TPR, FDR, AUC, Fail))
save_tab(tab, sprintf("baseline_summary_P%d", P_SHOW))
print(knitr::kable(tab, caption = sprintf("Baseline, P = %d, alpha = 1%%: estimate (Monte Carlo SE)", P_SHOW)))
