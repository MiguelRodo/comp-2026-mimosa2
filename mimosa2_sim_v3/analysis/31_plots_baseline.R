# =============================================================================
# 31_plots_baseline.R : figures and tables for Study 1 (baseline)
# v3.2 (14 Oct 2026): thesis figures only. See analysis_functions.R for the style.
# =============================================================================
# Figures made (placement agreed 14 Oct):
#   Section 4.1 : base_plot_clean, ROC_plot_clean, base_dAUC_summary, fdr_plot_clean
#   Appendix    : base_plot, base_auc, base_dAUC
# No longer made (not used in the thesis): base_specificity_*, base_fdr_*,
#   base_dAUC_clean. Their numbers are still in _tables/baseline_performance.csv.
# All tables (CSV) are made exactly as before.
# =============================================================================
source("analysis/analysis_functions.R")
res  <- load_results("baseline")
perf <- performance_table(res) %>%
  mutate(Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")))
add_strips <- function(d) d %>% mutate(Res_lab = math_factor(Res_prop, rho_math),
                                       P_lab   = math_factor(P, P_math))
main <- perf %>% filter(Res_prop != "Prop_0.00",
                        (Method == "MIMOSA2" & Rule == "BFDR") | (Method == "DiD" & Rule == "BH")) %>%
  mutate(Label = plot_method(Method, Rule)) %>% add_strips()
effects <- sort(unique(main$Effect))
P_SHOW  <- if (50 %in% main$P) 50 else max(main$P)        # single P for ROC / calibration figures
CLEAN_P   <- c(10, 50, 100)
CLEAN_RHO <- c("Prop_0.10", "Prop_0.50", "Prop_0.90")
FG <- facet_grid(Res_lab ~ P_lab, labeller = label_parsed)

band_plot <- function(d, y, se, ylab, ylim = c(0, 1)) {
  ggplot(d, aes(x = Effect, y = .data[[y]], colour = Cell_range, fill = Cell_range,
                linetype = Label, group = interaction(Cell_range, Label))) +
    geom_ribbon(aes(ymin = pmax(.data[[y]] - 1.96 * .data[[se]], ylim[1]),
                    ymax = pmin(.data[[y]] + 1.96 * .data[[se]], ylim[2])), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
    scale_x_log10(breaks = effects, labels = effect_axis) +
    scale_cell() + scale_method() +
    coord_cartesian(ylim = ylim) +
    labs(x = X_EFFECT, y = ylab) +
    theme_mimosa() + tilt_x()
}

# ---- 1. Sensitivity at 1% nominal FDR ---------------------------------------------
d <- main %>% filter(Alpha == 0.01)
save_fig(band_plot(d, "TPR", "TPR_mcse", "True positive rate") + FG, "base_plot", height = 10.5)    # appendix
d2 <- d %>% filter(P %in% CLEAN_P, Res_prop %in% CLEAN_RHO)
save_fig(band_plot(d2, "TPR", "TPR_mcse", "True positive rate") + FG, "base_plot_clean", height = 8)  # 4.1

# ---- 2. Mean per-dataset AUC (appendix) ----------------------------------------------
# The AUC does not depend on the calling rule; the unadjusted DiD rows carry it.
dA <- perf %>% filter(Res_prop != "Prop_0.00", Alpha == 0.01,
                      (Method == "MIMOSA2" & Rule == "BFDR") | (Method == "DiD" & Rule == "unadjusted")) %>%
  mutate(Label = plot_method(Method, Rule, auc = TRUE)) %>% add_strips()
save_fig(band_plot(dA, "AUC", "AUC_mcse", "Average AUC", ylim = c(0.5, 1)) + FG, "base_auc", height = 10.5)

# ---- 3. Paired Delta AUC, every scenario (appendix) ------------------------------------
dd <- paired_auc_diff(res, "MIMOSA2", "DiD") %>% filter(Res_prop != "Prop_0.00") %>%
  mutate(Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low"))) %>% add_strips()
save_tab(dd, "baseline_dAUC_table")
p <- ggplot(dd, aes(Effect, dAUC_mean, colour = Cell_range)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey20") +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.15),
                  size = 0.2, linewidth = 0.4) +
  scale_x_log10(breaks = effects, labels = effect_axis) +
  scale_cell(aesthetics = "colour") + FG +
  labs(x = X_EFFECT, y = DAUC_LAB(paste("MIMOSA2 -", LAB_DID_AUC))) +
  theme_mimosa() + tilt_x()
save_fig(p, "base_dAUC", height = 10.5)

# ---- 4. One-panel summary of the paired Delta AUC (4.1) ------------------------------
# Median over the 30 scenarios (6 values of P x 5 values of rho) for each effect
# size and cell count; bars are the interquartile range (NOT MCSE: say so in the caption).
dd_sum <- dd %>%
  group_by(Effect, Cell_range) %>%
  summarise(dAUC_med = median(dAUC_mean), q25 = quantile(dAUC_mean, 0.25),
            q75 = quantile(dAUC_mean, 0.75), n_scen = n(), .groups = "drop")
save_tab(dd_sum, "baseline_dAUC_summary")
p <- ggplot(dd_sum, aes(Effect, dAUC_med, colour = Cell_range, group = Cell_range)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey20") +
  geom_line(position = position_dodge(width = 0.12), linewidth = 0.7) +
  geom_pointrange(aes(ymin = q25, ymax = q75), position = position_dodge(width = 0.12),
                  size = 0.35, linewidth = 0.6) +
  scale_x_log10(breaks = effects, labels = effect_axis) +      # one panel: labels fit untilted
  scale_cell(aesthetics = "colour") +
  labs(x = X_EFFECT, y = DAUC_LAB(paste("MIMOSA2 -", LAB_DID_AUC))) +
  theme_mimosa()
save_fig(p, "base_dAUC_summary", height = 5)

# ---- 5. Vertically averaged ROC curves at P = P_SHOW (4.1) ------------------------------
sel_eff <- effects[unique(c(1, ceiling(length(effects) / 2), length(effects)))]
keep <- res$design %>% filter(P == P_SHOW, Effect %in% sel_eff, Res_prop %in% CLEAN_RHO) %>%
  pull(Scenario_ID) %>% unique()
roc <- roc_curves(res, c(MIMOSA2 = "MIMOSA2_prob", DiD = "DiD_z"), keep_scenarios = keep) %>%
  mutate(Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")),
         Label   = plot_method(Method, "unadjusted", auc = TRUE),
         Res_lab = math_factor(Res_prop, rho_math),
         Eff_lab = math_factor(Effect, delta_math))
p <- ggplot(roc, aes(FPR, TPR_mean, colour = Cell_range, linetype = Label,
                     group = interaction(Label, Cell_range))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey10", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  scale_cell(aesthetics = "colour") + scale_method() +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  facet_grid(Res_lab ~ Eff_lab, labeller = label_parsed) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_mimosa()
save_fig(p, "ROC_plot_clean", height = 8.6)

# ---- 6. FDR calibration at P = P_SHOW (4.1) ---------------------------------------------
fc <- fdr_curve(res) %>% filter(P == P_SHOW, Res_prop %in% CLEAN_RHO, Effect %in% sel_eff) %>%
  mutate(Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")),
         Label   = method_factor(ifelse(Method == "DiD (BH)", LAB_DID_BH, Method)),
         Res_lab = math_factor(Res_prop, rho_math),
         Eff_lab = math_factor(Effect, delta_math))
save_tab(fc, "baseline_fdr_calibration")
p <- ggplot(fc, aes(Alpha, FDR, colour = Cell_range, fill = Cell_range, linetype = Label,
                    group = interaction(Label, Cell_range))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey10", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmax(FDR - 1.96 * FDR_mcse, 0), ymax = FDR + 1.96 * FDR_mcse),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) +
  scale_cell() + scale_method() +
  # Zoomed to the levels discussed in the text; the full 0.01-0.20 range is in
  # _tables/baseline_fdr_calibration.csv. To show it all: xlim = c(0.01, 0.20), ylim = c(0, 0.35).
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10)) +
  coord_cartesian(xlim = c(0.01, 0.10), ylim = c(0, 0.25)) +
  facet_grid(Res_lab ~ Eff_lab, labeller = label_parsed) +
  labs(x = "Nominal FDR", y = "Observed FDR") +
  theme_mimosa()
save_fig(p, "fdr_plot_clean", height = 8.6)

# ---- 7. Null scenarios (rho = 0) ------------------------------------------------------
null <- perf %>% filter(Res_prop == "Prop_0.00") %>%
  transmute(P, Cell_range, Method = method_label(Method, Rule), Alpha,
            `P(any false call)` = fmt_mcse(Any_FP, Any_FP_mcse),
            `Per-subject FPR` = fmt_mcse(FPR, FPR_mcse), `Mean calls` = round(Calls, 2), Fail_rate) %>%
  arrange(Alpha, Method, Cell_range, P)
save_tab(null, "baseline_null_scenarios")
print(knitr::kable(null))

# ---- 8. Main summary table at P = P_SHOW (estimate (MCSE)) -----------------------------
tab <- perf %>% filter(P == P_SHOW, Res_prop != "Prop_0.00", Alpha == 0.01,
                       (Method == "MIMOSA2" & Rule == "BFDR") | (Method == "DiD" & Rule == "BH")) %>%
  arrange(Res_prop, desc(Effect), Cell_range) %>%
  transmute(Res_prop, Effect = effect_lab(Effect), Cell_range, Method = method_label(Method, Rule),
            TPR = fmt_mcse(TPR, TPR_mcse), FDR = fmt_mcse(FDR, FDR_mcse), AUC = fmt_mcse(AUC, AUC_mcse),
            Fail = sprintf("%.1f%%", 100 * Fail_rate)) %>%
  pivot_wider(names_from = Method, values_from = c(TPR, FDR, AUC, Fail))
save_tab(tab, sprintf("baseline_summary_P%d", P_SHOW))
print(knitr::kable(tab, caption = sprintf("Baseline, P = %d, alpha = 1%%: estimate (Monte Carlo SE)", P_SHOW)))
