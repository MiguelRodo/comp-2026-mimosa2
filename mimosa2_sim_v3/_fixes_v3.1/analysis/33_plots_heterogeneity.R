# =============================================================================
# 33_plots_heterogeneity.R : figures and tables for Study 3 (effect heterogeneity)
# Replaces Effect_heterogeneity_plot.R (version 2).
# =============================================================================
# [CHANGE] Performance is measured on the FOCAL (weak-effect) subgroup, which
# now really has the smaller effect (bug fix, scenarios.R). Measures are
# per-dataset with MCSE; the key comparison is the PAIRED difference
# pooled fit - separate fit on the same datasets. AUC is computed within
# each dataset with the direction fixed (higher score = responder);
# version 2 used pROC::roc() with its default direction = "auto", which
# flips any curve below the diagonal and so can only make AUC look better.
# =============================================================================
source("analysis/analysis_functions.R")
res <- load_results("heterogeneity")
MCOL <- c("MIMOSA2 (pooled fit)" = "deeppink", "MIMOSA2 (separate fit)" = "steelblue3", "DiD" = "grey30")
mlab <- function(m) recode(m, MIMOSA2_pooled = "MIMOSA2 (pooled fit)", MIMOSA2_separate = "MIMOSA2 (separate fit)", DiD = "DiD")

perf <- performance_table(res) %>%
  filter(Rule %in% c("BFDR", "BH")) %>%
  mutate(Model = mlab(Method),
         Split = factor(sprintf("%d focal : %d other", P_focal, P_other),
                        levels = sprintf("%d focal : %d other", c(20, 50, 80), c(80, 50, 20))),
         Cell_range = factor(Cell_range, levels = CELL_LEVELS))
eff_other <- sort(unique(perf$Effect_other))

plot_measure <- function(d, y, se, ylab, title, ylim, hline = NULL) {
  # [CHANGE v3.1] discrete x axis: on a log axis 6.25e-4, 8e-4 and 1e-3 sat on top of each other
  d <- d %>% mutate(Eff_f = factor(Effect_other, levels = eff_other, labels = effect_lab(eff_other)))
  p <- ggplot(d, aes(Eff_f, .data[[y]], colour = Model, fill = Model, group = Model)) +
    geom_ribbon(aes(ymin = pmax(.data[[y]] - 1.96 * .data[[se]], ylim[1]),
                    ymax = pmin(.data[[y]] + 1.96 * .data[[se]], ylim[2])), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
    scale_colour_manual(values = MCOL, name = NULL) + scale_fill_manual(values = MCOL, name = NULL) +
    coord_cartesian(ylim = ylim) +
    facet_grid(Cell_range ~ Split) +
    labs(title = title, x = "Effect size of the other subgroup", y = ylab,
         caption = "Measured on the focal (weak-effect) subjects only. Bands: +/- 1.96 Monte Carlo SE.") +
    theme_mimosa()
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed")
  p
}

for (rp in sort(unique(perf$Res_prop))) for (ef in sort(unique(perf$Effect_focal))) {
  d <- perf %>% filter(Res_prop == rp, Effect_focal == ef)
  tag <- sprintf("rho%02d_focal%s", round(100 * res_rho(rp)), gsub("[.+-]", "", effect_lab(ef)))
  sub <- sprintf("%s, focal effect %s.", res_lab(rp), effect_lab(ef))
  save_fig(plot_measure(d %>% filter(Alpha == 0.01), "AUC", "AUC_mcse", "Mean per-dataset AUC (focal subjects)",
                        "Effect-size heterogeneity: discrimination of weak responders.", c(0.4, 1)) + labs(subtitle = sub),
           paste0("het_auc_", tag), 10, 8)
  save_fig(plot_measure(d %>% filter(Alpha == 0.01), "TPR", "TPR_mcse", "TPR at 1% nominal FDR (focal subjects)",
                        "Effect-size heterogeneity: sensitivity for weak responders.", c(0, 1)) + labs(subtitle = sub),
           paste0("het_tpr_", tag), 10, 8)
  save_fig(plot_measure(d %>% filter(Alpha == 0.05), "FDR", "FDR_mcse", "Observed FDR at 5% (focal subjects)",
                        "Effect-size heterogeneity: false discoveries among weak responders.", c(0, 0.4), hline = 0.05) +
             labs(subtitle = sub), paste0("het_fdr_", tag), 10, 8)
}

# ---- Paired difference: pooled - separate ------------------------------------------
dd <- paired_auc_diff(res, "MIMOSA2_pooled", "MIMOSA2_separate") %>%
  mutate(Split = factor(sprintf("%d focal : %d other", P_focal, P_other),
                        levels = sprintf("%d focal : %d other", c(20, 50, 80), c(80, 50, 20))),
         Cell_range = factor(Cell_range, levels = CELL_LEVELS),
         Row_lab = paste0("focal ", effect_lab(Effect_focal), ", ", res_lab(Res_prop)))
save_tab(dd, "heterogeneity_dAUC_pooled_minus_separate")
p <- ggplot(dd %>% mutate(Eff_f = factor(Effect_other, levels = eff_other, labels = effect_lab(eff_other))),
            aes(Eff_f, dAUC_mean, colour = Cell_range)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.4), size = 0.25) +
  scale_colour_manual(values = CELL_COLOURS, name = "Cell count") +
  facet_grid(Row_lab ~ Split) +
  labs(title = "Cost of sharing information with strong responders.",
       subtitle = "Paired difference in AUC on the focal subjects: pooled fit - separate fit.",
       x = "Effect size of the other subgroup", y = "Difference in AUC (pooled - separate)",   # [v3.1] plain text: the Greek Delta did not render in the cluster PDFs
       caption = "Below 0 = pooling with strong responders hurts the weak responders. Bars: +/- 1.96 MCSE.") +
  theme_mimosa()
save_fig(p, "het_dAUC_pooled_vs_separate", 10, 3 + 2.5 * length(unique(dd$Row_lab)))

tab <- perf %>% filter(Alpha == 0.01) %>%
  arrange(Res_prop, Effect_focal, Cell_range, P_focal, Effect_other) %>%
  transmute(Res_prop, Effect_focal = effect_lab(Effect_focal), Cell_range, Split,
            Effect_other = effect_lab(Effect_other), Model, AUC = fmt_mcse(AUC, AUC_mcse),
            TPR = fmt_mcse(TPR, TPR_mcse), Fail = sprintf("%.1f%%", 100 * Fail_rate)) %>%
  pivot_wider(names_from = Model, values_from = c(AUC, TPR, Fail))
save_tab(tab, "heterogeneity_summary")
print(knitr::kable(tab, caption = "Effect heterogeneity, focal subjects, alpha = 1%: estimate (MCSE)"))
