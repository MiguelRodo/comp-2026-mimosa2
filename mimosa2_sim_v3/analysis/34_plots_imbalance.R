# =============================================================================
# 34_plots_imbalance.R : figures and tables for Study 4 (cell-count imbalance)
# v3.2 (14 Oct 2026): thesis figures only. See analysis_functions.R for the style.
# =============================================================================
# Figures made (placement agreed 14 Oct):
#   Section 4.4 : cell_imbalance_clean  (s1, u1, s0, u0 in a 2 x 2 grid)
#   Appendix    : cell_imbalance_plot   (all 14 depletion patterns)
# No longer made: imbalance_vs_balanced (the same differences are in the table
#   of 4.4) and auc_or_plot (the odds ratios stay in _tables/imbalance_summary.csv).
# Changes to the AUC figure:
#  * The x axis is categorical with evenly spaced fractions (on the log scale
#    75%, 90% and 100% were printed on top of each other).
#  * The reference (150,000 cells in every assay) is a horizontal line across
#    each panel; it was a single point at 100% before, so it was easy to miss.
# =============================================================================
source("analysis/analysis_functions.R")
res <- load_results("imbalance")
perf <- performance_table(res) %>%
  filter(Alpha == 0.01, Rule %in% c("BFDR", "unadjusted")) %>%
  mutate(Model = ifelse(Method == "DiD", "DiD", "MIMOSA2"))

ref <- perf %>% filter(Type == "Reference") %>%
  select(Model, AUC_ref = AUC, AUC_ref_mcse = AUC_mcse, TPR_ref = TPR, TPR_ref_mcse = TPR_mcse)
bal <- perf %>% filter(Type == "Balanced") %>%
  select(Model, Fraction, N_depleted, AUC_bal = AUC, AUC_bal_mcse = AUC_mcse, TPR_bal = TPR, TPR_bal_mcse = TPR_mcse)
imb <- perf %>% filter(Type == "Imbalanced") %>%
  left_join(bal, by = c("Model", "Fraction", "N_depleted")) %>%
  left_join(ref, by = "Model") %>%
  mutate(dAUC_vs_bal = AUC - AUC_bal, dAUC_vs_bal_mcse = sqrt(AUC_mcse^2 + AUC_bal_mcse^2),
         dTPR_vs_bal = TPR - TPR_bal, dTPR_vs_bal_mcse = sqrt(TPR_mcse^2 + TPR_bal_mcse^2),
         OR = (AUC / (1 - AUC)) / (AUC_ref / (1 - AUC_ref)),
         se_logOR = sqrt((AUC_mcse / (AUC * (1 - AUC)))^2 + (AUC_ref_mcse / (AUC_ref * (1 - AUC_ref)))^2),
         OR_lower = exp(log(OR) - 1.96 * se_logOR), OR_upper = exp(log(OR) + 1.96 * se_logOR),
         Pattern = factor(Depleted, levels = IMB_PATTERNS))
save_tab(imb, "imbalance_vs_balanced_and_reference")

# ---- AUC against the fraction kept: imbalanced, balanced control, reference -------------
fr <- sort(unique(imb$Fraction))
frac_f <- function(x) factor(x, levels = fr, labels = paste0(100 * fr, "%"))
pat_lab <- function(x) gsub("_", " + ", x)
SERIES <- c("Imbalanced", "Balanced (same total)", "Reference (150,000 in every assay)")
SERIES_COL <- setNames(c(PAL_GREEN, PAL_ORCHID, PAL_REF), SERIES)
mod_f <- function(m) plot_method(m, "unadjusted", auc = TRUE)
curves <- bind_rows(
  imb %>% transmute(Pattern, Model, Fraction, Series = SERIES[1], AUC, AUC_mcse),
  imb %>% distinct(Pattern, Model, Fraction, AUC_bal, AUC_bal_mcse) %>%
    transmute(Pattern, Model, Fraction, Series = SERIES[2], AUC = AUC_bal, AUC_mcse = AUC_bal_mcse)) %>%
  mutate(Series = factor(Series, levels = SERIES), Label = mod_f(Model), Frac = frac_f(Fraction))
refl <- imb %>% distinct(Pattern, Model, AUC_ref) %>%
  mutate(Series = factor(SERIES[3], levels = SERIES), Label = mod_f(Model))

auc_plot <- function(pats, ncol) {
  cv <- curves %>% filter(Pattern %in% pats); rl <- refl %>% filter(Pattern %in% pats)
  ggplot(cv, aes(Frac, AUC, colour = Series, linetype = Label, group = interaction(Series, Label))) +
    geom_hline(data = rl, aes(yintercept = AUC_ref, colour = Series, linetype = Label), linewidth = 0.6) +
    geom_ribbon(aes(ymin = AUC - 1.96 * AUC_mcse, ymax = AUC + 1.96 * AUC_mcse, fill = Series),
                alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
    scale_colour_manual(values = SERIES_COL, name = NULL, drop = FALSE) +
    scale_fill_manual(values = SERIES_COL, name = NULL, drop = FALSE, guide = "none") +
    scale_method() +
    facet_wrap(~ Pattern, ncol = ncol, labeller = labeller(Pattern = pat_lab)) +
    labs(x = "Depleted assay cell count", y = "Average AUC") +
    theme_mimosa() + guides(colour = guide_legend(order = 1, nrow = 2), linetype = guide_legend(order = 2))
}
# 4.4: the four single-assay patterns
save_fig(auc_plot(c("s1", "u1", "s0", "u0"), ncol = 2), "cell_imbalance_clean", height = 7)
# Appendix: all patterns (x labels tilted: seven labels in a narrow panel)
save_fig(auc_plot(levels(imb$Pattern), ncol = 4) + tilt_x(), "cell_imbalance_plot", height = 10.5)

# ---- Summary table (estimate (MCSE)), incl. the AUC odds ratio vs the reference -------------
tab <- imb %>% arrange(Pattern, desc(Fraction)) %>%
  transmute(Pattern, Fraction, Model, AUC = fmt_mcse(AUC, AUC_mcse), AUC_balanced = fmt_mcse(AUC_bal, AUC_bal_mcse),
            dAUC_vs_balanced = fmt_mcse(dAUC_vs_bal, dAUC_vs_bal_mcse),
            OR_vs_ref = sprintf("%.2f [%.2f, %.2f]", OR, OR_lower, OR_upper), Fail = sprintf("%.1f%%", 100 * Fail_rate))
save_tab(tab, "imbalance_summary")
print(knitr::kable(tab, caption = "Cell-count imbalance: estimate (MCSE)"))
