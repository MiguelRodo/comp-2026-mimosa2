# =============================================================================
# 34_plots_imbalance.R : figures and tables for Study 4 (cell-count imbalance)
# Replaces Count_Imbalanced_Plot.R (version 2).
# =============================================================================
# [CHANGE] What changed relative to Count_Imbalanced_Plot.R:
#  * BUG: scenario names end in "_Imb_s1", "_Bal_1", ... but the script tested
#    grepl("_Imb$") and grepl("_Bal$"). In the ROC figure every scenario was
#    therefore labelled "Balanced"; in the AUC odds-ratio figure every
#    non-reference scenario got Balance = NA and was dropped. The scenario
#    type, depleted assays and matched control are now explicit columns.
#  * AUC is computed per dataset and averaged; its MCSE replaces the DeLong SE
#    of a pooled AUC (DeLong treats all pooled subjects as independent, but
#    subjects in one dataset share one fitted model, so that SE is too small).
#  * Each imbalanced scenario is compared with its MATCHED balanced control
#    (same total cells: k depleted assays <-> Balanced k), which isolates the
#    effect of imbalance at a fixed budget, as Chapter 3 describes.
#  * The AUC odds ratio vs the reference is kept, with a delta-method CI
#    built from the two MCSEs.
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
  mutate(# effect of imbalance at a fixed total number of cells
         dAUC_vs_bal = AUC - AUC_bal, dAUC_vs_bal_mcse = sqrt(AUC_mcse^2 + AUC_bal_mcse^2),
         dTPR_vs_bal = TPR - TPR_bal, dTPR_vs_bal_mcse = sqrt(TPR_mcse^2 + TPR_bal_mcse^2),
         # AUC odds ratio vs the 150,000-cell reference (version-2 measure) + delta-method CI
         OR = (AUC / (1 - AUC)) / (AUC_ref / (1 - AUC_ref)),
         se_logOR = sqrt((AUC_mcse / (AUC * (1 - AUC)))^2 + (AUC_ref_mcse / (AUC_ref * (1 - AUC_ref)))^2),
         OR_lower = exp(log(OR) - 1.96 * se_logOR), OR_upper = exp(log(OR) + 1.96 * se_logOR),
         Pattern = factor(Depleted, levels = IMB_PATTERNS))
save_tab(imb, "imbalance_vs_balanced_and_reference")
fr <- sort(unique(imb$Fraction))   # depleted fractions only (1 = reference added separately)

# ---- 1. AUC vs fraction: imbalanced, matched balanced, reference ----------------
curves <- bind_rows(
  imb %>% transmute(Pattern, Model, Fraction, Series = "Imbalanced", AUC, AUC_mcse),
  imb %>% distinct(Pattern, Model, Fraction, AUC_bal, AUC_bal_mcse) %>%
    transmute(Pattern, Model, Fraction, Series = "Balanced (same total)", AUC = AUC_bal, AUC_mcse = AUC_bal_mcse),
  imb %>% distinct(Pattern, Model, AUC_ref, AUC_ref_mcse) %>%
    transmute(Pattern, Model, Fraction = 1, Series = "Reference (150k each)", AUC = AUC_ref, AUC_mcse = AUC_ref_mcse))
p <- ggplot(curves, aes(Fraction, AUC, colour = Series, linetype = Model, group = interaction(Series, Model))) +
  geom_ribbon(aes(ymin = AUC - 1.96 * AUC_mcse, ymax = AUC + 1.96 * AUC_mcse, fill = Series), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  scale_x_log10(breaks = c(fr, 1), labels = c(paste0(100 * fr, "%"), "100%")) +
  scale_colour_manual(values = c("Imbalanced" = "darkgreen", "Balanced (same total)" = "orchid", "Reference (150k each)" = "black"), name = NULL) +
  scale_fill_manual(values = c("Imbalanced" = "darkgreen", "Balanced (same total)" = "orchid", "Reference (150k each)" = "black"), name = NULL) +
  scale_linetype_manual(values = c(MIMOSA2 = "solid", DiD = "dotted"), name = "Model") +
  facet_wrap(~ Pattern, labeller = labeller(Pattern = function(x) paste("depleted:", gsub("_", " + ", x)))) +
  labs(title = "MIMOSA2 performance under parent cell count imbalance.",
       subtitle = "Mean per-dataset AUC; depleted assay(s) reduced to a fraction of 150,000 cells.",
       x = "Depleted assay cell count (fraction of 150,000)", y = "Mean per-dataset AUC",
       caption = "Balanced: the same total number of cells spread equally over the four assays. Bands: +/- 1.96 MCSE.") +
  theme_mimosa()
save_fig(p, "cell_imbalance_plot", 11, 3 + 2.6 * ceiling(length(unique(imb$Pattern)) / 3))

# ---- 2. Effect of imbalance at a fixed total (imbalanced - balanced) ---------------
p <- ggplot(imb, aes(Fraction, dAUC_vs_bal, colour = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = dAUC_vs_bal - 1.96 * dAUC_vs_bal_mcse, ymax = dAUC_vs_bal + 1.96 * dAUC_vs_bal_mcse),
                  position = position_dodge(width = 0.08), size = 0.25) +
  scale_x_log10(breaks = fr, labels = paste0(100 * fr, "%")) +
  scale_colour_manual(values = c(MIMOSA2 = "deeppink", DiD = "grey30"), name = "Model") +
  facet_wrap(~ Pattern, labeller = labeller(Pattern = function(x) paste("depleted:", gsub("_", " + ", x)))) +
  labs(title = "Cost of imbalance at a fixed total number of cells.",
       subtitle = "AUC(imbalanced) - AUC(balanced control with the same total).",
       x = "Depleted assay cell count (fraction of 150,000)", y = "Difference in AUC (imbalanced - balanced)",   # [v3.1] plain text: the Greek Delta did not render in the cluster PDFs
       caption = "Below 0 = imbalance itself (not just fewer cells) reduces discrimination. Bars: +/- 1.96 MCSE.") +
  theme_mimosa()
save_fig(p, "imbalance_vs_balanced", 11, 3 + 2.6 * ceiling(length(unique(imb$Pattern)) / 3))

# ---- 3. AUC odds ratio vs reference (version-2 figure, corrected) ------------------
p <- ggplot(imb, aes(Fraction, OR, colour = Pattern, linetype = Model, group = interaction(Pattern, Model))) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(linewidth = 0.6) +
  geom_pointrange(aes(ymin = OR_lower, ymax = OR_upper), size = 0.2) +
  scale_x_log10(breaks = fr, labels = paste0(100 * fr, "%")) + scale_y_log10() +
  scale_linetype_manual(values = c(MIMOSA2 = "solid", DiD = "dotted"), name = "Model") +
  labs(title = "MIMOSA2 performance under parent cell count imbalance.",
       subtitle = "AUC odds ratio relative to the 150,000-cell reference.",
       x = "Depleted assay cell count (fraction of 150,000)", y = "AUC odds ratio (log scale)", colour = "Depleted",
       caption = "OR = [AUC/(1-AUC)] / [AUC_ref/(1-AUC_ref)]; 95% CI by the delta method from the Monte Carlo SEs.") +
  theme_mimosa()
save_fig(p, "auc_or_plot", 9, 6)

tab <- imb %>% arrange(Pattern, desc(Fraction)) %>%
  transmute(Pattern, Fraction, Model, AUC = fmt_mcse(AUC, AUC_mcse), AUC_balanced = fmt_mcse(AUC_bal, AUC_bal_mcse),
            dAUC_vs_balanced = fmt_mcse(dAUC_vs_bal, dAUC_vs_bal_mcse),
            OR_vs_ref = sprintf("%.2f [%.2f, %.2f]", OR, OR_lower, OR_upper), Fail = sprintf("%.1f%%", 100 * Fail_rate))
save_tab(tab, "imbalance_summary")
print(knitr::kable(tab, caption = "Cell-count imbalance: estimate (MCSE)"))
