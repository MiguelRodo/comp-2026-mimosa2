# =============================================================================
# 33_plots_heterogeneity.R : figures and tables for Study 3 (effect heterogeneity)
# v3.2 (14 Oct 2026): replaces the 18 per-scenario figures (het_auc/tpr/fdr_*)
# and het_dAUC_pooled_vs_separate with two summary figures and one table.
# =============================================================================
# Performance is measured on the FOCAL (weak-effect) subjects only. The key
# comparison is the PAIRED difference pooled fit - separate fit on the same datasets.
#
#   het_dAUC_cost_focal5e04    (4.3)  Cost of pooling: Delta AUC (pooled - separate)
#                                     against cell count, one line per delta_o;
#                                     rows rho = 0.2, 0.5, 0.8; columns = split.
#   het_dAUC_cost_focal625e04  (App.) The same for delta_f = 6.25e-4 (almost identical).
#   het_auc_levels             (4.3)  AUC of the pooled fit, separate fit and DiD
#                                     against cell count, for mild (2x) and strong
#                                     (100x) heterogeneity; rho = 0.5, delta_f = 5e-4.
#   _tables/heterogeneity_dAUC_wide_rho50_focal50e04.csv  table for 4.3.
# =============================================================================
source("analysis/analysis_functions.R")
res <- load_results("heterogeneity")
SPLIT_LEVELS <- sprintf("%d focal : %d other", c(20, 50, 80), c(80, 50, 20))
add_split <- function(d) d %>% mutate(Split = factor(sprintf("%d focal : %d other", P_focal, P_other),
                                                     levels = SPLIT_LEVELS),
                                      Cell_range = factor(Cell_range, levels = CELL_LEVELS))
CELL_X <- scale_x_discrete(labels = CELL_LABELS)
RHO_SHOW <- "Prop_0.50"; EF_SHOW <- 5e-4

# Colour per value of delta_o: blues = mild (at most 2x delta_f), oranges = strong
DO_COLOURS <- c("6.25e-04" = "#9ECAE1", "8.00e-04" = "#6BAED6", "1.00e-03" = "#2171B5",
                "1.00e-02" = "#E6550D", "5.00e-02" = "#7F2704")
do_key <- function(x) sprintf("%.2e", x)
do_labels <- function(eo, ef) {                       # "6.25 x 10^-4 (1.25 delta_f)"
  parse(text = sprintf('%s ~ "(" * %s * delta[f] * ")"', sci_math(eo), signif(eo / ef, 3)))
}
DO_TITLE <- expression("Effect of the other subgroup, " * delta[o])

# ---- 1. Cost of pooling ---------------------------------------------------------------
dd <- paired_auc_diff(res, "MIMOSA2_pooled", "MIMOSA2_separate") %>% add_split() %>%
  mutate(Res_lab = math_factor(Res_prop, rho_math), Eo = do_key(Effect_other))
save_tab(dd, "heterogeneity_dAUC_pooled_minus_separate")
for (ef in sort(unique(dd$Effect_focal))) {
  d <- dd %>% filter(Effect_focal == ef)
  eo <- sort(unique(d$Effect_other))
  p <- ggplot(d, aes(Cell_range, dAUC_mean, colour = Eo, group = Eo)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey20") +
    geom_line(position = position_dodge(width = 0.35), linewidth = 0.6) +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.35),
                    size = 0.2, linewidth = 0.4) +
    scale_colour_manual(values = DO_COLOURS, labels = do_labels(eo, ef), breaks = do_key(eo), name = DO_TITLE) +
    CELL_X +
    facet_grid(Res_lab ~ Split, labeller = labeller(Res_lab = label_parsed, Split = label_value)) +
    labs(x = CELL_TITLE, y = DAUC_LAB("pooled - separate")) +
    theme_mimosa() + tilt_x() + guides(colour = guide_legend(nrow = 2))
  save_fig(p, sprintf("het_dAUC_cost_focal%se04", gsub(".", "", format(ef * 1e4), fixed = TRUE)),
           height = 8.5)                 # het_dAUC_cost_focal5e04, het_dAUC_cost_focal625e04
}

# ---- 2. AUC levels: pooled vs separate vs DiD, mild and strong heterogeneity ------------
MCOL <- c("MIMOSA2 (pooled fit)" = PAL_ORCHID, "MIMOSA2 (separate fit)" = PAL_GREEN)
MCOL[LAB_DID_AUC] <- PAL_DID
MLTY <- setNames(c("solid", "solid", "dotted"), names(MCOL))
mlab <- function(m) factor(recode(m, MIMOSA2_pooled = "MIMOSA2 (pooled fit)",
                                  MIMOSA2_separate = "MIMOSA2 (separate fit)", DiD = LAB_DID_AUC),
                           levels = names(MCOL))
EO_SHOW <- c(1e-3, 5e-2)
perf <- performance_table(res) %>%
  filter(Alpha == 0.01, Rule %in% c("BFDR", "BH")) %>%     # AUC does not depend on the rule
  add_split() %>% mutate(Model = mlab(Method))
lv <- perf %>% filter(Res_prop == RHO_SHOW, Effect_focal == EF_SHOW,
                      sapply(Effect_other, function(x) any(abs(x - EO_SHOW) < 1e-12))) %>%
  mutate(Eo_lab = math_factor(Effect_other, function(x)
    sprintf('delta[o] == %s ~ "(" * %s * delta[f] * ")"', sci_math(x), signif(x / EF_SHOW, 3))))
p <- ggplot(lv, aes(Cell_range, AUC, colour = Model, fill = Model, linetype = Model, group = Model)) +
  geom_ribbon(aes(ymin = AUC - 1.96 * AUC_mcse, ymax = AUC + 1.96 * AUC_mcse), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  scale_colour_manual(values = MCOL, name = NULL) + scale_fill_manual(values = MCOL, name = NULL) +
  scale_linetype_manual(values = MLTY, name = NULL) +
  CELL_X + coord_cartesian(ylim = c(0.45, 1)) +
  facet_grid(Eo_lab ~ Split, labeller = labeller(Eo_lab = label_parsed, Split = label_value)) +
  labs(x = CELL_TITLE, y = "Average AUC (focal subjects)") +
  theme_mimosa() + tilt_x()
save_fig(p, "het_auc_levels", height = 6.5)

# ---- 3. Tables ---------------------------------------------------------------------------
# Wide table for 4.3: Delta AUC (pooled - separate), rho = 0.5, delta_f = 5e-4,
# rows cell count x split, columns delta_o. The largest MCSE goes in the caption.
wide <- dd %>% filter(Res_prop == RHO_SHOW, Effect_focal == EF_SHOW) %>%
  arrange(Cell_range, Split, Effect_other) %>%
  transmute(Cell_range, Split, Effect_other = effect_lab(Effect_other), dAUC = sprintf("%.3f", dAUC_mean)) %>%
  pivot_wider(names_from = Effect_other, values_from = dAUC)
save_tab(wide, "heterogeneity_dAUC_wide_rho50_focal50e04")
message("largest MCSE in that table: ",
        round(max(dd$dAUC_mcse[dd$Res_prop == RHO_SHOW & dd$Effect_focal == EF_SHOW]), 3))

tab <- perf %>%
  arrange(Res_prop, Effect_focal, Cell_range, P_focal, Effect_other) %>%
  transmute(Res_prop, Effect_focal = effect_lab(Effect_focal), Cell_range, Split,
            Effect_other = effect_lab(Effect_other), Model = recode(Method, MIMOSA2_pooled = "MIMOSA2 (pooled fit)",
            MIMOSA2_separate = "MIMOSA2 (separate fit)"), AUC = fmt_mcse(AUC, AUC_mcse),
            TPR = fmt_mcse(TPR, TPR_mcse), Fail = sprintf("%.1f%%", 100 * Fail_rate)) %>%
  pivot_wider(names_from = Model, values_from = c(AUC, TPR, Fail))
save_tab(tab, "heterogeneity_summary")
