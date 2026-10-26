# =============================================================================
# 32_plots_prior.R : figures and tables for Study 2 (generative distribution)
# v3.2 (14 Oct 2026): thesis figures only. See analysis_functions.R for the style.
# =============================================================================
# Figures made (placement agreed 14 Oct):
#   Chapter 3   : prior_density
#   Section 4.2 : prior_plot_P50, prior_auc_P50, prior_dAUC_vs_beta_P50
#   Appendix    : prior_plot_P20/P100, prior_auc_P20/P100, prior_dAUC_vs_beta_P20/P100,
#                 fdr_prior_plot_P50 (recommended: the FDR result is quoted in 4.2)
# No longer made: prior_v2_vs_v3 (and its moments table), prior_ROC_comp_plot,
#   fdr_prior_plot_P20/P100.
# =============================================================================
source("analysis/analysis_functions.R")
source("R/dgm.R")
FAM_LEVELS <- c("Beta", "SX", "BB", "EG", "LN")
FAM_LABELS <- c("Beta", "Simplex-type", "Bimodal Beta", "Exponential-gamma", "Logit-normal")
fam_factor <- function(x) factor(x, levels = FAM_LEVELS, labels = FAM_LABELS)
# [v3.4] FAM_COLOURS is now defined in analysis_functions.R (palette v2)
stopifnot(all(FAM_LABELS %in% names(FAM_COLOURS)))

# Legend inside the figure, in the empty panel slot (x, y in 0-1 figure units)
legend_inside <- function(x, y) {
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    theme(legend.position = "inside", legend.position.inside = c(x, y),
          legend.justification = c(0.5, 0.5), legend.box = "vertical", legend.direction = "vertical")
  else theme(legend.position = c(x, y), legend.box = "vertical", legend.direction = "vertical")
}
# Five family panels in a 2 x 3 grid with the EMPTY slot at the TOP RIGHT, so the
# legend sits there and the bottom row holds only panels (with the x-axis labels).
# facet_wrap(as.table = FALSE) fills the grid from the bottom-left, so the levels
# are given in the order bottom row (left to right), then top row.
FAM_GRID_ORDER <- c("Bimodal Beta", "Exponential-gamma", "Logit-normal", "Beta", "Simplex-type")
fam_grid <- function(x) factor(as.character(x), levels = FAM_GRID_ORDER)

# ---- 0. Density figure (Chapter 3; does not need simulation results) ----------------
# Each panel uses the condition it is labelled with (v3.1 fix). Layout 1 + 2: u0 / u1
# top left, legend top right, s0 and s1 in the bottom row (as.table = FALSE fills the
# grid from the bottom-left, hence the level order s0, s1, u0 / u1).
# Axis numbers in scientific notation, not tilted.
set.seed(2026)
panels <- data.frame(
  cond = c("u0", "s0", "s1"),
  mu   = c(1e-4, 3.5e-4, 1.1e-3),
  lab  = c("u[0] / u[1]", "s[0]", "s[1]"))          # [v3.3] subscripts, as in the text
panels$strip <- sprintf('%s ~ "(" * mu == %s * ")"', panels$lab, sci_math(panels$mu))
dens <- bind_rows(lapply(seq_len(nrow(panels)), function(i) {
  # [v3.3, R7] BB is only used for stimulated assays; in u0/u1 draw_prop() falls
  # back to Beta, so the BB curve there was a duplicate of Beta. Leave it out.
  fams <- if (panels$cond[i] %in% c("s0", "s1")) FAM_LEVELS else setdiff(FAM_LEVELS, "BB")
  bind_rows(lapply(fams, function(f) {
    x <- draw_prop(4e5, panels$mu[i], panels$cond[i], f, 10000)
    dd <- density(x, from = 0, to = quantile(x, 0.995), n = 512)
    data.frame(Family = f, x = dd$x, y = dd$y, Panel = panels$strip[i], mu = panels$mu[i])
  }))
})) %>% mutate(Family = fam_factor(Family), Panel = factor(Panel, levels = panels$strip[c(2, 3, 1)]))
p <- ggplot(dens, aes(x, y, colour = Family)) +
  geom_vline(data = distinct(dens, Panel, mu), aes(xintercept = mu),
             linetype = "dashed", colour = PAL_MU, linewidth = 0.4, show.legend = FALSE) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Panel, ncol = 2, scales = "free", labeller = label_parsed, as.table = FALSE) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3), labels = sci_axis(2)) +
  scale_y_continuous(breaks = NULL) +                 # [v3.3, R7] density values carry no meaning here
  scale_colour_manual(values = FAM_COLOURS, name = "Generative distribution") +
  guides(colour = guide_legend(override.aes = list(linewidth = 1))) +
  labs(x = "Proportion of positive cells", y = "Density") +   # [v3.3] key term
  theme_mimosa() + legend_inside(0.76, 0.77)
save_fig(p, "prior_density", height = 6.5)

# ---- Simulation results ---------------------------------------------------------------
res <- load_results("prior")
perf <- performance_table(res) %>%
  mutate(Fam = fam_factor(Distribution), Fam_grid = fam_grid(Fam),
         Cell_range = factor(Cell_range, levels = CELL_LEVELS))
effects <- sort(unique(perf$Effect))
calls <- perf %>% filter((Method == "MIMOSA2" & Rule == "BFDR") | (Method == "DiD" & Rule == "BH")) %>%
  mutate(Label = plot_method(Method, Rule))
aucs  <- perf %>% filter(Alpha == 0.01, (Method == "MIMOSA2" & Rule == "BFDR") | (Method == "DiD" & Rule == "unadjusted")) %>%
  mutate(Label = plot_method(Method, Rule, auc = TRUE))

band <- function(d, y, se, ylab, ylim = c(0, 1), hline = NULL) {
  p <- ggplot(d, aes(Effect, .data[[y]], colour = Cell_range, fill = Cell_range, linetype = Label,
                     group = interaction(Cell_range, Label))) +
    geom_ribbon(aes(ymin = pmax(.data[[y]] - 1.96 * .data[[se]], ylim[1]),
                    ymax = pmin(.data[[y]] + 1.96 * .data[[se]], ylim[2])), alpha = BAND_ALPHA, colour = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
    scale_x_log10(breaks = effects, labels = effect_axis) +
    scale_cell() + scale_method() +
    coord_cartesian(ylim = ylim) +
    facet_wrap(~ Fam_grid, ncol = 3, as.table = FALSE) +
    labs(x = X_EFFECT, y = ylab) +
    theme_mimosa() + legend_inside(0.84, 0.76)          # [v3.3] untilted (R6)
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", colour = PAL_REFLINE)
  p
}
for (pp in sort(unique(perf$P))) {
  save_fig(band(calls %>% filter(P == pp, Alpha == 0.01), "TPR", "TPR_mcse", "Sensitivity"),
           sprintf("prior_plot_P%d", pp), height = 6.5)
  save_fig(band(aucs %>% filter(P == pp), "AUC", "AUC_mcse", "Average AUC", ylim = c(0.5, 1)),
           sprintf("prior_auc_P%d", pp), height = 6.5)
}
P_FDR <- if (50 %in% perf$P) 50 else max(perf$P)
save_fig(band(calls %>% filter(P == P_FDR, Alpha == 0.05), "FDR", "FDR_mcse", "Observed FDR",
              ylim = c(0, 0.3), hline = 0.05),
         sprintf("fdr_prior_plot_P%d", P_FDR), height = 6.5)

# ---- Difference from the Beta reference (same P, effect, cell count) -------------------
# Independent simulations, so MCSE = sqrt(MCSE_1^2 + MCSE_2^2). Rows: model; columns: family.
ref <- perf %>% filter(Distribution == "Beta") %>%
  select(P, Effect, Cell_range, Method, Rule, Alpha, AUC_ref = AUC, AUC_ref_mcse = AUC_mcse,
         TPR_ref = TPR, TPR_ref_mcse = TPR_mcse)
vs <- perf %>% filter(Distribution != "Beta", Alpha == 0.01, Rule %in% c("BFDR", "unadjusted")) %>%
  left_join(ref, by = c("P", "Effect", "Cell_range", "Method", "Rule", "Alpha")) %>%
  mutate(dAUC = AUC - AUC_ref, dAUC_mcse = sqrt(AUC_mcse^2 + AUC_ref_mcse^2),
         dTPR = TPR - TPR_ref, dTPR_mcse = sqrt(TPR_mcse^2 + TPR_ref_mcse^2),
         Model = ifelse(Method == "DiD", "DiD", "MIMOSA2"),        # table column (unchanged)
         Label = plot_method(Method, Rule, auc = TRUE))
save_tab(vs, "prior_difference_from_beta")
for (pp in sort(unique(vs$P))) {
  p <- ggplot(vs %>% filter(P == pp), aes(Effect, dAUC, colour = Cell_range)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = PAL_REFLINE) +
    geom_pointrange(aes(ymin = dAUC - 1.96 * dAUC_mcse, ymax = dAUC + 1.96 * dAUC_mcse),
                    position = position_dodge(width = 0.2), size = 0.2, linewidth = 0.4) +
    scale_x_log10(breaks = effects, labels = effect_axis) +
    scale_cell(aesthetics = "colour") +
    facet_grid(Label ~ Fam) +
    labs(x = X_EFFECT, y = DAUC_LAB("alternative - Beta")) +
    theme_mimosa()
  save_fig(p, sprintf("prior_dAUC_vs_beta_P%d", pp), height = 6)
}

# ---- Summary table (estimate (MCSE)) -----------------------------------------------------
P_SHOW <- if (50 %in% perf$P) 50 else max(perf$P)
tab <- perf %>% filter(P == P_SHOW, Alpha == 0.01, Rule %in% c("BFDR", "unadjusted")) %>%
  arrange(Fam, desc(Effect), Cell_range) %>%
  transmute(Distribution = Fam, Effect = effect_lab(Effect), Cell_range,
            Model = ifelse(Method == "DiD", "DiD", "MIMOSA2"), AUC = fmt_mcse(AUC, AUC_mcse)) %>%
  pivot_wider(names_from = Model, values_from = AUC)
save_tab(tab, sprintf("prior_auc_table_P%d", P_SHOW))
print(knitr::kable(tab, caption = sprintf("Mean per-dataset AUC (MCSE), P = %d", P_SHOW)))
