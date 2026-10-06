# =============================================================================
# 32_plots_prior.R : figures and tables for Study 2 (generative distribution)
# Replaces Prior_simulation_plots.R, Variance_Solver.R (figures),
# Bimodal_prior.R and Wonky_prior.R (version 2).
# =============================================================================
# [CHANGE] What changed relative to Prior_simulation_plots.R:
#  * Its factor levels for Cell_range omitted "Extremely Low", so those
#    scenarios became NA in every figure; levels now come from scenarios.R.
#  * Results are shown for one P at a time (version 2 pooled P = 20, 50, 100
#    in the ROC curves and averaged over P in the AUROC tables; audit A6).
#  * Family effects are shown as DIFFERENCES from the Beta reference in the
#    same (P, effect, cell count), with MCSE = sqrt(MCSE_1^2 + MCSE_2^2)
#    (scenarios are simulated independently).
#  * The table marked "# INCORRECT" (method columns swapped) is replaced.
#  * The density figure is drawn from the samplers actually used in the
#    simulation (Bimodal_prior.R drew a DIFFERENT bimodal distribution, with
#    mean 6e-4, variance 2e-8 and a 50:50 mixture, from the one simulated).
#    A second panel shows the version-2 vs version-3 EG and LN distributions
#    at the stimulated responder mean, for the write-up.
# =============================================================================
source("analysis/analysis_functions.R")
source("R/dgm.R")
FAM_LEVELS <- c("Beta", "SX", "BB", "EG", "LN")
FAM_LABELS <- c("Beta", "Inverse-Gaussian odds (simplex-type)", "Bimodal Beta", "Exponential-gamma", "Logit-normal")
fam_factor <- function(x) factor(x, levels = FAM_LEVELS, labels = FAM_LABELS)
FAM_COLOURS <- setNames(c("black", "deeppink", "steelblue3", "orange", "purple"), FAM_LABELS)

# ---- 0. Density figure (does not need simulation results) -----------------------
set.seed(2026)
# [CHANGE v3.1, 6 Oct] v3.0 drew EVERY panel with cond = "s1", so the
# "unstimulated background" panel showed a bimodal Bimodal-Beta curve although
# the simulation draws unstimulated proportions from the Beta (dgm.R,
# draw_prop). Each panel now uses the condition it is labelled with, and the
# stimulated baseline mean (3.5e-4, where the bimodal family IS used) is added.
mu_list <- c("Mean 1e-4 (unstimulated, u0 / u1)" = 1e-4,
             "Mean 3.5e-4 (stimulated baseline, s0)" = 3.5e-4,
             "Mean 1.1e-3 (stimulated, responders, delta = 1e-3)" = 1.1e-3)
cond_list <- c("u0", "s0", "s1")
dens <- bind_rows(lapply(seq_along(mu_list), function(i) {
  nm <- names(mu_list)[i]; mu <- mu_list[[i]]
  bind_rows(lapply(FAM_LEVELS, function(f) {
    x <- draw_prop(4e5, mu, cond_list[i], f, 10000)
    dd <- density(x, from = 0, to = quantile(x, 0.995), n = 512)
    data.frame(Family = f, x = dd$x, y = dd$y, mean = mean(x), sd = sd(x), Mean_label = nm)
  }))
})) %>% mutate(Family = fam_factor(Family), Mean_label = factor(Mean_label, levels = names(mu_list)))
p <- ggplot(dens, aes(x, y, colour = Family)) + geom_line(linewidth = 0.8) +
  geom_vline(aes(xintercept = mean, colour = Family), linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ Mean_label, scales = "free") +
  scale_colour_manual(values = FAM_COLOURS, name = NULL) +
  labs(title = "Generative distributions of the cell proportions.",
       subtitle = "Every family has the mean and variance of Beta(phi = 10,000) at each condition mean.",
       x = "Proportion of cytokine-positive cells", y = "Density",
       caption = "4e5 draws per family from the samplers in R/dgm.R; dashed lines: means.") +
  theme_mimosa() + guides(colour = guide_legend(nrow = 2))
save_fig(p, "prior_density", 14, 5)

mu <- 1.1e-3
old_new <- bind_rows(
  data.frame(Version = "Version 2 (calibrated at 1e-4 only)", Family = "Exponential-gamma",
             x = exp(-rgamma(4e5, 113, 113 / (113 * (mu^(-1/113) - 1))))),
  data.frame(Version = "Version 2 (calibrated at 1e-4 only)", Family = "Logit-normal",
             x = plogis(rnorm(4e5, qlogis(mu), 1 / sqrt(2.05)))),
  data.frame(Version = "Version 3 (calibrated at each mean)", Family = "Exponential-gamma", x = draw_prop(4e5, mu, "s1", "EG", 10000)),
  data.frame(Version = "Version 3 (calibrated at each mean)", Family = "Logit-normal", x = draw_prop(4e5, mu, "s1", "LN", 10000)),
  data.frame(Version = "Reference", Family = "Beta", x = rbeta(4e5, mu * 1e4, (1 - mu) * 1e4)))
sumtab <- old_new %>% group_by(Version, Family) %>% summarise(mean = mean(x), sd = sd(x), .groups = "drop") %>%
  mutate(mean_ratio = mean / mu, sd_ratio = sd / sqrt(beta_var(mu, 1e4)))
save_tab(sumtab, "prior_v2_vs_v3_moments_at_MS1")
print(knitr::kable(sumtab, digits = 3, caption = "Moments at mu = 1.1e-3 (stimulated responder mean, delta = 1e-3)"))
p <- ggplot(old_new, aes(x, colour = Version)) + geom_density(linewidth = 0.7) +
  facet_wrap(~ Family, scales = "free_y") + coord_cartesian(xlim = c(0, 0.004)) +
  scale_colour_manual(values = c("Reference" = "black", "Version 2 (calibrated at 1e-4 only)" = "red3",
                                 "Version 3 (calibrated at each mean)" = "steelblue3"), name = NULL) +
  labs(title = "Why the prior-family study was re-run.",
       subtitle = "Stimulated responder mean 1.1e-3: version-2 EG and LN were 2-3x more variable than Beta.",
       x = "Proportion", y = "Density") + theme_mimosa() + guides(colour = guide_legend(nrow = 2))
save_fig(p, "prior_v2_vs_v3", 10, 5)

# ---- Simulation results -------------------------------------------------------------
res <- tryCatch(load_results("prior"), error = function(e) { message(conditionMessage(e)); NULL })
if (is.null(res)) stop("No prior-family results yet: the density figures above were made; run sims/21_sim_prior.R for the rest.")
perf <- performance_table(res) %>%
  mutate(Label = method_label(Method, Rule), Fam = fam_factor(Distribution),
         Cell_range = factor(Cell_range, levels = CELL_LEVELS),
         P_lab = factor(paste0("P = ", P), levels = paste0("P = ", sort(unique(P)))))
effects <- sort(unique(perf$Effect))
sel <- function(d) d %>% filter((Method == "MIMOSA2" & Rule == "BFDR") | (Method == "DiD" & Rule == "BH"))

band <- function(d, y, se, ylab, title, subtitle, ylim = c(0, 1), hline = NULL) {
  p <- ggplot(d, aes(Effect, .data[[y]], colour = Cell_range, fill = Cell_range, linetype = Label,
                     group = interaction(Cell_range, Label))) +
    geom_ribbon(aes(ymin = pmax(.data[[y]] - 1.96 * .data[[se]], ylim[1]),
                    ymax = pmin(.data[[y]] + 1.96 * .data[[se]], ylim[2])), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
    scale_x_log10(breaks = effects, labels = effect_lab(effects)) +
    scale_colour_manual(values = CELL_COLOURS, name = "Cell count") +
    scale_fill_manual(values = CELL_COLOURS, name = "Cell count") +
    scale_linetype_manual(values = c("MIMOSA2" = "solid", "DiD (BH)" = "dotted"), name = "Model") +
    coord_cartesian(ylim = ylim) +
    labs(title = title, subtitle = subtitle, x = "Effect size", y = ylab,
         caption = "Bands: +/- 1.96 Monte Carlo SE.") + theme_mimosa()
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed")
  p
}
for (pp in sort(unique(perf$P))) {
  d <- sel(perf) %>% filter(P == pp)
  save_fig(band(d %>% filter(Alpha == 0.01), "TPR", "TPR_mcse", "TPR", "Sensitivity analysis of MIMOSA2.",
                sprintf("TPR at 1%% nominal FDR, P = %d.", pp)) + facet_wrap(~ Fam, nrow = 1),
           sprintf("prior_plot_P%d", pp), 13, 4.5)
  save_fig(band(d %>% filter(Alpha == 0.05), "FDR", "FDR_mcse", "Observed FDR", "False discovery rate.",
                sprintf("Nominal 5%% (dashed), P = %d.", pp), ylim = c(0, 0.3), hline = 0.05) + facet_wrap(~ Fam, nrow = 1),
           sprintf("fdr_prior_plot_P%d", pp), 13, 4.5)
  dA <- perf %>% filter(P == pp, Alpha == 0.01, (Method == "MIMOSA2" & Rule == "BFDR") | (Method == "DiD" & Rule == "unadjusted")) %>%
    mutate(Label = ifelse(Method == "DiD", "DiD (BH)", "MIMOSA2"))
  save_fig(band(dA, "AUC", "AUC_mcse", "Mean per-dataset AUC", "Discrimination by generative distribution.",
                sprintf("P = %d.", pp), ylim = c(0.5, 1)) + facet_wrap(~ Fam, nrow = 1),
           sprintf("prior_auc_P%d", pp), 13, 4.5)
}

# ---- Difference from the Beta reference (same P, effect, cell count) -----------------
ref <- perf %>% filter(Distribution == "Beta") %>%
  select(P, Effect, Cell_range, Method, Rule, Alpha, AUC_ref = AUC, AUC_ref_mcse = AUC_mcse,
         TPR_ref = TPR, TPR_ref_mcse = TPR_mcse)
vs <- perf %>% filter(Distribution != "Beta", Alpha == 0.01, Rule %in% c("BFDR", "unadjusted")) %>%
  left_join(ref, by = c("P", "Effect", "Cell_range", "Method", "Rule", "Alpha")) %>%
  mutate(dAUC = AUC - AUC_ref, dAUC_mcse = sqrt(AUC_mcse^2 + AUC_ref_mcse^2),
         dTPR = TPR - TPR_ref, dTPR_mcse = sqrt(TPR_mcse^2 + TPR_ref_mcse^2),
         Model = ifelse(Method == "DiD", "DiD", "MIMOSA2"))
save_tab(vs, "prior_difference_from_beta")
for (pp in sort(unique(vs$P))) {
  p <- ggplot(vs %>% filter(P == pp), aes(Effect, dAUC, colour = Cell_range, shape = Model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(ymin = dAUC - 1.96 * dAUC_mcse, ymax = dAUC + 1.96 * dAUC_mcse),
                    position = position_dodge(width = 0.2), size = 0.25) +
    scale_x_log10(breaks = effects, labels = effect_lab(effects)) +
    scale_colour_manual(values = CELL_COLOURS, name = "Cell count") +
    facet_wrap(~ Fam, nrow = 1) +
    labs(title = "Change in AUC relative to the Beta distribution.", subtitle = sprintf("P = %d.", pp),
         x = "Effect size", y = expression(AUC[family] - AUC[Beta]),
         caption = "Bars: +/- 1.96 MCSE. Below 0 = the method does worse than under its own (Beta) assumption.") +
    theme_mimosa()
  save_fig(p, sprintf("prior_dAUC_vs_beta_P%d", pp), 12, 4.5)
}

# ---- Vertically averaged ROC curves, Beta vs each family ---------------------------
P_SHOW <- if (50 %in% perf$P) 50 else max(perf$P)
E_SHOW <- effects[ceiling(length(effects) / 2)]
keep <- res$design %>% filter(P == P_SHOW, Effect == E_SHOW) %>% pull(Scenario_ID) %>% unique()
roc <- roc_curves(res, c(MIMOSA2 = "MIMOSA2_prob", DiD = "DiD_z"), keep_scenarios = keep) %>%
  mutate(Fam = fam_factor(Distribution), Cell_range = factor(Cell_range, levels = CELL_LEVELS))
p <- ggplot(roc, aes(FPR, TPR_mean, colour = Fam, linetype = Method, group = interaction(Fam, Method))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey10", linewidth = 0.4) +
  geom_line(linewidth = 0.7) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +   # [v3.1] full 0-1 axes
  scale_colour_manual(values = FAM_COLOURS, name = "Distribution") +
  scale_linetype_manual(values = c(MIMOSA2 = "solid", DiD = "dotted"), name = "Model") +
  facet_wrap(~ Cell_range, nrow = 1) +
  labs(title = "ROC analysis across generative distributions.",
       subtitle = sprintf("Vertically averaged ROC curves, P = %d, effect %s.", P_SHOW, effect_lab(E_SHOW)),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_mimosa() + guides(colour = guide_legend(nrow = 2))
save_fig(p, "prior_ROC_comp_plot", 13, 5)

# ---- Summary table (estimate (MCSE)) --------------------------------------------------
tab <- perf %>% filter(P == P_SHOW, Alpha == 0.01, Rule %in% c("BFDR", "unadjusted")) %>%
  arrange(Fam, desc(Effect), Cell_range) %>%
  transmute(Distribution = Fam, Effect = effect_lab(Effect), Cell_range,
            Model = ifelse(Method == "DiD", "DiD", "MIMOSA2"), AUC = fmt_mcse(AUC, AUC_mcse)) %>%
  pivot_wider(names_from = Model, values_from = AUC)
save_tab(tab, sprintf("prior_auc_table_P%d", P_SHOW))
print(knitr::kable(tab, caption = sprintf("Mean per-dataset AUC (MCSE), P = %d", P_SHOW)))
