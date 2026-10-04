# =============================================================================
# 40_legacy_plots.R : plots of the LEGACY (version-2) output, in your style
#   Count_Imbalance.Rdata            (Count_Imbalanced_Sim.R, 120 scenarios x 100)
#   EffectSize_Heterogeneity.Rdata   (Heterogeneous_effect_sim.R, 20,250 datasets)
# =============================================================================
# Same look as Count_Imbalanced_Plot.R and Effect_heterogeneity_plot.R
# (theme_bw, orchid / darkgreen, deeppink / steelblue3 / orange, ROC panels,
# AUC odds-ratio plot), with these corrections:
#
#  * Scenario labels are parsed from the real names ("_Imb_s1", "_Bal_1", ...).
#    The old grepl("_Imb$") / grepl("_Bal$") matched nothing, so every scenario
#    was "Balanced" in the ROC plot and dropped from the OR plot.
#  * Each imbalanced pattern is compared with its MATCHED balanced control
#    (k depleted assays <-> "_Bal_k", same total number of cells).
#  * ROC curves are averaged over datasets (TPR of every dataset at a fixed
#    grid of FPR values, then the mean; "vertical averaging"), and AUC is the
#    mean of the per-dataset AUCs, with a Monte Carlo SE. The old code pooled
#    all subjects of all 100 datasets into one ROC curve, and the DeLong SE
#    treated them as independent (too small).
#  * Heterogeneity: in version 2 the column "Small_effect" holds the LARGER
#    effect, and that group is the one evaluated (bug in the old sim script).
#    Columns are renamed Eval_effect / Other_effect so the plots say what was
#    actually simulated.
#  * Caveat printed on every figure: version-2 DiD scores gave the tie value
#    0.05 to every negative / non-converged subject (33-48% of subjects), which
#    pulls the DiD curves and AUCs down.
#
# Run locally in RStudio (only needs dplyr, tidyr, ggplot2):
#   setwd("C:/Github/comp-2026-mimosa2")      # folder that contains _simulations/
#   source("path/to/mimosa2_sim_v3/analysis/40_legacy_plots.R")
# Figures go to _fig/legacy/, tables to _tables/legacy/.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

# ---- settings you may want to change ------------------------------------------
IMB_FILE <- "_simulations/Count_Imbalance.Rdata"
HET_FILE <- "_simulations/EffectSize_Heterogeneity.Rdata"
FIG <- "_fig/legacy"; TAB <- "_tables/legacy"
IMB_PATTERN <- "s1"            # depleted assay(s) for the main ROC figure (also loops over all below)
HET_RHO     <- "Prop_0.50"     # heterogeneity ROC figure: responder proportion
HET_P_EVAL  <- 20              #   size of the evaluated group (20, 50 or 80)
HET_CELL    <- "Medium"        #   High, Medium, Low, Very Low, Extremely Low
HET_OTHER   <- 5e-4            #   effect of the other group (5e-4 or 6.25e-4)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
CAVEAT <- "Legacy version-2 output. Curves: mean of per-dataset ROC curves; AUC: mean per-dataset AUC. DiD scores contain ties at 0.05 (v2 scoring)."

# ---- helpers ------------------------------------------------------------------------
auc_mw <- function(score, truth) {                       # AUC within ONE dataset
  ok <- !is.na(score) & !is.na(truth); score <- score[ok]; truth <- truth[ok]
  n1 <- sum(truth == 1); n0 <- sum(truth == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
FPR_GRID <- seq(0, 1, by = 0.01)
roc_at_grid <- function(score, truth, grid = FPR_GRID) { # TPR of ONE dataset at fixed FPRs
  ok <- !is.na(score) & !is.na(truth); score <- score[ok]; truth <- truth[ok]
  n1 <- sum(truth == 1); n0 <- sum(truth == 0)
  if (n1 == 0 || n0 == 0) return(rep(NA_real_, length(grid)))
  o <- order(score, decreasing = TRUE); s <- score[o]; t <- truth[o]
  last <- c(s[-1] != s[-length(s)], TRUE)                # end of each block of tied scores
  tp <- c(0, cumsum(t == 1)[last] / n1)
  fp <- c(0, cumsum(t == 0)[last] / n0)
  tp[findInterval(grid, fp)]
}
mcse <- function(x) { x <- x[!is.na(x)]; if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x)) }
theme_mimosa <- function() {
  theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
          plot.caption = element_text(size = 7, color = "grey40", hjust = 0),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text = element_text(face = "bold", size = 9),
          strip.background = element_rect(fill = "grey92", color = NA),
          panel.grid.minor = element_blank(),
          legend.position = "bottom", legend.box = "horizontal")
}

# =============================================================================
# PART 1: CELL-COUNT IMBALANCE
# =============================================================================
load(IMB_FILE)                                   # -> results_continuous
imb <- results_continuous %>%
  mutate(Fraction = as.numeric(sub("_Count_.*", "", Cell_range)),
         Type = case_when(grepl("_Ref$", Cell_range) ~ "Reference",
                          grepl("_Bal_", Cell_range) ~ "Balanced",
                          grepl("_Imb_", Cell_range) ~ "Imbalanced"),
         Pattern = ifelse(Type == "Imbalanced", sub(".*_Imb_", "", Cell_range), NA_character_),
         k = case_when(Type == "Balanced" ~ suppressWarnings(as.integer(sub(".*_Bal_", "", Cell_range))),
                       Type == "Imbalanced" ~ lengths(strsplit(Pattern, "_")),
                       TRUE ~ 0L))
stopifnot(!anyNA(imb$Type), !anyNA(imb$Fraction))
message("Imbalance scenarios: ", n_distinct(imb$Cell_range), "; datasets: ", n_distinct(imb$Run_ID))

# ---- per-dataset AUC ----
imb_ds <- imb %>% group_by(Run_ID, Cell_range, Type, Pattern, k, Fraction) %>%
  summarise(MIMOSA2 = auc_mw(MIMOSA2_prob, Truth), `DiD Baseline` = auc_mw(DiD_GLM_prob, Truth),
            .groups = "drop") %>%
  pivot_longer(c(MIMOSA2, `DiD Baseline`), names_to = "Method", values_to = "AUC")
imb_auc <- imb_ds %>% group_by(Cell_range, Type, Pattern, k, Fraction, Method) %>%
  summarise(n = sum(!is.na(AUC)), AUC_se = mcse(AUC), AUC = mean(AUC, na.rm = TRUE), .groups = "drop")
write.csv(imb_auc, file.path(TAB, "legacy_imbalance_auc_by_scenario.csv"), row.names = FALSE)

level_lab <- function(f) factor(paste0(100 * f, "% target"),
                                levels = paste0(100 * sort(unique(imb$Fraction), decreasing = TRUE), "% target"))

# ---- ROC figure for one depletion pattern (your cell_imbalance_plot) ----
imbalance_roc_plot <- function(pattern) {
  kk <- length(strsplit(pattern, "_")[[1]])
  sub_ids <- imb %>%
    filter(Type == "Reference" | (Type == "Imbalanced" & Pattern == pattern) | (Type == "Balanced" & k == kk))
  roc <- sub_ids %>%
    group_by(Run_ID, Fraction, Type) %>%
    group_modify(~ bind_rows(
      data.frame(Method = "MIMOSA2",      FPR = FPR_GRID, TPR = roc_at_grid(.x$MIMOSA2_prob, .x$Truth)),
      data.frame(Method = "DiD Baseline", FPR = FPR_GRID, TPR = roc_at_grid(.x$DiD_GLM_prob, .x$Truth)))) %>%
    ungroup() %>%
    mutate(Condition = ifelse(Type == "Imbalanced", "Imbalanced", "Balanced")) %>%
    group_by(Fraction, Condition, Method, FPR) %>%
    summarise(TPR = mean(TPR, na.rm = TRUE), .groups = "drop") %>%
    mutate(Count_Level = level_lab(Fraction))
  labs_df <- imb_auc %>%
    filter(Type == "Reference" | (Type == "Imbalanced" & Pattern == pattern) | (Type == "Balanced" & k == kk)) %>%
    mutate(Condition = ifelse(Type == "Imbalanced", "Imbalanced", "Balanced"),
           Count_Level = level_lab(Fraction),
           label = sprintf("%s %s: %.3f", ifelse(Method == "MIMOSA2", "MIMOSA2", "DiD"), substr(Condition, 1, 3), AUC)) %>%
    group_by(Count_Level) %>% arrange(Condition, desc(Method), .by_group = TRUE) %>%
    mutate(x = 0.40, y = seq(0.30, 0.05, length.out = n())) %>% ungroup()
  ggplot(roc, aes(FPR, TPR, color = Condition, linetype = Method, group = interaction(Method, Condition))) +
    geom_line(linewidth = 1, alpha = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey10", linewidth = 0.5) +
    geom_text(data = labs_df, aes(x = x, y = y, label = label, color = Condition), inherit.aes = FALSE,
              size = 2.4, hjust = 0, show.legend = FALSE) +
    facet_wrap(~ Count_Level, nrow = 2, ncol = 4) +
    scale_color_manual(values = c("Balanced" = "orchid", "Imbalanced" = "darkgreen")) +
    scale_linetype_manual(values = c("MIMOSA2" = "solid", "DiD Baseline" = "dotted")) +
    guides(linetype = guide_legend(order = 1, title = "Model", override.aes = list(linewidth = 1.2)),
           color = guide_legend(order = 2, title = "Parent cell count", override.aes = list(linewidth = 1.2))) +
    labs(title = "MIMOSA2 performance under parent cell count imbalance.",
         subtitle = sprintf("ROC analysis. Imbalanced: %s depleted; balanced: same total cells spread over all four assays.",
                            gsub("_", " + ", pattern)),
         x = "1-Specificity", y = "Sensitivity", caption = CAVEAT) +
    theme_mimosa()
}
p <- imbalance_roc_plot(IMB_PATTERN); print(p)
ggsave(file.path(FIG, "cell_imbalance_plot.pdf"), plot = p, width = 8, height = 6)
for (pt in sort(unique(na.omit(imb$Pattern))))
  ggsave(file.path(FIG, sprintf("cell_imbalance_plot_%s.pdf", pt)), plot = imbalance_roc_plot(pt), width = 8, height = 6)

# ---- AUC odds ratio vs the 150k reference (your auc_or_plot), all patterns ----
ref <- imb_auc %>% filter(Type == "Reference") %>% select(Method, Ref_AUC = AUC, Ref_se = AUC_se)
or_df <- bind_rows(
  imb_auc %>% filter(Type == "Imbalanced"),
  # balanced control repeated for each pattern with the same number of depleted assays
  imb_auc %>% filter(Type == "Balanced") %>% select(-Pattern) %>%
    inner_join(imb_auc %>% filter(Type == "Imbalanced") %>% distinct(Pattern, k), by = "k")) %>%
  left_join(ref, by = "Method") %>%
  mutate(Balance = factor(ifelse(Type == "Imbalanced", "Imbalanced", "Balanced"), levels = c("Balanced", "Imbalanced")),
         Method = factor(Method, levels = c("MIMOSA2", "DiD Baseline")),
         AUC_OR = (AUC / (1 - AUC)) / (Ref_AUC / (1 - Ref_AUC)),
         SE_log_OR = sqrt((AUC_se / (AUC * (1 - AUC)))^2 + (Ref_se / (Ref_AUC * (1 - Ref_AUC)))^2),
         OR_Lower = exp(log(AUC_OR) - 1.96 * SE_log_OR), OR_Upper = exp(log(AUC_OR) + 1.96 * SE_log_OR),
         Count_Level = factor(paste0(100 * Fraction, "% target"),
                              levels = paste0(100 * sort(unique(Fraction), decreasing = TRUE), "% target")),
         x_numeric = as.numeric(Count_Level),
         x_jitter = x_numeric + case_when(Method == "MIMOSA2" & Balance == "Balanced" ~ -0.12,
                                          Method == "MIMOSA2" & Balance == "Imbalanced" ~ -0.04,
                                          Method == "DiD Baseline" & Balance == "Balanced" ~ 0.04,
                                          TRUE ~ 0.12),
         Pattern_lab = factor(paste("depleted:", gsub("_", " + ", Pattern)),
                              levels = paste("depleted:", gsub("_", " + ", c("s1", "u1", "s0", "u0", "s1_s0", "u1_u0", "s1_u1", "s0_u0",
                                                                         "s1_u0", "s0_u1", "s1_u1_u0", "s1_u1_s0", "s1_u0_s0", "u1_u0_s0")))))
write.csv(or_df, file.path(TAB, "legacy_imbalance_auc_or.csv"), row.names = FALSE)
lvl <- levels(or_df$Count_Level)
or_plot <- function(d, title_extra = "") {
  ggplot(d, aes(x_jitter, AUC_OR, group = interaction(Method, Balance), colour = Balance, linetype = Method)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey10", linewidth = 0.5) +
    geom_line(linewidth = 1, alpha = 0.85) +
    geom_errorbar(aes(ymin = OR_Lower, ymax = OR_Upper), width = 0, linetype = "solid", linewidth = 0.4) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = seq_along(lvl), labels = lvl) + scale_y_log10() +
    scale_colour_manual(values = c("Balanced" = "orchid", "Imbalanced" = "darkgreen")) +
    scale_linetype_manual(values = c("MIMOSA2" = "solid", "DiD Baseline" = "dashed")) +
    guides(linetype = guide_legend(order = 1, title = "Model", override.aes = list(linewidth = 1.2)),
           colour = guide_legend(order = 2, title = "Parent cell balance", override.aes = list(linewidth = 1.2))) +
    labs(title = "MIMOSA2 performance under parent cell count imbalance.",
         subtitle = paste0("AUC odds ratios relative to 100% target count reference.", title_extra),
         x = "Target downsampled parent cell count", y = "AUC odds ratio", caption = paste(CAVEAT, "Bars: 95% CI from Monte Carlo SEs.")) +
    theme_mimosa() + theme(legend.key.width = unit(1, "cm"))
}
p <- or_plot(filter(or_df, Pattern == IMB_PATTERN), sprintf(" Depleted: %s.", IMB_PATTERN)); print(p)
ggsave(file.path(FIG, "auc_or_plot.pdf"), plot = p, width = 8, height = 6)
ggsave(file.path(FIG, "auc_or_plot_all_patterns.pdf"),
       plot = or_plot(or_df) + facet_wrap(~ Pattern_lab, ncol = 4) + theme(axis.text.x = element_text(size = 6)),
       width = 13, height = 12)

# =============================================================================
# PART 2: EFFECT-SIZE HETEROGENEITY
# =============================================================================
load(HET_FILE)                                   # -> results_continuous
het <- results_continuous %>%
  rename(Eval_effect = Small_effect, Other_effect = Large_effect, P_eval = P_small, P_other = P_big)
message("Heterogeneity datasets: ", n_distinct(het$Run_ID),
        " (NB: the evaluated group has the LARGER effect in version 2)")
MODEL_COLS <- c("MIMOSA2_prob_independent" = "Independent", "MIMOSA2_prob_combined" = "Combined", "DiD_GLM_prob" = "DiD GLM")
het_long <- het %>%
  pivot_longer(names(MODEL_COLS), names_to = "Model", values_to = "Probability") %>%
  mutate(Model_label = factor(MODEL_COLS[Model], levels = MODEL_COLS))

# ---- per-dataset AUC for every scenario ----
het_ds <- het_long %>%
  group_by(Run_ID, Res_prop, P_eval, P_other, Eval_effect, Other_effect, Cell_range, Model_label) %>%
  summarise(AUC = auc_mw(Probability, Truth), .groups = "drop")
het_auc <- het_ds %>%
  group_by(Res_prop, P_eval, P_other, Eval_effect, Other_effect, Cell_range, Model_label) %>%
  summarise(n = sum(!is.na(AUC)), AUC_se = mcse(AUC), AUC = mean(AUC, na.rm = TRUE), .groups = "drop")
het_diff <- het_ds %>% filter(Model_label != "DiD GLM") %>%
  pivot_wider(names_from = Model_label, values_from = AUC) %>%
  mutate(dAUC = Combined - Independent) %>%
  group_by(Res_prop, P_eval, P_other, Eval_effect, Other_effect, Cell_range) %>%
  summarise(n = sum(!is.na(dAUC)), dAUC_se = mcse(dAUC), dAUC = mean(dAUC, na.rm = TRUE), .groups = "drop")
write.csv(het_auc, file.path(TAB, "legacy_heterogeneity_auc.csv"), row.names = FALSE)
write.csv(het_diff, file.path(TAB, "legacy_heterogeneity_dAUC_combined_minus_independent.csv"), row.names = FALSE)

HET_COLS <- c("Independent" = "deeppink", "Combined" = "steelblue3", "DiD GLM" = "orange")
HET_LTY  <- c("Independent" = "solid", "Combined" = "dashed", "DiD GLM" = "dotdash")
eff_lab <- function(e) factor(paste0("Effect: ", formatC(e, format = "e", digits = 1)),
                              levels = paste0("Effect: ", formatC(sort(unique(e)), format = "e", digits = 1)))

# ---- ROC panels by evaluated-group effect (your roc_plot_2x2) ----
sel <- het_long %>% filter(Res_prop == HET_RHO, P_eval == HET_P_EVAL, Cell_range == HET_CELL,
                           abs(Other_effect - HET_OTHER) < 1e-12)
stopifnot(nrow(sel) > 0)
roc_df <- sel %>% group_by(Run_ID, Eval_effect, Model_label) %>%
  group_modify(~ data.frame(FPR = FPR_GRID, TPR = roc_at_grid(.x$Probability, .x$Truth))) %>%
  ungroup() %>% group_by(Eval_effect, Model_label, FPR) %>%
  summarise(TPR = mean(TPR, na.rm = TRUE), .groups = "drop") %>%
  mutate(Effect_clean = eff_lab(Eval_effect))
auc_txt <- het_auc %>% filter(Res_prop == HET_RHO, P_eval == HET_P_EVAL, Cell_range == HET_CELL,
                              abs(Other_effect - HET_OTHER) < 1e-12) %>%
  mutate(Effect_clean = eff_lab(Eval_effect), label = sprintf("%s: %.3f", Model_label, AUC)) %>%
  group_by(Effect_clean) %>% arrange(Model_label, .by_group = TRUE) %>%
  mutate(x = 0.45, y = seq(0.25, 0.05, length.out = n())) %>% ungroup()
roc_plot <- ggplot(roc_df, aes(FPR, TPR, color = Model_label, linetype = Model_label)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  geom_text(data = auc_txt, aes(x = x, y = y, label = label, color = Model_label), inherit.aes = FALSE,
            size = 2.6, hjust = 0, show.legend = FALSE) +
  facet_wrap(~ Effect_clean, ncol = 3) +
  scale_color_manual(values = HET_COLS) + scale_linetype_manual(values = HET_LTY) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.01, 0.01)) +
  labs(title = "ROC curve comparison across effect sizes.",
       subtitle = sprintf("Independent vs combined MIMOSA2 fit vs DiD GLM. %s, %d evaluated : %d other subjects, %s cells, other-group effect %s.",
                          sub("Prop_", "rho = ", HET_RHO), HET_P_EVAL, 100 - HET_P_EVAL, HET_CELL, formatC(HET_OTHER, format = "e", digits = 2)),
       x = "1 - Specificity (False Positive Rate)", y = "Sensitivity (True Positive Rate)",
       color = "Model", linetype = "Model",
       caption = paste("Panel effect = effect of the EVALUATED group (v2 column 'Small_effect', which holds the larger effect).", CAVEAT)) +
  theme_mimosa() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))
print(roc_plot)
ggsave(file.path(FIG, "roc_2x2_3models.pdf"), plot = roc_plot, width = 9, height = 7)

# ---- AUC across all cell counts and group sizes ----
d <- het_auc %>% filter(Res_prop == HET_RHO, abs(Other_effect - HET_OTHER) < 1e-12) %>%
  mutate(Split = factor(sprintf("%d evaluated : %d other", P_eval, P_other),
                        levels = sprintf("%d evaluated : %d other", c(20, 50, 80), c(80, 50, 20))),
         Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low", "Very Low", "Extremely Low")),
         Effect_clean = eff_lab(Eval_effect))
p <- ggplot(d, aes(Effect_clean, AUC, color = Model_label, fill = Model_label, linetype = Model_label, group = Model_label)) +
  geom_ribbon(aes(ymin = AUC - 1.96 * AUC_se, ymax = AUC + 1.96 * AUC_se), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
  facet_grid(Split ~ Cell_range) +
  scale_color_manual(values = HET_COLS) + scale_fill_manual(values = HET_COLS) + scale_linetype_manual(values = HET_LTY) +
  labs(title = "Effect-size heterogeneity: discrimination of the evaluated group.",
       subtitle = sprintf("%s, other-group effect %s. Mean per-dataset AUC +/- 1.96 MCSE.", sub("Prop_", "rho = ", HET_RHO),
                          formatC(HET_OTHER, format = "e", digits = 2)),
       x = "Effect size of the evaluated group", y = "Mean per-dataset AUC", color = "Model", fill = "Model", linetype = "Model",
       caption = CAVEAT) +
  theme_mimosa() + theme(axis.text.x = element_text(size = 6))
print(p); ggsave(file.path(FIG, "het_auc_by_cell_and_split.pdf"), plot = p, width = 12, height = 8)

# ---- paired difference combined - independent ----
dd <- het_diff %>% filter(Res_prop == HET_RHO) %>%
  mutate(Split = factor(sprintf("%d evaluated : %d other", P_eval, P_other),
                        levels = sprintf("%d evaluated : %d other", c(20, 50, 80), c(80, 50, 20))),
         Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low", "Very Low", "Extremely Low")),
         Effect_clean = eff_lab(Eval_effect),
         Other = paste("other-group effect", formatC(Other_effect, format = "e", digits = 2)))
p <- ggplot(dd, aes(Effect_clean, dAUC, color = Other, group = Other)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = dAUC - 1.96 * dAUC_se, ymax = dAUC + 1.96 * dAUC_se),
                  position = position_dodge(width = 0.4), size = 0.2) +
  facet_grid(Split ~ Cell_range) +
  scale_color_manual(values = c("deeppink", "steelblue3"), name = NULL) +
  labs(title = "Cost of fitting MIMOSA2 to all subjects together.",
       subtitle = sprintf("Paired difference in AUC of the evaluated group: combined - independent fit. %s.", sub("Prop_", "rho = ", HET_RHO)),
       x = "Effect size of the evaluated group", y = "AUC(combined) - AUC(independent)",
       caption = paste("Below 0 = fitting together hurts the evaluated group. Bars: +/- 1.96 MCSE.", CAVEAT)) +
  theme_mimosa() + theme(axis.text.x = element_text(size = 6))
print(p); ggsave(file.path(FIG, "het_dAUC_combined_vs_independent.pdf"), plot = p, width = 12, height = 8)

message("Done. Figures in ", FIG, ", tables in ", TAB,
        ". Check: reference AUC MIMOSA2 ~0.940, DiD ~0.897 (Python check of the same files).")
