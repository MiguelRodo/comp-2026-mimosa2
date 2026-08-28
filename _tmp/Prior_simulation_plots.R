# MIMOSA2 Simulation plots
# Isabella Lethbridge and Tayyeb Abrahams 
# 23 August 2026

#-------------------------------------------------------------------------------
# Load required libraries: 
library(tidyverse)
library(ggplot2)
library(cowplot)
library(plotROC)
library(Hmisc) # Required for mean_cl_boot
library(dplyr)
library(tidyr)
library(knitr)

# Load data:
setwd("C:/Github/comp-2026-mimosa2")
load('_simulations/Simulation_3.0.Rdata')

#-------------------------------------------------------------------------------
# Sensitivity analysis of MIMOSA2
# TPR at 1% FDR 
# 'prior_plot' 
#-------------------------------------------------------------------------------
# Filter data:
plot_data = subset(results_summary, Status == 'Success')

# Fact cell count: 
plot_data$Cell_range = factor(
  plot_data$Cell_range,
  levels = c("High", "Medium", "Low", "Very Low")
)

# Format effect sizes:
unique_effects = sort(unique(plot_data$Effect), decreasing = TRUE)

plot_data$Effect_fact = factor(
  plot_data$Effect,
  levels = unique_effects,
  labels = formatC(unique_effects, format = "e", digits = 1)
)

# Relabel Distribution: 
plot_data$Distribution_clean = factor(
  plot_data$Distribution,
  levels = c("SX", "BB", "Beta", "EG", "LN"),
  labels = c("Simplex", "Bimodal Beta", "Beta", "Exponential Gamma", "Logit Normal")
)

# Plot: 
prior_plot = ggplot(data = plot_data,
                   mapping = aes(x = Effect_fact, 
                                 y = TPR_001, 
                                 group = Cell_range,
                                 color = Cell_range,
                                 fill = Cell_range)) +
  # Add 95% Bootstrap confidence interval band:
  stat_summary(
    fun.data = mean_cl_boot, 
    geom = "ribbon", 
    alpha = 0.25, 
    color = NA,
    fun.args = list(B = 2000, 
                    conf.int = 0.95)
  ) +
  # Add mean line:
  stat_summary(fun = mean, 
               geom = "line", 
               linewidth = 0.8) +
  # Add mean points:
  stat_summary(fun = mean, 
               geom = "point", 
               size = 2) +
  # Remove variable prefix on RHS strip labels using labeller: 
  facet_wrap(
    ~ Distribution_clean, 
    labeller = labeller(Distribution_clean = label_value)
  ) +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = c(0, 0.5, 1.0), 
                     expand = c(0.02, 0.02)) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "steelblue3", 
                                "Low" = "orange",
                                "Very Low" = "purple")) +
  scale_fill_manual(values = c("High" = "deeppink", 
                               "Medium" = "steelblue3", 
                               "Low" = "orange",
                               "Very Low" = "purple")) +
  labs(
    title = "Sensitivity analysis of MIMOSA2.",
    subtitle = "True positive rate (TPR) at 1% nominal false discovery rate (FDR) threshold.",
    x = "Effect size",
    y = "TPR",
    color = "Cell count",
    fill = "Cell count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(prior_plot)
# ggsave('_fig/prior_plot.pdf',plot=prior_plot,width=8,height=6)

#-------------------------------------------------------------------------------
# Simulation performance of MIMOSA2
# ROC analysis 
# 'prior_ROC_plot' 
#-------------------------------------------------------------------------------
# Prepare data: 
ROC_data_prepared = results_continuous |>
  filter(!is.na(MIMOSA2_prob), !is.na(DiD_GLM_prob)) |>
  pivot_longer(
    cols = c(MIMOSA2_prob, DiD_GLM_prob),
    names_to = "Method",
    values_to = "Score"
  ) |>
  mutate(
    Method = ifelse(Method == "MIMOSA2_prob", "MIMOSA2", "DiD Baseline"),
    Sample_Size = paste0("N: ", P),
    Effect_Label = paste0("Effect: ", Effect),
    # Fixed factor definition and relabeling for Cell_range
    Cell_range = factor(
      Cell_range,
      levels = c("High", "Medium", "Low", "Very Low"),
      labels = c("High", "Medium", "Low", "Very Low")
    ),
    # Fixed factor definition for Distribution
    Distribution_clean = factor(
      Distribution,
      levels = c("SX", "BB", "Beta", "EG", "LN"),
      labels = c("Simplex", "Bimodal Beta", "Beta", "Exponential Gamma", "Logit Normal")
    )
  )

# Plot: 
prior_ROC_plot = ggplot(
  data = ROC_data_prepared,
  mapping = aes(
    d = Truth, 
    m = Score, 
    color = Cell_range, 
    linetype = Method, 
    group = interaction(Method, Cell_range)
  )
) +
  geom_roc(n.cuts = 0, size = 0.8, linealpha = 0.9) +
  geom_abline(
    slope = 1, 
    intercept = 0, 
    linetype = "dashed", 
    color = "grey10", 
    linewidth = 0.5
  ) +
  facet_grid(
    Distribution_clean ~ Effect, 
    labeller = labeller(
      Distribution_clean = label_value, 
      Effect = label_both
    )
  ) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "steelblue3", 
                                "Low" = "orange",
                                "Very Low" = "purple")) +
  scale_linetype_manual(values = c("MIMOSA2" = "solid", 
                                   "DiD Baseline" = "dotted")) +
  guides(
    color = guide_legend(
      order = 2, 
      title = "Cell Count",
      override.aes = list(linewidth = 1.2, size = 1.5)
    ),
    linetype = guide_legend(
      order = 1, 
      title = "Model",
      override.aes = list(linewidth = 1.2, size = 2)
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Simulation performance of MIMOSA2.",
    subtitle = "ROC analysis.",
    x = "FPR",
    y = "TPR"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

print(prior_ROC_plot)
# ggsave('_fig/prior_ROC_plot.pdf',plot=prior_ROC_plot,width=10,height=10)

#-------------------------------------------------------------------------------
# Simulation performance of MIMOSA2
# ROC analysis 
# 'prior_ROC_comp_plot' 
#-------------------------------------------------------------------------------
# Clean data:
df_clean <- results_continuous |>
  filter(!is.na(MIMOSA2_prob), !is.na(DiD_GLM_prob)) |>
  pivot_longer(c(MIMOSA2_prob, DiD_GLM_prob), names_to = "Method", values_to = "Score") |>
  mutate(Method = ifelse(Method == "MIMOSA2_prob", "MIMOSA2", "DiD Baseline"))

# Relabel 'Distribution':
comparison_labels <- c(
  "SX" = "Beta vs. Simplex",
  "BB" = "Beta vs. Bimodal Beta",
  "EG" = "Beta vs. Exponential Gamma",
  "LN" = "Beta vs. Logit Normal"
)

df_paired <- bind_rows(lapply(names(comparison_labels), function(d) {
  df_clean |>
    filter(Distribution %in% c("Beta", d)) |>
    mutate(
      Comparison = factor(comparison_labels[d], levels = unname(comparison_labels)),
      Distribution_clean = factor(
        Distribution,
        levels = c("Beta", "SX", "BB", "EG", "LN"),
        labels = c("Beta", "Simplex", "Bimodal Beta", "Exponential Gamma", "Logit Normal")
      )
    )
}))

# Plot:
prior_ROC_comp_plot = ggplot(
  data = df_paired,
  mapping = aes(
    d = Truth, 
    m = Score, 
    color = Distribution_clean, 
    linetype = Method
  )
) +
  geom_roc(n.cuts = 0, size = 1.2, linealpha = 0.7) +
  scale_alpha_identity() + 
  geom_abline(
    slope = 1, 
    intercept = 0, 
    linetype = "dashed", 
    color = "grey10", 
    linewidth = 0.5
  ) +
  facet_wrap(
    ~ Comparison, 
    ncol = 2,
    labeller = labeller(Comparison = label_value)
  ) +
  scale_color_manual(values = c(
    "Beta"              = "grey10", 
    "Simplex"           = "deeppink", 
    "Bimodal Beta"      = "steelblue3", 
    "Exponential Gamma" = "orange", 
    "Logit Normal"      = "purple"
  )) +
  scale_linetype_manual(values = c(
    "MIMOSA2"      = "solid", 
    "DiD Baseline" = "dotted"
  )) +
  guides(
    color = guide_legend(
      order = 1, 
      title = "Distribution",
      override.aes = list(linewidth = 1.2, alpha = 1)
    ),
    linetype = guide_legend(
      order = 2, 
      title = "Model",
      override.aes = list(linewidth = 1.2, alpha = 1)
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Simulation performance of MIMOSA2.",
    subtitle = "ROC analysis.",
    x = "FPR",
    y = "TPR"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical"
  )

print(prior_ROC_comp_plot)
# ggsave('_fig/prior_ROC_comp_plot.pdf',plot=prior_ROC_comp_plot,width=8,height=6)

#-------------------------------------------------------------------------------
# AUROC Calculation
#-------------------------------------------------------------------------------
# Calculate AUROC per scenario:
AUROC = ROC_data_prepared |>
  group_by(Distribution, P, Cell_range, Effect, Method) |>
  filter(length(unique(Truth)) == 2) |>
  do(plotROC::calc_auc(
    ggplot(., aes(d = Truth, m = Score)) + geom_roc()
  )) |>
  ungroup() |>
  rename(AUROC = AUC)

write.table(AUROC, "auc_prior.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Summarise AUROC matrix across Distribution, Effect, and Cell Range:
AUROC_matrix = AUROC |>
  mutate(
    Distribution_clean = factor(
      Distribution,
      levels = c("Beta", "SX", "BB", "EG", "LN"),
      labels = c("Beta", "Simplex", "Bimodal Beta", "Exponential Gamma", "Logit Normal")
    ),
    Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low", "Very Low")),
    Effect_clean = paste0(formatC(Effect, format = "e", digits = 1))
  ) |>
  group_by(Distribution_clean, Effect_clean, Cell_range, Method) |>
  summarise(
    Mean_AUROC = sprintf("%.3f", mean(AUROC, na.rm = TRUE)), 
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = Method,
    values_from = Mean_AUROC
  ) |>
  arrange(Distribution_clean, Effect_clean, Cell_range)

# Save as file: 
write.table(
  AUROC_matrix, 
  "auc_prior_summary_matrix.txt", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

# Table: 
kable(
  AUROC_matrix, 
  caption = 'AUROC values across Prior Distributions (50% Responders)',
  col.names = c("Prior Distribution", "Effect Size", "Cell Count", "MIMOSA2", "DiD Baseline")
)

# INCORRECT
# Table: AUROC values across Prior Distributions (50% Responders)
# 
# |Prior Distribution |Effect Size |Cell Count |MIMOSA2 |DiD Baseline |
# |:------------------|:-----------|:----------|:-------|:------------|
# |Beta               |1.0e-03     |High       |0.956   |0.986        |
# |Beta               |1.0e-03     |Medium     |0.912   |0.963        |
# |Beta               |1.0e-03     |Low        |0.848   |0.933        |
# |Beta               |1.0e-03     |Very Low   |0.762   |0.881        |
# |Beta               |1.2e-04     |High       |0.903   |0.914        |
# |Beta               |1.2e-04     |Medium     |0.824   |0.857        |
# |Beta               |1.2e-04     |Low        |0.725   |0.746        |
# |Beta               |1.2e-04     |Very Low   |0.660   |0.692        |
# |Beta               |2.5e-04     |High       |0.896   |0.927        |
# |Beta               |2.5e-04     |Medium     |0.855   |0.908        |
# |Beta               |2.5e-04     |Low        |0.742   |0.785        |
# |Beta               |2.5e-04     |Very Low   |0.676   |0.732        |
# |Beta               |6.2e-05     |High       |0.882   |0.895        |
# |Beta               |6.2e-05     |Medium     |0.842   |0.882        |
# |Beta               |6.2e-05     |Low        |0.657   |0.699        |
# |Beta               |6.2e-05     |Very Low   |0.614   |0.647        |
# |Simplex            |1.0e-03     |High       |0.938   |0.967        |
# |Simplex            |1.0e-03     |Medium     |0.923   |0.968        |
# |Simplex            |1.0e-03     |Low        |0.860   |0.941        |
# |Simplex            |1.0e-03     |Very Low   |0.805   |0.912        |
# |Simplex            |1.2e-04     |High       |0.876   |0.896        |
# |Simplex            |1.2e-04     |Medium     |0.805   |0.843        |
# |Simplex            |1.2e-04     |Low        |0.661   |0.678        |
# |Simplex            |1.2e-04     |Very Low   |0.630   |0.675        |
# |Simplex            |2.5e-04     |High       |0.901   |0.939        |
# |Simplex            |2.5e-04     |Medium     |0.865   |0.898        |
# |Simplex            |2.5e-04     |Low        |0.735   |0.781        |
# |Simplex            |2.5e-04     |Very Low   |0.623   |0.707        |
# |Simplex            |6.2e-05     |High       |0.836   |0.855        |
# |Simplex            |6.2e-05     |Medium     |0.818   |0.828        |
# |Simplex            |6.2e-05     |Low        |0.627   |0.653        |
# |Simplex            |6.2e-05     |Very Low   |0.586   |0.626        |
# |Bimodal Beta       |1.0e-03     |High       |0.932   |0.963        |
# |Bimodal Beta       |1.0e-03     |Medium     |0.925   |0.972        |
# |Bimodal Beta       |1.0e-03     |Low        |0.843   |0.925        |
# |Bimodal Beta       |1.0e-03     |Very Low   |0.780   |0.900        |
# |Bimodal Beta       |1.2e-04     |High       |0.905   |0.934        |
# |Bimodal Beta       |1.2e-04     |Medium     |0.822   |0.865        |
# |Bimodal Beta       |1.2e-04     |Low        |0.691   |0.726        |
# |Bimodal Beta       |1.2e-04     |Very Low   |0.646   |0.693        |
# |Bimodal Beta       |2.5e-04     |High       |0.892   |0.931        |
# |Bimodal Beta       |2.5e-04     |Medium     |0.856   |0.903        |
# |Bimodal Beta       |2.5e-04     |Low        |0.702   |0.779        |
# |Bimodal Beta       |2.5e-04     |Very Low   |0.659   |0.728        |
# |Bimodal Beta       |6.2e-05     |High       |0.891   |0.901        |
# |Bimodal Beta       |6.2e-05     |Medium     |0.813   |0.848        |
# |Bimodal Beta       |6.2e-05     |Low        |0.689   |0.732        |
# |Bimodal Beta       |6.2e-05     |Very Low   |0.615   |0.645        |
# |Exponential Gamma  |1.0e-03     |High       |0.927   |0.959        |
# |Exponential Gamma  |1.0e-03     |Medium     |0.905   |0.938        |
# |Exponential Gamma  |1.0e-03     |Low        |0.826   |0.915        |
# |Exponential Gamma  |1.0e-03     |Very Low   |0.832   |0.903        |
# |Exponential Gamma  |1.2e-04     |High       |0.870   |0.885        |
# |Exponential Gamma  |1.2e-04     |Medium     |0.864   |0.882        |
# |Exponential Gamma  |1.2e-04     |Low        |0.678   |0.716        |
# |Exponential Gamma  |1.2e-04     |Very Low   |0.620   |0.675        |
# |Exponential Gamma  |2.5e-04     |High       |0.913   |0.924        |
# |Exponential Gamma  |2.5e-04     |Medium     |0.885   |0.910        |
# |Exponential Gamma  |2.5e-04     |Low        |0.757   |0.815        |
# |Exponential Gamma  |2.5e-04     |Very Low   |0.624   |0.674        |
# |Exponential Gamma  |6.2e-05     |High       |0.859   |0.880        |
# |Exponential Gamma  |6.2e-05     |Medium     |0.822   |0.828        |
# |Exponential Gamma  |6.2e-05     |Low        |0.699   |0.747        |
# |Exponential Gamma  |6.2e-05     |Very Low   |0.607   |0.665        |
# |Logit Normal       |1.0e-03     |High       |0.937   |0.967        |
# |Logit Normal       |1.0e-03     |Medium     |0.921   |0.955        |
# |Logit Normal       |1.0e-03     |Low        |0.816   |0.881        |
# |Logit Normal       |1.0e-03     |Very Low   |0.818   |0.908        |
# |Logit Normal       |1.2e-04     |High       |0.890   |0.908        |
# |Logit Normal       |1.2e-04     |Medium     |0.826   |0.844        |
# |Logit Normal       |1.2e-04     |Low        |0.686   |0.755        |
# |Logit Normal       |1.2e-04     |Very Low   |0.659   |0.760        |
# |Logit Normal       |2.5e-04     |High       |0.883   |0.918        |
# |Logit Normal       |2.5e-04     |Medium     |0.876   |0.899        |
# |Logit Normal       |2.5e-04     |Low        |0.734   |0.779        |
# |Logit Normal       |2.5e-04     |Very Low   |0.646   |0.756        |
# |Logit Normal       |6.2e-05     |High       |0.903   |0.908        |
# |Logit Normal       |6.2e-05     |Medium     |0.812   |0.840        |
# |Logit Normal       |6.2e-05     |Low        |0.685   |0.727        |
# |Logit Normal       |6.2e-05     |Very Low   |0.627   |0.705        |

AUROC_matrix_2 <- AUROC |>
  mutate(
    Distribution_clean = factor(
      Distribution,
      levels = c("Beta", "SX", "BB", "EG", "LN"),
      labels = c("Beta", "Simplex", "Bimodal Beta", "Exponential Gamma", "Logit Normal")
    ),
    Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low", "Very Low")),
    Effect_clean = formatC(Effect, format = "e", digits = 1),
    Method = factor(Method, levels = c("MIMOSA2", "DiD Baseline")) # Ensures MIMOSA2 is first
  ) |>
  group_by(Distribution_clean, Effect_clean, Cell_range, Method) |>
  summarise(
    Mean_AUROC = sprintf("%.3f", mean(AUROC, na.rm = TRUE)), 
    .groups = "drop"
  ) |>
  arrange(Distribution_clean, Effect_clean, Cell_range, Method) |>
  # Combine MIMOSA2 and DiD Baseline into "val1, val2"
  group_by(Distribution_clean, Effect_clean, Cell_range) |>
  summarise(
    Combined_AUROC = paste(Mean_AUROC, collapse = ", "),
    .groups = "drop"
  ) |>
  # Pivot Cell_range to columns
  pivot_wider(
    names_from = Cell_range,
    values_from = Combined_AUROC
  ) |>
  arrange(Distribution_clean, Effect_clean)

# Save as file: 
write.table(
  AUROC_matrix_2, 
  "auc_prior_summary_matrix_2.txt", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

# Render Table: 
kable(
  AUROC_matrix_2, 
  caption = 'AUROC values (MIMOSA2, DiD Baseline) across Prior Distributions (50% Responders)',
  col.names = c("Prior Distribution", "Effect Size", "High", "Medium", "Low", "Very Low")
)

# Table: AUROC values (MIMOSA2, DiD Baseline) across Prior Distributions (50% Responders)
# 
# |Prior Distribution |Effect Size |High         |Medium       |Low          |Very Low     |
# |:------------------|:-----------|:------------|:------------|:------------|:------------|
# |Beta               |1.0e-03     |0.986, 0.956 |0.963, 0.912 |0.933, 0.848 |0.881, 0.762 |
# |Beta               |1.2e-04     |0.914, 0.903 |0.857, 0.824 |0.746, 0.725 |0.692, 0.660 |
# |Beta               |2.5e-04     |0.927, 0.896 |0.908, 0.855 |0.785, 0.742 |0.732, 0.676 |
# |Beta               |6.2e-05     |0.895, 0.882 |0.882, 0.842 |0.699, 0.657 |0.647, 0.614 |
# |Simplex            |1.0e-03     |0.967, 0.938 |0.968, 0.923 |0.941, 0.860 |0.912, 0.805 |
# |Simplex            |1.2e-04     |0.896, 0.876 |0.843, 0.805 |0.678, 0.661 |0.675, 0.630 |
# |Simplex            |2.5e-04     |0.939, 0.901 |0.898, 0.865 |0.781, 0.735 |0.707, 0.623 |
# |Simplex            |6.2e-05     |0.855, 0.836 |0.828, 0.818 |0.653, 0.627 |0.626, 0.586 |
# |Bimodal Beta       |1.0e-03     |0.963, 0.932 |0.972, 0.925 |0.925, 0.843 |0.900, 0.780 |
# |Bimodal Beta       |1.2e-04     |0.934, 0.905 |0.865, 0.822 |0.726, 0.691 |0.693, 0.646 |
# |Bimodal Beta       |2.5e-04     |0.931, 0.892 |0.903, 0.856 |0.779, 0.702 |0.728, 0.659 |
# |Bimodal Beta       |6.2e-05     |0.901, 0.891 |0.848, 0.813 |0.732, 0.689 |0.645, 0.615 |
# |Exponential Gamma  |1.0e-03     |0.959, 0.927 |0.938, 0.905 |0.915, 0.826 |0.903, 0.832 |
# |Exponential Gamma  |1.2e-04     |0.885, 0.870 |0.882, 0.864 |0.716, 0.678 |0.675, 0.620 |
# |Exponential Gamma  |2.5e-04     |0.924, 0.913 |0.910, 0.885 |0.815, 0.757 |0.674, 0.624 |
# |Exponential Gamma  |6.2e-05     |0.880, 0.859 |0.828, 0.822 |0.747, 0.699 |0.665, 0.607 |
# |Logit Normal       |1.0e-03     |0.967, 0.937 |0.955, 0.921 |0.881, 0.816 |0.908, 0.818 |
# |Logit Normal       |1.2e-04     |0.908, 0.890 |0.844, 0.826 |0.755, 0.686 |0.760, 0.659 |
# |Logit Normal       |2.5e-04     |0.918, 0.883 |0.899, 0.876 |0.779, 0.734 |0.756, 0.646 |
# |Logit Normal       |6.2e-05     |0.908, 0.903 |0.840, 0.812 |0.727, 0.685 |0.705, 0.627 |

# Simulation count table:
Sim_count_matrix <- ROC_data_prepared |>
  mutate(
    Distribution_clean = factor(
      Distribution,
      levels = c("Beta", "SX", "BB", "EG", "LN"),
      labels = c("Beta", "Simplex", "Bimodal Beta", "Exponential Gamma", "Logit Normal")
    ),
    Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low", "Very Low")),
    Effect_clean = formatC(Effect, format = "e", digits = 1),
    Method = factor(Method, levels = c("MIMOSA2", "DiD Baseline")) # Ensures MIMOSA2 is first
  ) |>
  group_by(Distribution_clean, Effect_clean, Cell_range, Method) |>
  summarise(
    N_sims = n(), 
    .groups = "drop"
  ) |>
  arrange(Distribution_clean, Effect_clean, Cell_range, Method) |>
  # Combine MIMOSA2 and DiD Baseline simulation counts
  group_by(Distribution_clean, Effect_clean, Cell_range) |>
  summarise(
    Combined_Count = paste(N_sims, collapse = ", "),
    .groups = "drop"
  ) |>
  # Pivot Cell_range to columns
  pivot_wider(
    names_from = Cell_range,
    values_from = Combined_Count
  ) |>
  arrange(Distribution_clean, Effect_clean)

# Save as file:
write.table(
  Sim_count_matrix, 
  "simulation_counts_prior_matrix.txt", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

# Render Table:
kable(
  Sim_count_matrix, 
  caption = 'Simulation Counts (MIMOSA2, DiD Baseline) across Prior Distributions (50% Responders)',
  col.names = c("Prior Distribution", "Effect Size", "High", "Medium", "Low", "Very Low")
)

# Table: Simulation Counts (MIMOSA2, DiD Baseline) across Prior Distributions (50% Responders)
# 
# |Prior Distribution |Effect Size |High     |Medium   |Low      |Very Low |
# |:------------------|:-----------|:--------|:--------|:--------|:--------|
# |Beta               |1.0e-03     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Beta               |1.2e-04     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Beta               |2.5e-04     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Beta               |6.2e-05     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Simplex            |1.0e-03     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Simplex            |1.2e-04     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Simplex            |2.5e-04     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Simplex            |6.2e-05     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Bimodal Beta       |1.0e-03     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Bimodal Beta       |1.2e-04     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Bimodal Beta       |2.5e-04     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Bimodal Beta       |6.2e-05     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Exponential Gamma  |1.0e-03     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Exponential Gamma  |1.2e-04     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Exponential Gamma  |2.5e-04     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Exponential Gamma  |6.2e-05     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Logit Normal       |1.0e-03     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Logit Normal       |1.2e-04     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Logit Normal       |2.5e-04     |500, 500 |500, 500 |500, 500 |500, 500 |
# |Logit Normal       |6.2e-05     |500, 500 |500, 500 |500, 500 |500, 500 |

get_auc_diff_matrix <- function(auc_matrix, sim_matrix, r = 0.5, prop_pos = 0.5, digits = 3) {
  
  # Helper to compute Hanley & McNeil SE for a single AUC and sample size N
  calc_se <- function(auc, n) {
    n1 <- round(n * prop_pos)
    n0 <- n - n1
    q1 <- auc / (2 - auc)
    q2 <- (2 * (auc^2)) / (1 + auc)
    var_auc <- (auc * (1 - auc) + (n1 - 1) * (q1 - auc^2) + (n0 - 1) * (q2 - auc^2)) / (n1 * n0)
    return(sqrt(var_auc))
  }
  
  # Helper to process a single cell string pair
  process_cell <- function(auc_str, sim_str) {
    if (is.na(auc_str) || is.na(sim_str)) return(NA_character_)
    
    # Parse numbers from comma-separated strings
    aucs <- as.numeric(trimws(unlist(strsplit(as.character(auc_str), ","))))
    sims <- as.numeric(trimws(unlist(strsplit(as.character(sim_str), ","))))
    
    auc1 <- aucs[1] # MIMOSA2
    auc2 <- aucs[2] # DiD Baseline
    n1   <- sims[1]
    n2   <- sims[2]
    
    # Compute SEs
    se1 <- calc_se(auc1, n1)
    se2 <- calc_se(auc2, n2)
    
    # Paired standard error of difference
    se_diff <- sqrt(se1^2 + se2^2 - 2 * r * se1 * se2)
    diff <- auc1 - auc2
    
    ci_lower <- diff - 1.96 * se_diff
    ci_upper <- diff + 1.96 * se_diff
    
    # Format output as "Diff [Lower, Upper]"
    fmt <- paste0("%.", digits, "f")
    return(sprintf(paste0(fmt, " [", fmt, ", ", fmt, "]"), diff, ci_lower, ci_upper))
  }
  
  # Initialize output structure
  diff_matrix <- auc_matrix
  # Explicitly include "Very Low" in the value columns
  value_cols <- c("High", "Medium", "Low", "Very Low")
  
  # Loop over numeric columns
  for (col in value_cols) {
    diff_matrix[[col]] <- mapply(
      process_cell, 
      auc_matrix[[col]], 
      sim_matrix[[col]]
    )
  }
  
  return(diff_matrix)
}

# Generate the matrix:
AUC_Diff_matrix <- get_auc_diff_matrix(
  auc_matrix = AUROC_matrix_2, 
  sim_matrix = Sim_count_matrix, 
  r = 0.5
)

# Render formatted table:
kable(
  AUC_Diff_matrix, 
  caption = 'Δ AUROC (MIMOSA2 - DiD Baseline) with 95% Confidence Intervals across Prior Distributions',
  col.names = c("Prior Distribution", "Effect Size", "High", "Medium", "Low", "Very Low")
)

# Table: Δ AUROC (MIMOSA2 - DiD Baseline) with 95% Confidence Intervals across Prior Distributions
# 
# |Prior Distribution |Effect Size |High                  |Medium                |Low                   |Very Low              |
# |:------------------|:-----------|:---------------------|:---------------------|:---------------------|:---------------------|
# |Beta               |1.0e-03     |0.030 [0.014, 0.046]  |0.051 [0.028, 0.074]  |0.085 [0.055, 0.115]  |0.119 [0.082, 0.156]  |
# |Beta               |1.2e-04     |0.011 [-0.016, 0.038] |0.033 [-0.002, 0.068] |0.021 [-0.023, 0.065] |0.032 [-0.015, 0.079] |
# |Beta               |2.5e-04     |0.031 [0.004, 0.058]  |0.053 [0.022, 0.084]  |0.043 [0.001, 0.085]  |0.056 [0.011, 0.101]  |
# |Beta               |6.2e-05     |0.013 [-0.017, 0.043] |0.040 [0.007, 0.073]  |0.042 [-0.005, 0.089] |0.033 [-0.016, 0.082] |
# |Simplex            |1.0e-03     |0.029 [0.009, 0.049]  |0.045 [0.023, 0.067]  |0.081 [0.052, 0.110]  |0.107 [0.073, 0.141]  |
# |Simplex            |1.2e-04     |0.020 [-0.010, 0.050] |0.038 [0.001, 0.075]  |0.017 [-0.030, 0.064] |0.045 [-0.003, 0.093] |
# |Simplex            |2.5e-04     |0.038 [0.013, 0.063]  |0.033 [0.002, 0.064]  |0.046 [0.004, 0.088]  |0.084 [0.037, 0.131]  |
# |Simplex            |6.2e-05     |0.019 [-0.016, 0.054] |0.010 [-0.027, 0.047] |0.026 [-0.022, 0.074] |0.040 [-0.009, 0.089] |
# |Bimodal Beta       |1.0e-03     |0.031 [0.010, 0.052]  |0.047 [0.026, 0.068]  |0.082 [0.051, 0.113]  |0.120 [0.084, 0.156]  |
# |Bimodal Beta       |1.2e-04     |0.029 [0.004, 0.054]  |0.043 [0.008, 0.078]  |0.035 [-0.010, 0.080] |0.047 [0.000, 0.094]  |
# |Bimodal Beta       |2.5e-04     |0.039 [0.012, 0.066]  |0.047 [0.016, 0.078]  |0.077 [0.034, 0.120]  |0.069 [0.023, 0.115]  |
# |Bimodal Beta       |6.2e-05     |0.010 [-0.019, 0.039] |0.035 [-0.001, 0.071] |0.043 [-0.002, 0.088] |0.030 [-0.019, 0.079] |
# |Exponential Gamma  |1.0e-03     |0.032 [0.010, 0.054]  |0.033 [0.008, 0.058]  |0.089 [0.057, 0.121]  |0.071 [0.038, 0.104]  |
# |Exponential Gamma  |1.2e-04     |0.015 [-0.016, 0.046] |0.018 [-0.013, 0.049] |0.038 [-0.008, 0.084] |0.055 [0.007, 0.103]  |
# |Exponential Gamma  |2.5e-04     |0.011 [-0.014, 0.036] |0.025 [-0.003, 0.053] |0.058 [0.018, 0.098]  |0.050 [0.002, 0.098]  |
# |Exponential Gamma  |6.2e-05     |0.021 [-0.011, 0.053] |0.006 [-0.031, 0.043] |0.048 [0.004, 0.092]  |0.058 [0.010, 0.106]  |
# |Logit Normal       |1.0e-03     |0.030 [0.010, 0.050]  |0.034 [0.012, 0.056]  |0.065 [0.031, 0.099]  |0.090 [0.057, 0.123]  |
# |Logit Normal       |1.2e-04     |0.018 [-0.010, 0.046] |0.018 [-0.018, 0.054] |0.069 [0.025, 0.113]  |0.101 [0.056, 0.146]  |
# |Logit Normal       |2.5e-04     |0.035 [0.007, 0.063]  |0.023 [-0.007, 0.053] |0.045 [0.003, 0.087]  |0.110 [0.065, 0.155]  |
# |Logit Normal       |6.2e-05     |0.005 [-0.022, 0.032] |0.028 [-0.008, 0.064] |0.042 [-0.003, 0.087] |0.078 [0.031, 0.125]  |

# -----------------------------------------------------------------------------
# Simulation performance of MIMOSA2
# ROC analysis
# 'ROC_plot_prior_overlay' 
# -----------------------------------------------------------------------------
# Filter data: 
ROC_data_filtered <- ROC_data_prepared |>
  mutate(
    Distribution_clean = factor(
      Distribution,
      levels = c("SX", "BB", "Beta", "EG", "LN"),
      labels = c("Simplex", "Bimodal Beta", "Beta", "Exponential Gamma", "Logit Normal")
    ),
    Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low", "Very Low")),
    Effect_clean = paste0("Effect: ", formatC(Effect, format = "e", digits = 1))
  )

# Calculate AUC: 
auroc_df = ROC_data_filtered |>
  group_by(Distribution_clean, Effect_clean, Cell_range, Method) |>
  filter(length(unique(Truth)) == 2) |>
  do(plotROC::calc_auc(
    ggplot(., aes(d = Truth, m = Score)) + geom_roc()
  )) |>
  ungroup() |>
  rename(AUROC = AUC)

# Format AUC text overlay (Adjusted for 4 cell ranges x 2 methods = 8 text lines per panel):
auc_text_df = auroc_df |>
  mutate(
    Method_short = ifelse(Method == "DiD Baseline", "DiD", Method),
    auc_label = paste0(Method_short, ": ", sprintf("%.3f", AUROC))
  ) |>
  group_by(Distribution_clean, Effect_clean) |>
  arrange(desc(Method), desc(Cell_range), .by_group = TRUE) |>
  mutate(
    x = 0.50, 
    # Spread 8 labels cleanly down the panel:
    y = seq(0.40, 0.02, length.out = n()) 
  ) |>
  ungroup()

# Plot: 
ROC_plot_prior_overlay = ggplot(
  data = ROC_data_filtered,
  mapping = aes(
    d = Truth, 
    m = Score, 
    color = Cell_range, 
    linetype = Method, 
    group = interaction(Method, Cell_range)
  )
) +
  geom_roc(n.cuts = 0, size = 0.8, linealpha = 0.9) +
  geom_abline(
    slope = 1, 
    intercept = 0, 
    linetype = "dashed", 
    color = "grey10", 
    linewidth = 0.5
  ) +
  geom_text(
    data = auc_text_df,
    aes(x = x, y = y, label = auc_label, color = Cell_range),
    inherit.aes = FALSE,
    size = 2.8,
    hjust = 0,
    vjust = 1
  ) +
  facet_grid(Distribution_clean ~ Effect_clean) +
  scale_color_manual(values = c(
    "High"     = "deeppink", 
    "Medium"   = "steelblue3", 
    "Low"      = "orange",
    "Very Low" = "purple"
  )) +
  scale_linetype_manual(values = c(
    "MIMOSA2"      = "solid", 
    "DiD Baseline" = "dotted"
  )) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = sprintf("%.2f", seq(0, 1, 0.25))) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = sprintf("%.2f", seq(0, 1, 0.25))) +
  guides(
    color = guide_legend(
      order = 2, 
      title = "Cell Count",
      override.aes = list(linewidth = 1.2, size = 1.5, alpha = 1)
    ),
    linetype = guide_legend(
      order = 1, 
      title = "Model",
      override.aes = list(linewidth = 1.2, size = 2)
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Simulation performance of MIMOSA2.",
    subtitle = "ROC analysis.",
    x = "FPR",
    y = "TPR"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8.5),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

print(ROC_plot_prior_overlay) 
# ggsave('_fig/ROC_plot_prior_overlay.pdf', plot = ROC_plot_prior_overlay, width = 12, height = 12)

#------------------------------------------------------------------------------
# Simulation performance of MIMOSA2
# ROC analysis across prior distributions.
# 'ROC_plot_overlay' 
# Effect: small,med,large
# Res_prop: 10%,50%,90%
#------------------------------------------------------------------------------
# Filter data:
df_clean <- results_continuous |>
  filter(!is.na(MIMOSA2_prob), !is.na(DiD_GLM_prob)) |>
  pivot_longer(c(MIMOSA2_prob, DiD_GLM_prob), names_to = "Method", values_to = "Score") |>
  mutate(Method = ifelse(Method == "MIMOSA2_prob", "MIMOSA2", "DiD Baseline"))

# Relabel 'DIstribution':
comparison_labels <- c(
  "SX" = "Beta vs. Simplex",
  "BB" = "Beta vs. Bimodal Beta",
  "EG" = "Beta vs. Exponential Gamma",
  "LN" = "Beta vs. Logit Normal"
)

df_paired <- bind_rows(lapply(names(comparison_labels), function(d) {
  df_clean |>
    filter(Distribution %in% c("Beta", d)) |>
    mutate(
      Comparison = factor(comparison_labels[d], levels = unname(comparison_labels)),
      Distribution_clean = factor(
        Distribution,
        levels = c("Beta", "SX", "BB", "EG", "LN"),
        labels = c("Beta", "Simplex", "Bimodal Beta", "Exponential Gamma", "Logit Normal")
      )
    )
}))

# Calculate AUC values:
auc_comp_df <- df_paired |>
  group_by(Comparison, Distribution_clean, Method) |>
  filter(length(unique(Truth)) == 2) |>
  do(plotROC::calc_auc(
    ggplot(., aes(d = Truth, m = Score)) + geom_roc()
  )) |>
  ungroup() |>
  rename(AUROC = AUC)

# Fix formatting:
auc_comp_text <- auc_comp_df |>
  mutate(
    Method_short = ifelse(Method == "DiD Baseline", "DiD", Method),
    auc_label = paste0(Method_short, ": ", sprintf("%.3f", AUROC))
  ) |>
  group_by(Comparison) |>
  arrange(Distribution_clean, desc(Method), .by_group = TRUE) |>
  mutate(
    x = 0.60,
    y = seq(0.50, 0.12, length.out = n()) # Fits the 4 simplified lines
  ) |>
  ungroup()

# Plot:
prior_ROC_comp_auc_plot = ggplot(
  data = df_paired,
  mapping = aes(
    d = Truth, 
    m = Score, 
    color = Distribution_clean, 
    linetype = Method
  )
) +
  geom_roc(n.cuts = 0, size = 1.2, linealpha = 0.7) +
  geom_abline(
    slope = 1, 
    intercept = 0, 
    linetype = "dashed", 
    color = "grey10", 
    linewidth = 0.5
  ) +
  geom_text(
    data = auc_comp_text,
    aes(x = x, y = y, label = auc_label, color = Distribution_clean),
    inherit.aes = FALSE,
    size = 3.2,
    hjust = 0,
    vjust = 1
  ) +
  facet_wrap(
    ~ Comparison, 
    ncol = 2,
    labeller = labeller(Comparison = label_value)
  ) +
  scale_color_manual(values = c(
    "Beta"              = "grey10", 
    "Simplex"           = "deeppink", 
    "Bimodal Beta"      = "steelblue3", 
    "Exponential Gamma" = "orange", 
    "Logit Normal"      = "purple"
  )) +
  scale_linetype_manual(values = c(
    "MIMOSA2"      = "solid", 
    "DiD Baseline" = "dotted"
  )) +
  guides(
    color = guide_legend(
      order = 1, 
      title = "Distribution",
      override.aes = list(linewidth = 1.2, alpha = 1)
    ),
    linetype = guide_legend(
      order = 2, 
      title = "Model",
      override.aes = list(linewidth = 1.2, alpha = 1)
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Simulation performance of MIMOSA2.",
    subtitle = "ROC analysis across prior distributions.",
    x = "FPR",
    y = "TPR"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical"
  )

print(prior_ROC_comp_auc_plot)
# ggsave('_fig/prior_ROC_comp_auc_plot.pdf', plot = prior_ROC_comp_auc_plot, width = 8, height = 6)

# -----------------------------------------------------------------------------
# False discovery rate
# Points on or below the dashed line indicate properly controlled FDR
# 'fdr_plot_clean'
# Effect: small,med,large
# Res_prop: 10%,50%,90%
# -----------------------------------------------------------------------------
alpha_grid = seq(0.01, 0.20, by = 0.01)

# Compute FDR evaluation across prior distributions:
fdr_eval = map_dfr(alpha_grid, function(a) {
  results_continuous |>
    filter(!is.na(MIMOSA2_prob), !is.na(Truth)) |>
    group_by(Distribution, Effect, Cell_range) |>
    summarise(
      Nominal_FDR = a,
      Observed_FDR = {
        selected <- (1 - MIMOSA2_prob) <= a
        n_selected <- sum(selected, na.rm = TRUE)
        if (n_selected == 0) 0 else sum(selected & (Truth == 0), na.rm = TRUE) / n_selected
      },
      .groups = "drop"
    )
}) |>
  mutate(
    Distribution_clean = factor(
      Distribution,
      levels = c("SX", "BB", "Beta", "EG", "LN"),
      labels = c("Simplex", "Bimodal Beta", "Beta", "Exponential Gamma", "Logit Normal")
    ),
    Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low", "Very Low")),
    Effect_clean = paste0("Effect: ", formatC(Effect, format = "e", digits = 1))
  )

# Plot:
fdr_prior_plot = ggplot(
  data = fdr_eval, 
  mapping = aes(
    x = Nominal_FDR, 
    y = Observed_FDR, 
    color = Cell_range
  )
) +
  geom_abline(
    slope = 1, 
    intercept = 0, 
    linetype = "dashed", 
    color = "grey10", 
    linewidth = 0.5
  ) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  facet_grid(
    Distribution_clean ~ Effect_clean
  ) +
  scale_color_manual(values = c(
    "High"     = "deeppink", 
    "Medium"   = "steelblue3", 
    "Low"      = "orange",
    "Very Low" = "purple"
  )) +
  scale_x_continuous(
    breaks = seq(0, 0.20, 0.05), 
    labels = sprintf("%.2f", seq(0, 0.20, 0.05))
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.25, 0.05), 
    labels = sprintf("%.2f", seq(0, 0.25, 0.05))
  ) +
  # Use coord_cartesian to zoom without discarding summary data points:
  coord_cartesian(ylim = c(0, 0.25)) +
  guides(
    color = guide_legend(
      title = "Cell Count",
      override.aes = list(linewidth = 1.2, size = 2)
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "False discovery rate.",
    subtitle = "Points on or below the dashed line indicate properly controlled FDR.",
    x = "Nominal FDR",
    y = "Observed FDR"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8.5),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(fdr_prior_plot)
# ggsave('_fig/fdr_prior_plot.pdf',plot=fdr_prior_plot,width=10,height=10)
