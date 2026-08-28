# MIMOSA2 Simulation plots
# Isabella Lethbridge and Tayyeb Abrahams 
# 16 August 2026

#-------------------------------------------------------------------------------
# Load required libraries: 
library(tidyverse)
library(ggplot2)
library(cowplot)
library(plotROC)
library(Hmisc) 
library(dplyr)
library(tidyr)
library(knitr)

# Load data:
setwd("C:/Github/comp-2026-mimosa2")
load('_simulations/Simulation_2.0.Rdata')

#-------------------------------------------------------------------------------
# Sensitivity analysis of MIMOSA2
# TPR at 1% FDR 
# 'base_plot' 
# All simulation conditions 
#-------------------------------------------------------------------------------
# Filter data:
plot_data = subset(results_summary, Status == 'Success')

# Relabel 'Sparse' to 'Low': 
plot_data$Cell_range = factor(
  plot_data$Cell_range,
  levels = c("High", "Medium", "Sparse"),
  labels = c("High", "Medium", "Low")
)

# Format effect sizes:
unique_effects = sort(unique(plot_data$Effect), decreasing = TRUE)
plot_data$Effect_fact = factor(
  plot_data$Effect,
  levels = unique_effects,
  labels = formatC(unique_effects, format = "e", digits = 1)
)

# Relabel 'Res_prop': 
plot_data$Res_prop_clean = factor(
  plot_data$Res_prop,
  levels = c("Prop_0.10", "Prop_0.25", "Prop_0.50", "Prop_0.75", "Prop_0.90"),
  labels = c("10% responders", "25% responders", "50% responders", "75% responders", "90% responders")
)

# Plot: 
base_plot = ggplot(data = plot_data,
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
  facet_grid(
    Res_prop_clean ~ P, 
    labeller = labeller(Res_prop_clean = label_value, P = label_both)
  ) +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = c(0, 0.5, 1.0), 
                     expand = c(0.02, 0.02)) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "steelblue3", 
                                "Low" = "orange"),
                     breaks = c("Low", "Medium", "High")) +
  scale_fill_manual(values = c("High" = "deeppink", 
                               "Medium" = "steelblue3", 
                               "Low" = "orange"),
                    breaks = c("Low", "Medium", "High")) +
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

print(base_plot)
# ggsave('_fig/base_plot.pdf',plot=base_plot,width=10,height=10)

#-------------------------------------------------------------------------------
# Sensitivity analysis of MIMOSA2
# TPR at 1% FDR 
# 'base_plot' 
# P: 10,50,100
# Res_prop: 10%,50%,90%
#-------------------------------------------------------------------------------
# Filter data:
plot_data_clean = results_summary |>
  filter(Status == 'Success',
         P %in% c(10, 50, 100),       # Change P 
         Res_prop %in% c("Prop_0.10","Prop_0.50","Prop_0.90")  # Change Res_prop 
  )

# Relabel 'Sparse' to 'Low': 
plot_data_clean$Cell_range = factor(
  plot_data_clean$Cell_range,
  levels = c("High", "Medium", "Sparse"),
  labels = c("High", "Medium", "Low")
)

# Relabel 'Res_prop':
plot_data_clean$Res_prop_clean = factor(
  plot_data_clean$Res_prop,
  levels = c("Prop_0.10", "Prop_0.50", "Prop_0.90"),
  labels = c("10% Responders", "50% Responders", "90% Responders")
)

# Format effect size: 
unique_effects = sort(unique(plot_data_clean$Effect), decreasing = TRUE)
plot_data_clean$Effect_clean = factor(
  plot_data_clean$Effect,
  levels = unique_effects,
  labels = formatC(unique_effects, format = "e", digits = 1)
)

# Format sample size: 
plot_data_clean$P_label = factor(
  plot_data_clean$P,
  levels = c(10, 50, 100),
  labels = c("N = 10", "N = 50", "N = 100")
)

# Plot: 
base_plot_clean = ggplot(data = plot_data_clean,
                         mapping = aes(x = Effect_clean, 
                                       y = TPR_001, 
                                       group = Cell_range,
                                       color = Cell_range,
                                       fill = Cell_range)
                          ) +
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
  facet_grid(
    rows = vars(Res_prop_clean),
    cols = vars(P_label), 
    labeller = labeller(Res_prop_clean = label_value, P = label_value)
  ) +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = c(0, 0.5, 1.0), 
                     expand = c(0.02, 0.02)) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "steelblue3", 
                                "Low" = "orange"),
                     breaks = c("Low", "Medium", "High")) +
  scale_fill_manual(values = c("High" = "deeppink", 
                               "Medium" = "steelblue3", 
                               "Low" = "orange"),
                    breaks = c("Low", "Medium", "High")) +
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
 
print(base_plot_clean)
# ggsave('_fig/base_plot_clean.pdf',plot=base_plot_clean,width=8,height=6)

#-------------------------------------------------------------------------------
# Simulation performance of MIMOSA2
# ROC analysis 
# 'ROC_plot' 
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
    # Relabel 'Sparse' to 'Low':
    Cell_range = factor(
      Cell_range,
      levels = c("High", "Medium", "Sparse"),
      labels = c("High", "Medium", "Low")
    ),
    # Relabel 'Res_prop':
    Res_prop_clean = factor(
      Res_prop,
      levels = c("Prop_0.10", "Prop_0.25", "Prop_0.50", "Prop_0.75", "Prop_0.90"),
      labels = c("10% responders", "25% responders", "50% responders", "75% responders", "90% responders")
    )
  )

# Plot: 
ROC_plot = ggplot(
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
    Res_prop_clean ~ Effect, 
    labeller = labeller(
      Res_prop_clean = label_value, 
      Effect = label_both
    )
  ) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "steelblue3", 
                                "Low" = "orange")) +
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

print(ROC_plot)
# ggsave('_fig/ROC_plot.pdf',plot=ROC_plot,width=10,height=10)

#-------------------------------------------------------------------------------
# Simulation performance of MIMOSA2
# ROC analysis
# 'ROC_plot_clean' 
# Effect: small,med,large
# Res_prop: 10%,50%,90%
#-------------------------------------------------------------------------------
# Filter by effect size: 
all_effects = sort(unique(results_continuous$Effect))
selected_effects = all_effects[c(1, ceiling(length(all_effects)/2), length(all_effects))]

# Filter data: 
ROC_data_clean = results_continuous |>
  filter(
    !is.na(MIMOSA2_prob), 
    !is.na(DiD_GLM_prob),
    Res_prop %in% c("Prop_0.10","Prop_0.50","Prop_0.90"), 
    Effect %in% selected_effects        # Keep 3 representative effect sizes
  ) |>
  pivot_longer(
    cols = c(MIMOSA2_prob, DiD_GLM_prob),
    names_to = "Method",
    values_to = "Score"
  ) |>
  mutate(
    Method = ifelse(Method == "MIMOSA2_prob", "MIMOSA2", "DiD Baseline"),
    # Relabel 'Res_prop': 
    Res_prop_clean = factor(
      Res_prop,
      levels = c("Prop_0.10", "Prop_0.50", "Prop_0.90"),
      labels = c("10% Responders", "50% Responders", "90% Responders")
    ),
    # Relabel 'Sparse' to 'Low':
    Cell_range = factor(
      Cell_range,
      levels = c("High", "Medium", "Sparse"),
      labels = c("High", "Medium", "Low")
    ),
    # Format Effect labels cleanly:
    Effect_clean = factor(
      Effect,
      levels = selected_effects,
      labels = paste0("Effect: ", formatC(selected_effects, format = "e", digits = 1))
    )
  )

# Plot: 
ROC_plot_clean = ggplot(
  data = ROC_data_clean,
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
  # Facet across selected effect sizes in a single horizontal row:
  facet_grid(
    rows = vars(Res_prop_clean),
    cols = vars(Effect_clean)
  ) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "steelblue3", 
                                "Low" = "orange"),
                     breaks = c("Low", "Medium", "High")) +
  scale_linetype_manual(values = c("MIMOSA2" = "solid", 
                                   "DiD Baseline" = "dotted")) +
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
  theme_bw(base_size = 11) +
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

print(ROC_plot_clean)
# ggsave('_fig/ROC_plot_clean.pdf',plot=ROC_plot_clean,width=8,height=6)

#-------------------------------------------------------------------------------
# AUROC Calculation
#-------------------------------------------------------------------------------
# Calculate AUROC:
# AUROC = ROC_data_prepared |>
#   group_by(Res_prop, P, Cell_range, Effect, Method) |>
#   filter(length(unique(Truth)) == 2) |>
#   do(plotROC::calc_auc(
#     ggplot(., aes(d = Truth, m = Score)) + geom_roc()
#   )) |>
#   ungroup() |>
#   rename(AUROC = AUC)

# write.table(AUROC, "auc.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# View(AUROC)

# Summarise: 
# AUROC matrix across Effect and Cell Range: 
AUROC_matrix = AUROC |>
  filter(Res_prop %in% c("Prop_0.10", "Prop_0.50", "Prop_0.90") | 
           Res_prop %in% c("0.10", "0.50", "0.90") | 
           Res_prop %in% c(0.10, 0.50, 0.90)) |>
  mutate(
    Res_prop_clean = factor(
      Res_prop,
      levels = c("Prop_0.10","Prop_0.50","Prop_0.90"),
      labels = c("10% responders", "50% responders", "90% responders")
    ), 
    Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")),
    Effect_clean = paste0(formatC(Effect, format = "e", digits = 1))
  ) |>
  group_by(Res_prop_clean, Effect_clean, Cell_range, Method) |>
  summarise(Mean_AUROC = sprintf("%.3f", 
                                 mean(AUROC, na.rm = TRUE)), 
            .groups = "drop") |>
  pivot_wider(
    names_from = Method,
    values_from = Mean_AUROC
  ) |>
  arrange(Effect_clean, Cell_range)

# Save as file: 
# write.table(AUROC_matrix, 
#             "auc_summary_matrix.txt", 
#             sep = "\t", 
#             row.names = FALSE, 
#             quote = FALSE)

# Table: AUROC values
kable(AUROC_matrix, 
      caption = 'AUROC values',
      col.names = c("Responders","Effect Size", "Cell Count", "MIMOSA2", "DiD Baseline"))

# Table: AUROC values
# 
# |Effect Size |Cell Count |MIMOSA2 |DiD Baseline |
# |:-----------|:----------|:-------|:------------|
# |1.0e-03     |High       |0.916   |0.954        |
# |1.0e-03     |Medium     |0.903   |0.943        |
# |1.0e-03     |Low        |0.854   |0.903        |
# |1.2e-04     |High       |0.876   |0.909        |
# |1.2e-04     |Medium     |0.849   |0.881        |
# |1.2e-04     |Low        |0.759   |0.783        |
# |2.5e-04     |High       |0.892   |0.930        |
# |2.5e-04     |Medium     |0.869   |0.903        |
# |2.5e-04     |Low        |0.783   |0.820        |
# |5.0e-04     |High       |0.903   |0.944        |
# |5.0e-04     |Medium     |0.887   |0.929        |
# |5.0e-04     |Low        |0.811   |0.854        |
# |6.2e-05     |High       |0.872   |0.900        |
# |6.2e-05     |Medium     |0.845   |0.869        |
# |6.2e-05     |Low        |0.749   |0.775        |

# Table: AUROC values
# 
# |Responders     |Effect Size |Cell Count |MIMOSA2 |DiD Baseline |
# |:--------------|:-----------|:----------|:-------|:------------|
# |10% responders |1.0e-03     |High       |0.920   |0.956        |
# |50% responders |1.0e-03     |High       |0.910   |0.950        |
# |90% responders |1.0e-03     |High       |0.916   |0.952        |
# |10% responders |1.0e-03     |Medium     |0.918   |0.947        |
# |50% responders |1.0e-03     |Medium     |0.899   |0.942        |
# |90% responders |1.0e-03     |Medium     |0.898   |0.940        |
# |10% responders |1.0e-03     |Low        |0.871   |0.895        |
# |50% responders |1.0e-03     |Low        |0.842   |0.899        |
# |90% responders |1.0e-03     |Low        |0.853   |0.906        |
# |10% responders |1.2e-04     |High       |0.876   |0.916        |
# |50% responders |1.2e-04     |High       |0.875   |0.906        |
# |90% responders |1.2e-04     |High       |0.872   |0.900        |
# |10% responders |1.2e-04     |Medium     |0.848   |0.880        |
# |50% responders |1.2e-04     |Medium     |0.837   |0.877        |
# |90% responders |1.2e-04     |Medium     |0.861   |0.889        |
# |10% responders |1.2e-04     |Low        |0.765   |0.781        |
# |50% responders |1.2e-04     |Low        |0.772   |0.794        |
# |90% responders |1.2e-04     |Low        |0.759   |0.781        |
# |10% responders |2.5e-04     |High       |0.895   |0.930        |
# |50% responders |2.5e-04     |High       |0.885   |0.928        |
# |90% responders |2.5e-04     |High       |0.893   |0.932        |
# |10% responders |2.5e-04     |Medium     |0.868   |0.908        |
# |50% responders |2.5e-04     |Medium     |0.871   |0.903        |
# |90% responders |2.5e-04     |Medium     |0.874   |0.901        |
# |10% responders |2.5e-04     |Low        |0.782   |0.803        |
# |50% responders |2.5e-04     |Low        |0.778   |0.821        |
# |90% responders |2.5e-04     |Low        |0.778   |0.824        |
# |10% responders |5.0e-04     |High       |0.914   |0.952        |
# |50% responders |5.0e-04     |High       |0.899   |0.945        |
# |90% responders |5.0e-04     |High       |0.900   |0.941        |
# |10% responders |5.0e-04     |Medium     |0.883   |0.930        |
# |50% responders |5.0e-04     |Medium     |0.884   |0.927        |
# |90% responders |5.0e-04     |Medium     |0.890   |0.928        |
# |10% responders |5.0e-04     |Low        |0.815   |0.841        |
# |50% responders |5.0e-04     |Low        |0.805   |0.860        |
# |90% responders |5.0e-04     |Low        |0.823   |0.867        |
# |10% responders |6.2e-05     |High       |0.874   |0.901        |
# |50% responders |6.2e-05     |High       |0.870   |0.898        |
# |90% responders |6.2e-05     |High       |0.876   |0.902        |
# |10% responders |6.2e-05     |Medium     |0.840   |0.865        |
# |50% responders |6.2e-05     |Medium     |0.847   |0.874        |
# |90% responders |6.2e-05     |Medium     |0.847   |0.870        |
# |10% responders |6.2e-05     |Low        |0.745   |0.761        |
# |50% responders |6.2e-05     |Low        |0.749   |0.772        |
# |90% responders |6.2e-05     |Low        |0.759   |0.784        |

# Table: AUROC values (MIMOSA2, DiD Baseline)
AUROC_matrix_2 <- AUROC |>
  filter(Res_prop %in% c("Prop_0.10", "Prop_0.50", "Prop_0.90") | 
           Res_prop %in% c("0.10", "0.50", "0.90") | 
           Res_prop %in% c(0.10, 0.50, 0.90)) |>
  mutate(
    Res_prop_clean = factor(
      Res_prop,
      levels = c("Prop_0.10", "Prop_0.50", "Prop_0.90", "0.10", "0.50", "0.90", 0.10, 0.50, 0.90),
      labels = c("10% responders", "50% responders", "90% responders",
                 "10% responders", "50% responders", "90% responders",
                 "10% responders", "50% responders", "90% responders")
    ),
    Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")),
    Effect_clean = formatC(Effect, format = "e", digits = 1),
    Method = factor(Method, levels = c("MIMOSA2", "DiD Baseline")) # Ensures MIMOSA2 is first
  ) |>
  group_by(Res_prop_clean, Effect_clean, Cell_range, Method) |>
  summarise(
    Mean_AUROC = sprintf("%.3f", mean(AUROC, na.rm = TRUE)), 
    .groups = "drop"
  ) |>
  arrange(Res_prop_clean, Effect_clean, Cell_range, Method) |>
  # Combine MIMOSA2 and DiD Baseline into "val1, val2"
  group_by(Res_prop_clean, Effect_clean, Cell_range) |>
  summarise(
    Combined_AUROC = paste(Mean_AUROC, collapse = ", "),
    .groups = "drop"
  ) |>
  # Pivot Cell_range to columns
  pivot_wider(
    names_from = Cell_range,
    values_from = Combined_AUROC
  ) |>
  arrange(Effect_clean,Res_prop_clean)

# Save as file: 
# write.table(
#   AUROC_matrix_2, 
#   "auc_summary_matrix_2.txt", 
#   sep = "\t", 
#   row.names = FALSE, 
#   quote = FALSE
# )

# Render Table: 
kable(
  AUROC_matrix_2, 
  caption = 'AUROC values (MIMOSA2, DiD Baseline)',
  col.names = c("Responders", "Effect Size", "High", "Medium", "Low")
)

# Table: AUROC values (MIMOSA2, DiD Baseline)
# 
# |Effect Size |High         |Medium       |Low          |
# |:-----------|:------------|:------------|:------------|
# |1.0e-03     |0.954, 0.916 |0.943, 0.903 |0.903, 0.854 |
# |1.2e-04     |0.909, 0.876 |0.881, 0.849 |0.783, 0.759 |
# |2.5e-04     |0.930, 0.892 |0.903, 0.869 |0.820, 0.783 |
# |5.0e-04     |0.944, 0.903 |0.929, 0.887 |0.854, 0.811 |
# |6.2e-05     |0.900, 0.872 |0.869, 0.845 |0.775, 0.749 |

# Table: AUROC values (MIMOSA2, DiD Baseline)
# 
# |Responders     |Effect Size |High         |Medium       |Low          |
# |:--------------|:-----------|:------------|:------------|:------------|
# |10% responders |1.0e-03     |0.956, 0.920 |0.947, 0.918 |0.895, 0.871 |
# |50% responders |1.0e-03     |0.950, 0.910 |0.942, 0.899 |0.899, 0.842 |
# |90% responders |1.0e-03     |0.952, 0.916 |0.940, 0.898 |0.906, 0.853 |
# |10% responders |1.2e-04     |0.916, 0.876 |0.880, 0.848 |0.781, 0.765 |
# |50% responders |1.2e-04     |0.906, 0.875 |0.877, 0.837 |0.794, 0.772 |
# |90% responders |1.2e-04     |0.900, 0.872 |0.889, 0.861 |0.781, 0.759 |
# |10% responders |2.5e-04     |0.930, 0.895 |0.908, 0.868 |0.803, 0.782 |
# |50% responders |2.5e-04     |0.928, 0.885 |0.903, 0.871 |0.821, 0.778 |
# |90% responders |2.5e-04     |0.932, 0.893 |0.901, 0.874 |0.824, 0.778 |
# |10% responders |5.0e-04     |0.952, 0.914 |0.930, 0.883 |0.841, 0.815 |
# |50% responders |5.0e-04     |0.945, 0.899 |0.927, 0.884 |0.860, 0.805 |
# |90% responders |5.0e-04     |0.941, 0.900 |0.928, 0.890 |0.867, 0.823 |
# |10% responders |6.2e-05     |0.901, 0.874 |0.865, 0.840 |0.761, 0.745 |
# |50% responders |6.2e-05     |0.898, 0.870 |0.874, 0.847 |0.772, 0.749 |
# |90% responders |6.2e-05     |0.902, 0.876 |0.870, 0.847 |0.784, 0.759 |

# Simulation count table:
Sim_count_matrix <- ROC_data_prepared |>
  filter(Res_prop %in% c("Prop_0.10", "Prop_0.50", "Prop_0.90", "0.10", "0.50", "0.90", 0.10, 0.50, 0.90)) |>
  mutate(
    Res_prop_clean = factor(
      Res_prop,
      levels = c("Prop_0.10", "Prop_0.50", "Prop_0.90", "0.10", "0.50", "0.90", 0.10, 0.50, 0.90),
      labels = c("10% responders", "50% responders", "90% responders",
                 "10% responders", "50% responders", "90% responders",
                 "10% responders", "50% responders", "90% responders")
    ),
    Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")),
    Effect_clean = formatC(Effect, format = "e", digits = 1),
    Method = factor(Method, levels = c("MIMOSA2", "DiD Baseline")) # Ensures MIMOSA2 is first
  ) |>
  group_by(Res_prop_clean, Effect_clean, Cell_range, Method) |>
  summarise(
    N_sims = n(), 
    .groups = "drop"
  ) |>
  arrange(Res_prop_clean, Effect_clean, Cell_range, Method) |>
  # Combine MIMOSA2 and DiD Baseline simulation counts
  group_by(Res_prop_clean, Effect_clean, Cell_range) |>
  summarise(
    Combined_Count = paste(N_sims, collapse = ", "),
    .groups = "drop"
  ) |>
  # Pivot Cell_range to columns
  pivot_wider(
    names_from = Cell_range,
    values_from = Combined_Count
  ) |>
  arrange(Res_prop_clean, Effect_clean)

# Save as file:
# write.table(
#   Sim_count_matrix, 
#   "simulation_counts_matrix.txt", 
#   sep = "\t", 
#   row.names = FALSE, 
#   quote = FALSE
# )

# Render Table:
kable(
  Sim_count_matrix, 
  caption = 'Simulation Counts (MIMOSA2, DiD Baseline)',
  col.names = c("Responders", "Effect Size", "High", "Medium", "Low")
)

# Table: Simulation Counts (MIMOSA2, DiD Baseline)
# 
# |Effect Size |High         |Medium       |Low          |
# |:-----------|:------------|:------------|:------------|
# |1.0e-03     |42750, 42750 |42750, 42750 |42750, 42750 |
# |1.2e-04     |42750, 42750 |42750, 42750 |42750, 42750 |
# |2.5e-04     |42750, 42750 |42750, 42750 |42740, 42740 |
# |5.0e-04     |42750, 42750 |42750, 42750 |42750, 42750 |
# |6.2e-05     |42750, 42750 |42750, 42750 |42750, 42750 |


# Table: Simulation Counts (MIMOSA2, DiD Baseline)
# 
# |Responders     |Effect Size |High       |Medium     |Low        |
# |:--------------|:-----------|:----------|:----------|:----------|
# |10% responders |1.0e-03     |8550, 8550 |8550, 8550 |8550, 8550 |
# |10% responders |1.2e-04     |8550, 8550 |8550, 8550 |8550, 8550 |
# |10% responders |2.5e-04     |8550, 8550 |8550, 8550 |8540, 8540 |
# |10% responders |5.0e-04     |8550, 8550 |8550, 8550 |8550, 8550 |
# |10% responders |6.2e-05     |8550, 8550 |8550, 8550 |8550, 8550 |
# |50% responders |1.0e-03     |8550, 8550 |8550, 8550 |8550, 8550 |
# |50% responders |1.2e-04     |8550, 8550 |8550, 8550 |8550, 8550 |
# |50% responders |2.5e-04     |8550, 8550 |8550, 8550 |8550, 8550 |
# |50% responders |5.0e-04     |8550, 8550 |8550, 8550 |8550, 8550 |
# |50% responders |6.2e-05     |8550, 8550 |8550, 8550 |8550, 8550 |
# |90% responders |1.0e-03     |8550, 8550 |8550, 8550 |8550, 8550 |
# |90% responders |1.2e-04     |8550, 8550 |8550, 8550 |8550, 8550 |
# |90% responders |2.5e-04     |8550, 8550 |8550, 8550 |8550, 8550 |
# |90% responders |5.0e-04     |8550, 8550 |8550, 8550 |8550, 8550 |
# |90% responders |6.2e-05     |8550, 8550 |8550, 8550 |8550, 8550 |

# Individual confidence intervals: 
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
  # value_cols <- setdiff(colnames(auc_matrix), colnames(auc_matrix)[1])
  value_cols <- c("High", "Medium", "Low")
  
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
  caption = 'Δ AUROC (MIMOSA2 - DiD Baseline) with 95% Confidence Intervals',
  col.names = c("Responders", "Effect Size", "High", "Medium", "Low")
)

# Table: Δ AUROC (MIMOSA2 - DiD Baseline) with 95% Confidence Intervals
# 
# |Effect Size |High                 |Medium               |Low                  |
# |:-----------|:--------------------|:--------------------|:--------------------|
# |1.0e-03     |0.038 [0.036, 0.040] |0.040 [0.037, 0.043] |0.049 [0.046, 0.052] |
# |1.2e-04     |0.033 [0.030, 0.036] |0.032 [0.028, 0.036] |0.024 [0.020, 0.028] |
# |2.5e-04     |0.038 [0.035, 0.041] |0.034 [0.031, 0.037] |0.037 [0.033, 0.041] |
# |5.0e-04     |0.041 [0.038, 0.044] |0.042 [0.039, 0.045] |0.043 [0.039, 0.047] |
# |6.2e-05     |0.028 [0.025, 0.031] |0.024 [0.020, 0.028] |0.026 [0.021, 0.031] |


# Table: Δ AUROC (MIMOSA2 - DiD Baseline) with 95% Confidence Intervals
# 
# |Responders     |Effect Size |High                 |Medium               |Low                  |
# |:--------------|:-----------|:--------------------|:--------------------|:--------------------|
# |10% responders |1.0e-03     |0.036 [0.031, 0.041] |0.029 [0.023, 0.035] |0.024 [0.017, 0.031] |
# |50% responders |1.0e-03     |0.040 [0.034, 0.046] |0.043 [0.037, 0.049] |0.057 [0.049, 0.065] |
# |90% responders |1.0e-03     |0.036 [0.030, 0.042] |0.042 [0.036, 0.048] |0.053 [0.046, 0.060] |
# |10% responders |1.2e-04     |0.040 [0.033, 0.047] |0.032 [0.024, 0.040] |0.016 [0.006, 0.026] |
# |50% responders |1.2e-04     |0.031 [0.024, 0.038] |0.040 [0.032, 0.048] |0.022 [0.012, 0.032] |
# |90% responders |1.2e-04     |0.028 [0.021, 0.035] |0.028 [0.020, 0.036] |0.022 [0.012, 0.032] |
# |10% responders |2.5e-04     |0.035 [0.029, 0.041] |0.040 [0.033, 0.047] |0.021 [0.011, 0.031] |
# |50% responders |2.5e-04     |0.043 [0.036, 0.050] |0.032 [0.025, 0.039] |0.043 [0.034, 0.052] |
# |90% responders |2.5e-04     |0.039 [0.033, 0.045] |0.027 [0.020, 0.034] |0.046 [0.037, 0.055] |
# |10% responders |5.0e-04     |0.038 [0.032, 0.044] |0.047 [0.040, 0.054] |0.026 [0.017, 0.035] |
# |50% responders |5.0e-04     |0.046 [0.040, 0.052] |0.043 [0.036, 0.050] |0.055 [0.046, 0.064] |
# |90% responders |5.0e-04     |0.041 [0.035, 0.047] |0.038 [0.031, 0.045] |0.044 [0.036, 0.052] |
# |10% responders |6.2e-05     |0.027 [0.020, 0.034] |0.025 [0.017, 0.033] |0.016 [0.006, 0.026] |
# |50% responders |6.2e-05     |0.028 [0.021, 0.035] |0.027 [0.019, 0.035] |0.023 [0.013, 0.033] |
# |90% responders |6.2e-05     |0.026 [0.019, 0.033] |0.023 [0.015, 0.031] |0.025 [0.015, 0.035] |
  
# -----------------------------------------------------------------------------
# Simulation performance of MIMOSA2
# ROC analysis
# 'ROC_plot_overlay' 
# Effect: small,med,large
# Res_prop: 10%,50%,90%
# -----------------------------------------------------------------------------
# Filter data: 
ROC_data_filtered <- ROC_data_clean |>
  filter(Res_prop %in% c("Prop_0.10", "Prop_0.50", "Prop_0.90", "0.10", "0.50", "0.90", 0.10, 0.50, 0.90)) |>
  mutate(
    Res_prop_clean = factor(
      Res_prop,
      levels = c("Prop_0.10", "Prop_0.50", "Prop_0.90", "0.10", "0.50", "0.90", 0.10, 0.50, 0.90),
      labels = c("10% responders", "50% responders", "90% responders",
                 "10% responders", "50% responders", "90% responders",
                 "10% responders", "50% responders", "90% responders")
    )
  )

# Calculate AUC: 
auroc_df = ROC_data_filtered |>
  group_by(Res_prop_clean, Effect_clean, Cell_range, Method) |>
  filter(length(unique(Truth)) == 2) |>
  do(plotROC::calc_auc(
    ggplot(., aes(d = Truth, m = Score)) + geom_roc()
  )) |>
  ungroup() |>
  rename(AUROC = AUC)

# Format AUC text overlay (Template: MIMOSA2: 0.XXX / DiD: 0.XXX)
auc_text_df = auroc_df |>
  mutate(
    # Clean model label (DiD Baseline -> DiD)
    Method_short = ifelse(Method == "DiD Baseline", "DiD", Method),
    # Format string to 3 decimal places
    auc_label = paste0(Method_short, ": ", sprintf("%.3f", AUROC))
  ) |>
  group_by(Res_prop_clean, Effect_clean) |>
  arrange(desc(Method), desc(Cell_range), .by_group = TRUE) |>
  mutate(
    x = 0.60,                             # This shift text box left/right                           
    y = seq(0.52, 0.08, length.out = n()) # This shifts text box up/down   
  ) |>
  ungroup()

# Plot: 
ROC_plot_overlay = ggplot(
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
  facet_grid(Res_prop_clean ~ Effect_clean) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "steelblue3", 
                                "Low" = "orange")) +
  scale_linetype_manual(values = c("MIMOSA2" = "solid", 
                                   "DiD Baseline" = "dotted")) +
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

print(ROC_plot_overlay)
# ggsave('_fig/ROC_plot_overlay.pdf',plot=ROC_plot_overlay,width=8,height=6)

# -----------------------------------------------------------------------------
# False discovery rate
# Points on or below the dashed line indicate properly controlled FDR
# 'fdr_plot' 
# -----------------------------------------------------------------------------
# Define grid: 
alpha_grid = seq(0.01, 0.20, by = 0.01)

fdr_eval = map_dfr(alpha_grid, function(a) {
  results_continuous |>
    filter(!is.na(MIMOSA2_prob), !is.na(Truth)) |>
    group_by(Res_prop, Effect, Cell_range) |>
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
    # Relabel 'Sparse' -> 'Low':
    Cell_range = factor(
      Cell_range,
      levels = c("High", "Medium", "Sparse"),
      labels = c("High", "Medium", "Low")
    ),
    # Clean responder proportion labels:
    Res_prop_clean = factor(
      Res_prop,
      levels = c("Prop_0.10", "Prop_0.25", "Prop_0.50", "Prop_0.75", "Prop_0.90"),
      labels = c("10% responders", "25% responders", "50% responders", "75% responders", "90% responders")
    )
  )

# Plot: 
fdr_plot = ggplot(fdr_eval, 
                  aes(x = Nominal_FDR, 
                      y = Observed_FDR, 
                      color = Cell_range)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey10", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  facet_grid(
    Res_prop_clean ~ Effect,
    labeller = labeller(Res_prop_clean = label_value, Effect = label_both)
  ) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "steelblue3", 
                                "Low" = "orange")) +
  scale_x_continuous(breaks = seq(0, 0.20, 0.05), 
                     labels = sprintf("%.2f", seq(0, 0.20, 0.05))) +
  scale_y_continuous(limits = c(0, 0.25), 
                     breaks = seq(0, 0.25, 0.05), 
                     labels = sprintf("%.2f", seq(0, 0.25, 0.05))) +
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
    x = "Nominal FDR ",
    y = "Observed FDF"
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

print(fdr_plot)
# ggsave('_fig/fdr_plot.pdf',plot=fdr_plot,width=10,height=10)

#-------------------------------------------------------------------------------
# False discovery rate
# Points on or below the dashed line indicate properly controlled FDR
# 'fdr_plot_clean'
# Effect: small,med,large
# Res_prop: 10%,50%,90%
#-------------------------------------------------------------------------------
# Filter data: 
all_effects = sort(unique(fdr_eval$Effect))
selected_effects = all_effects[c(1, ceiling(length(all_effects)/2), 
                                 length(all_effects))]

fdr_eval_clean = fdr_eval |>
  filter(
    Res_prop %in% c("Prop_0.10", "Prop_0.50", "Prop_0.90"),
    Effect %in% selected_effects
  ) |>
  mutate(
    Res_prop_clean = factor(
      Res_prop,
      levels = c("Prop_0.10", "Prop_0.50", "Prop_0.90", "0.10", "0.50", "0.90", 0.10, 0.50, 0.90),
      labels = c("10% responders", "50% responders", "90% responders",
                 "10% responders", "50% responders", "90% responders",
                 "10% responders", "50% responders", "90% responders")
    ),
    Effect_clean = factor(
      Effect,
      levels = selected_effects,
      labels = paste0("Effect: ", formatC(selected_effects, format = "e", digits = 1))
    )
  )

# Plot:
fdr_plot_clean = ggplot(fdr_eval_clean, 
                        aes(x = Nominal_FDR, 
                            y = Observed_FDR, 
                            color = Cell_range)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey10", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  facet_grid(Res_prop_clean ~ Effect_clean) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "steelblue3", 
                                "Low" = "orange")) +
  scale_x_continuous(breaks = seq(0, 0.20, 0.05), 
                     labels = sprintf("%.2f", seq(0, 0.20, 0.05))) +
  scale_y_continuous(limits = c(0, 0.25), 
                     breaks = seq(0, 0.25, 0.05), 
                     labels = sprintf("%.2f", seq(0, 0.25, 0.05))) +
  guides(
    color = guide_legend(
      title = "Cell Count",
      override.aes = list(linewidth = 1.2, size = 2)
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "False discovery rate.",
    subtitle = "Points on or below the dashed line indicate properly controlled FDR",
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

print(fdr_plot_clean)
# ggsave('_fig/fdr_plot_clean.pdf',plot=fdr_plot_clean,width=8,height=8)
