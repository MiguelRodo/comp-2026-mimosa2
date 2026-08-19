# MIMOSA2 Simulation plots
# Isabella Lethbridge and Tayyeb Abrahams 
# 16 August 2026

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
load('_simulations/Simulation_2.0.Rdata')

#-------------------------------------------------------------------------------
# TPR (sensitivity) vs 1% Nominal FDR 
# This gives a 5x6 grid 
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

# Relabel Res_prop: 
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
                                "Medium" = "blue", 
                                "Low" = "purple")) +
  scale_fill_manual(values = c("High" = "deeppink", 
                               "Medium" = "blue", 
                               "Low" = "purple")) +
  labs(
    title = "Sensitivity Analysis of MIMOSA2.",
    subtitle = "True positive rate (TPR) at 1% nominal false discovery rate (FDR) threshold.",
    x = "Effect size",
    y = "Sensitivity (TPR)",
    color = "Cell Count",
    fill = "Cell Count"
  ) +
  theme_minimal(base_size = 11) +
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
# TPR (sensitivity) vs 1% Nominal FDR 
# This plots a 1x3 grid 
# Plots graph for only P=10,50,100
# Plots for only 50% responders 
#-------------------------------------------------------------------------------
# Filter data:
plot_data_clean = results_summary |>
  filter(Status == 'Success',
         P %in% c(10, 50, 100),       # Change P 
         Res_prop == "Prop_0.50"      # Change Res_prop 
  )

# Relabel 'Sparse' to 'Low': 
plot_data_clean$Cell_range = factor(
  plot_data_clean$Cell_range,
  levels = c("High", "Medium", "Sparse"),
  labels = c("High", "Medium", "Low")
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
  labels = c("P: 10", "P: 50", "P: 100")
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
    Res_prop_clean ~ P, 
    labeller = labeller(Res_prop_clean = label_value, P = label_both)
  ) +
  # Facet horizontally across sample sizes (1 x 3 grid): 
  facet_wrap(~ P_label, 
             nrow = 1) +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = c(0, 0.5, 1.0), 
                     expand = c(0.02, 0.02)) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "blue", 
                                "Low" = "purple")) +
  scale_fill_manual(values = c("High" = "deeppink", 
                               "Medium" = "blue", 
                               "Low" = "purple")) +
  labs(
    title = "Sensitivity Analysis of MIMOSA2.",
    subtitle = "True positive rate (TPR) at 1% nominal FDR for 50% repsonders.",
    x = "Effect size",
    y = "Sensitivity (TPR)",
    color = "Cell Count",
    fill = "Cell Count"
  ) +
  theme_minimal(base_size = 11) +
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
# These dimensions need to be changed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ggsave('_fig/base_plot_clean.pdf',plot=base_plot_clean,width=10,height=10)

#-------------------------------------------------------------------------------
# ROC plots
# This gives a 5x5 grid 
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
  geom_roc(n.cuts = 0, size = 0.8, alpha = 0.7) +
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
                                "Medium" = "blue", 
                                "Low" = "purple")) +
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
  theme_bw(base_size = 11) +
  labs(
    title = "MIMOSA2 simulation performance.",
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
# ROC plot
# This gives 1x3 grid
# Effect sizes: 6.2e-05,2.5e-04,1.0e-03 
# Res_prop = 0.5
#-------------------------------------------------------------------------------
# Filter by effect size: 
all_effects = sort(unique(results_continuous$Effect))
selected_effects = all_effects[c(1, ceiling(length(all_effects)/2), length(all_effects))]

# Filter data: 
ROC_data_clean = results_continuous |>
  filter(
    !is.na(MIMOSA2_prob), 
    !is.na(DiD_GLM_prob),
    Res_prop == "Prop_0.50",            # Fix to baseline 50% responders
    Effect %in% selected_effects        # Keep 3 representative effect sizes
  ) |>
  pivot_longer(
    cols = c(MIMOSA2_prob, DiD_GLM_prob),
    names_to = "Method",
    values_to = "Score"
  ) |>
  mutate(
    Method = ifelse(Method == "MIMOSA2_prob", "MIMOSA2", "DiD Baseline"),
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
  geom_roc(n.cuts = 0, size = 0.8, alpha = 0.7) +
  geom_abline(
    slope = 1, 
    intercept = 0, 
    linetype = "dashed", 
    color = "grey10", 
    linewidth = 0.5
  ) +
  # Facet across selected effect sizes in a single horizontal row:
  facet_wrap(~ Effect_clean, nrow = 1) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "blue", 
                                "Low" = "purple")) +
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
    title = "MIMOSA2 Simulation Performance",
    subtitle = "ROC Analysis for 50% responders.",
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
# These dimensions need to be changed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ggsave('_fig/ROC_plot_clean.pdf',plot=ROC_plot_clean,width=10,height=10)

#-------------------------------------------------------------------------------
# AUROC Calculation
#-------------------------------------------------------------------------------
# Calculate AUROC:
AUROC = ROC_data_prepared |>
  group_by(Res_prop, P, Cell_range, Effect, Method) |>
  filter(length(unique(Truth)) == 2) |>
  do(plotROC::calc_auc(
    ggplot(., aes(d = Truth, m = Score)) + geom_roc()
  )) |>
  ungroup() |>
  rename(AUROC = AUC)

write.table(AUROC, "auc.txt", sep = "\t", row.names = FALSE, quote = FALSE)
View(AUROC)

# Summarise: 
# AUROC matrix across Effect and Cell Range: 
AUROC_matrix = AUROC |>
  mutate(Cell_range = factor(Cell_range, 
                             levels = c("High", "Medium", "Low")),
    Effect_clean = paste0(formatC(Effect, format = "e", digits = 1))
  ) |>
  group_by(Effect_clean, Cell_range, Method) |>
  summarise(Mean_AUROC = sprintf("%.3f", 
                                 mean(AUROC, na.rm = TRUE)), 
            .groups = "drop") |>
  pivot_wider(
    names_from = Method,
    values_from = Mean_AUROC
  ) |>
  arrange(Effect_clean, Cell_range)

# Save as file: 
write.table(AUROC_matrix, 
            "auc_summary_matrix.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Table: 
kable(AUROC_matrix, 
      caption = 'AUROC values',
      col.names = c("Effect Size", "Cell Count", "MIMOSA2", "DiD Baseline"))

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
AUROC_matrix_2 <- AUROC |>
  mutate(
    Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")),
    Effect_clean = formatC(Effect, format = "e", digits = 1),
    Method = factor(Method, levels = c("MIMOSA2", "DiD Baseline")) # Ensures MIMOSA2 is first
  ) |>
  group_by(Effect_clean, Cell_range, Method) |>
  summarise(
    Mean_AUROC = sprintf("%.3f", mean(AUROC, na.rm = TRUE)), 
    .groups = "drop"
  ) |>
  arrange(Effect_clean, Cell_range, Method) |>
  # Combine MIMOSA2 and DiD Baseline into "val1, val2"
  group_by(Effect_clean, Cell_range) |>
  summarise(
    Combined_AUROC = paste(Mean_AUROC, collapse = ", "),
    .groups = "drop"
  ) |>
  # Pivot Cell_range to columns
  pivot_wider(
    names_from = Cell_range,
    values_from = Combined_AUROC
  ) |>
  arrange(Effect_clean)

# Save as file: 
write.table(
  AUROC_matrix_2, 
  "auc_summary_matrix_2.txt", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

# Render Table: 
kable(
  AUROC_matrix_2, 
  caption = 'AUROC values (MIMOSA2, DiD Baseline)',
  col.names = c("Effect Size", "High", "Medium", "Low")
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

# Simulation count table:
Sim_count_matrix <- ROC_data_prepared |>
  mutate(
    Cell_range = factor(Cell_range, levels = c("High", "Medium", "Low")),
    Effect_clean = formatC(Effect, format = "e", digits = 1),
    Method = factor(Method, levels = c("MIMOSA2", "DiD Baseline")) # Ensures MIMOSA2 is first
  ) |>
  group_by(Effect_clean, Cell_range, Method) |>
  summarise(
    N_sims = n(), 
    .groups = "drop"
  ) |>
  arrange(Effect_clean, Cell_range, Method) |>
  # Combine MIMOSA2 and DiD Baseline simulation counts
  group_by(Effect_clean, Cell_range) |>
  summarise(
    Combined_Count = paste(N_sims, collapse = ", "),
    .groups = "drop"
  ) |>
  # Pivot Cell_range to columns
  pivot_wider(
    names_from = Cell_range,
    values_from = Combined_Count
  ) |>
  arrange(Effect_clean)

# Save as file:
write.table(
  Sim_count_matrix, 
  "simulation_counts_matrix.txt", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

# Render Table:
kable(
  Sim_count_matrix, 
  caption = 'Simulation Counts (MIMOSA2, DiD Baseline)',
  col.names = c("Effect Size", "High", "Medium", "Low")
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
  value_cols <- setdiff(colnames(auc_matrix), colnames(auc_matrix)[1])
  
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
  col.names = c("Effect Size", "High", "Medium", "Low")
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
  
# -----------------------------------------------------------------------------
# Overlay AUC values to clean plot:
# -----------------------------------------------------------------------------
# Calculate AUC: 
auroc_df = ROC_data_clean |>
  group_by(Effect_clean, Cell_range, Method) |>
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
  group_by(Effect_clean) |>
  arrange(desc(Method), desc(Cell_range), .by_group = TRUE) |>
  mutate(
    x = 0.55,                                      
    y = seq(0.42, 0.07, length.out = n())          
  ) |>
  ungroup()

# Plot: 
ROC_plot_overlay = ggplot(
  data = ROC_data_clean,
  mapping = aes(
    d = Truth, 
    m = Score, 
    color = Cell_range, 
    linetype = Method, 
    group = interaction(Method, Cell_range)
  )
) +
  geom_roc(n.cuts = 0, size = 0.8, alpha = 0.7) +
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
    size = 3.2,
    hjust = 0,
    vjust = 1
  ) +
  facet_wrap(~ Effect_clean, nrow = 1) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "blue", 
                                "Low" = "purple")) +
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
  theme_bw(base_size = 11) +
  labs(
    title = "MIMOSA2 Simulation Performance",
    subtitle = "AUROC Analysis for 50% responders.",
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
# These dimensions need to be changed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Not happy with how the AUC values look but that is a quick fix 
ggsave('_fig/ROC_plot_overlay.pdf',plot=ROC_plot_overlay,width=10,height=10)

# -----------------------------------------------------------------------------
# Nominal vs. observed FDR plot
# This gives a 5x5 grid
# -----------------------------------------------------------------------------
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
                                "Medium" = "blue", 
                                "Low" = "purple")) +
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
  theme_bw(base_size = 11) +
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
# Nominal vs Observed FDR plots
# This gives a 1x3 plot
#-------------------------------------------------------------------------------
all_effects = sort(unique(fdr_eval$Effect))
selected_effects = all_effects[c(1, ceiling(length(all_effects)/2), 
                                 length(all_effects))]

fdr_eval_clean = fdr_eval |>
  filter(
    Res_prop == "Prop_0.50",
    Effect %in% selected_effects
  ) |>
  mutate(
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
  facet_wrap(~ Effect_clean, nrow = 1) +
  scale_color_manual(values = c("High" = "deeppink", 
                                "Medium" = "blue", 
                                "Low" = "purple")) +
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
  theme_bw(base_size = 11) +
  labs(
    title = "False discovery rate.",
    subtitle = "Nominal vs. observed FDR for 50% responders.",
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
# ggsave('_fig/fdr_plot_clean.pdf',plot=fdr_plot_clean,width=10,height=10)
