# MIMOSA2 effect size heterogeneity plots
# Isabella Lethbridge and Tayyeb Abrahams 
# 28 August 2026

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
library(purrr)
library(pROC)

# Load data:
setwd("C:/Github/comp-2026-mimosa2")
load('_simulations/Simulation_Small_Vs_Combined_Continuous.Rdata')

# 1. Reshape data for all 3 models and format Large_effect labels
roc_data_long <- results_continuous %>%
  pivot_longer(
    cols = c(MIMOSA2_prob_independent, MIMOSA2_prob_combined, DiD_GLM_prob),
    names_to = "Model",
    values_to = "Probability"
  ) %>%
  mutate(
    Model_label = factor(
      recode(Model,
             "MIMOSA2_prob_independent" = "Independent",
             "MIMOSA2_prob_combined"    = "Combined",
             "DiD_GLM_prob"             = "DiD GLM"
      ),
      levels = c("Independent", "Combined", "DiD GLM")
    ),
    # Format Large_effect values nicely for panel headers
    Large_effect_clean = factor(
      Large_effect,
      levels = sort(unique(Large_effect)),
      labels = paste0("Effect: ", formatC(sort(unique(Large_effect)), format = "e", digits = 1))
    )
  )

# 2. Compute ROC coordinates and AUCs per panel per model
roc_df <- roc_data_long %>%
  group_by(Large_effect_clean, Model_label) %>%
  group_modify(~ {
    roc_obj <- roc(.x$Truth, .x$Probability, quiet = TRUE)
    data.frame(
      Specificity = roc_obj$specificities,
      Sensitivity = roc_obj$sensitivities,
      AUC = round(as.numeric(auc(roc_obj)), 3)
    )
  }) %>%
  ungroup()

# 3. Build 2x2 Faceted Plot for all 3 models
roc_plot_2x2 <- ggplot(roc_df, aes(
  x = 1 - Specificity, 
  y = Sensitivity, 
  color = Model_label,
  linetype = Model_label
)) +
  geom_path(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  
  # 2x2 Layout by Large Effect Size
  facet_wrap(~ Large_effect_clean, ncol = 2) +
  
  # Color and Linetype mapping for 3 models
  scale_color_manual(values = c(
    "Independent" = "deeppink", 
    "Combined"    = "steelblue3", 
    "DiD GLM"     = "orange"
  )) +
  scale_linetype_manual(values = c(
    "Independent" = "solid", 
    "Combined"    = "dashed", 
    "DiD GLM"     = "dotdash"
  )) +
  
  # Format axes & limits
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.01, 0.01)) +
  labs(
    title = "ROC Curve Comparison Across Effect Sizes",
    subtitle = "Comparing MIMOSA2 (Independent & Combined) vs. DiD GLM",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = "Model",
    linetype = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    legend.position = "bottom"
  )

print(roc_plot_2x2)
# ggsave('_fig/roc_2x2_3models.pdf', plot = roc_plot_2x2, width = 8, height = 7)

#-------------------------------------------------------------------------------
# Load data:
setwd("C:/Github/comp-2026-mimosa2")
load('_simulations/Count_Imbalance.Rdata')

View(results_continuous)
