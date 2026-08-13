library(ggplot2)
library(cowplot)
library(Hmisc) # Required for mean_cl_boot

setwd("GitHub/comp-2026-mimosa2")
setwd("C:/Github/comp-2026-mimosa2") # for Bella 
load('_simulations/Simulation_2.0.Rdata')

# 1. Filter for successful runs (DO NOT aggregate yet)
plot_data = subset(results_summary, Status == "Success")

# 2. Clean up the Effect Factor labels safely on raw data
unique_effects = sort(unique(plot_data$Effect), decreasing = TRUE)
plot_data$Effect_fact = factor(plot_data$Effect, 
                               levels = unique_effects,
                               labels = paste0("E: ", formatC(unique_effects, format = "e", digits = 1)))

# 3. Clean up the Proportion labels
plot_data$Res_prop_clean = factor(plot_data$Res_prop,
                                  levels = c("Prop_0.10", "Prop_0.25", "Prop_0.50", "Prop_0.75", "Prop_0.90"),
                                  labels = c("10% Resp", "25% Resp", "50% Resp", "75% Resp", "90% Resp"))

# 4. Clean, Professional Line Plot with Bootstrap CIs
sim2.0_plot = ggplot(data = plot_data,
                     mapping = aes(x = Effect_fact,
                                   y = TPR_001,
                                   group = Cell_range, 
                                   color = Cell_range,
                                   fill = Cell_range)) +
  
  # Option A: Shaded Confidence Ribbon (Recommended for clean trend lines)
  stat_summary(fun.data = mean_cl_boot, 
               geom = "ribbon", 
               alpha = 0.5, 
               color = NA,
               fun.args = list(B = 2000, conf.int = 0.95)) +
  
  # Mean line tracking mean across bootstrap samples
  stat_summary(fun = mean, 
               geom = "line", 
               linewidth = 0.5) +
  
  # Mean points
  stat_summary(fun = mean, 
               geom = "point", 
               size = 2) +
  
  # Clean faceting layout with fixed labels
  facet_grid(Res_prop_clean ~ P, labeller = label_both) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  
  # Styling
  cowplot::theme_cowplot(font_size = 11) +
  cowplot::background_grid(major = "xy") +
  scale_color_brewer(palette = "Set1") + 
  scale_fill_brewer(palette = "Set1") + 
  labs(
    title = "MIMOSA2 Sensitivity Analysis (Mean & 95% Bootstrap CI)",
    subtitle = "True Positive Rate (TPR) at 1% Nominal FDR Threshold",
    x = "Simulated Effect Size",
    y = "Mean Sensitivity (TPR)",
    color = "Cell Count Range",
    fill = "Cell Count Range"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text.y = element_text(angle = 0, face = "bold", size = 8),
    strip.text.x = element_text(face = "bold", size = 9),
    panel.spacing = unit(0.4, "lines"),
    legend.position = "bottom"
  )

print(sim2.0_plot)

#==============================================================
#DiD comparison
#==============================================================

library(tidyverse)
library(plotROC)

# 1. Load your dataset


# 2. Prepare data for plotting (drop rows where MIMOSA failed/returned NaN)
ROC_data_prepared <- results_continuous %>%
  filter(!is.na(MIMOSA2_prob), !is.na(DiD_GLM_prob)) %>%
  pivot_longer(
    cols      = c(MIMOSA2_prob, DiD_GLM_prob),
    names_to  = "Method",
    values_to = "Score"
  ) %>%
  mutate(
    Method = case_when(
      Method == "MIMOSA2_prob" ~ "MIMOSA2",
      Method == "DiD_GLM_prob"  ~ "DiD Baseline"
    ),
    # Create clean factor labels for plotting
    Sample_Size = paste0("N: ", P),
    Effect_Label = paste0("Effect: ", Effect)
  )

# 3. Generate the ROC Plot grid 
# (You can swap the facet variables depending on which slice you want to look at)

did_cutoff_points <- results_continuous %>%
  filter(!is.na(DiD_GLM_prob), !is.na(Truth)) %>%
  group_by(Res_prop, Effect, Cell_range) %>%
  summarise(
    # Rejection decision at p <= 0.01
    pred_positive = DiD_GLM_prob >= 0.99,
    
    # Calculate rates
    FPR = mean(pred_positive[Truth == 0], na.rm = TRUE),
    TPR = mean(pred_positive[Truth == 1], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Method = "DiD Baseline" # Matches your linetype label
  )

ggplot(data = ROC_data_prepared,
       mapping = aes(d = Truth,
                     m = Score,
                     colour = Cell_range,
                     linetype = Method,
                     group = interaction(Method, Cell_range))) + 
  
  geom_roc(n.cuts = 0, size = 0.9) +
  
  # Reference line
  geom_abline(slope = 1, intercept = 0, linetype = 'dotted', colour = 'grey40', linewidth = 0.6) +
  
  # --- ADD THIS LAYER FOR THE THRESHOLD POINTS ---
  geom_point(
    data = did_cutoff_points,
    aes(x = FPR, y = TPR, fill = Cell_range),
    inherit.aes = FALSE, # Prevent inheritance of 'd' and 'm' from geom_roc
    shape = 21,
    color = "black",
    size = 2.5,
    stroke = 0.8,
    show.legend = FALSE
  ) +
  # -----------------------------------------------

facet_grid(Res_prop ~ Effect, labeller = label_both) + 
  
  scale_colour_manual(
    values = c(
      "Sparse" = "magenta",
      "Medium" = "violet",
      "High"   = "purple"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Sparse" = "magenta",
      "Medium" = "violet",
      "High"   = "purple"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "MIMOSA2"      = "solid",
      "DiD Baseline" = "dashed"
    )
  ) +
  
  theme_bw() +
  labs(
    title    = 'ROC Performance across Simulation Parameters',
    x        = 'False Positive Rate (1 - Specificity)',
    y        = 'True Positive Rate (Sensitivity)',
    colour   = 'Cell Count Range',
    linetype = 'Model Framework'
  ) +
  theme(
    legend.position  = 'bottom',
    legend.box       = 'vertical', 
    plot.title       = element_text(face = 'bold', hjust = 0.5),
    strip.text.x     = element_text(size = 9, face = "bold"),
    strip.text.y     = element_text(size = 9, face = "bold", angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey95"),
    plot.margin      = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")
  )

# 4. Calculate AUROC values dynamically
AUROC <- ROC_data_prepared %>%
  group_by(Res_prop, P, Cell_range, Effect, Method) %>%
  filter(length(unique(Truth)) == 2) %>% 
  do(plotROC::calc_auc(
    ggplot(., aes(d = Truth, m = Score)) + geom_roc()
  )) %>%
  ungroup() %>%
  rename(AUROC = AUC)

write.table(AUROC,"auc.text")

print(head(AUROC))

print(ROC_plot)

View(AUROC[AUROC$Cell_range=='High',])
