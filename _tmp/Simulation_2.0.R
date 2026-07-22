library(MIMOSA2)
library(ggplot2)

setwd("~/GitHub/comp-2026-mimosa2/")
load('_simulations/Simulation_2.0.Rdata')

# 1. Filter for successful runs
plot_data = subset(results_summary, Status == "Success")

# 2. Aggregate across the 5 replications (Calculate mean TPR per scenario)
plot_aggregated = aggregate(TPR_001 ~ Res_prop + P + Effect + Cell_range, 
                            data = plot_data, 
                            FUN = mean)

# 3. Clean up the Effect Factor labels safely using actual unique values
unique_effects = sort(unique(plot_aggregated$Effect), decreasing = TRUE)
plot_aggregated$Effect_fact = factor(plot_aggregated$Effect, 
                                     levels = unique_effects,
                                     labels = paste0("E: ", formatC(unique_effects, format = "e", digits = 1)))

# 4. Clean up the Proportion labels so they fit nicely on the facet strips
plot_aggregated$Res_prop_clean = factor(plot_aggregated$Res_prop,
                                        levels = c("Prop_0.10", "Prop_0.25", "Prop_0.50", "Prop_0.75", "Prop_0.90"),
                                        labels = c("10% Resp", "25% Resp", "50% Resp", "75% Resp", "90% Resp"))

# 5. Clean, Professional Line Plot
sim2.0_plot = ggplot(data = plot_aggregated,
                     mapping = aes(x = Effect_fact,
                                   y = TPR_001,
                                   group = Cell_range, 
                                   color = Cell_range)) +
  # Single clean line tracking the average across replicates
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  
  # Clean faceting layout with fixed labels
  facet_grid(Res_prop_clean ~ P, labeller = label_both) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  
  # Styling
  # theme_bw(base_size = 11) +
  cowplot::theme_cowplot(font_size = 11) +
  cowplot::background_grid(major = "xy") +
  scale_color_brewer(palette = "Set1") + 
  labs(
    title = "MIMOSA2 Sensitivity Analysis (Mean Across Replicates)",
    subtitle = "True Positive Rate (TPR) at 1% Nominal FDR Threshold",
    x = "Simulated Effect Size",
    y = "Mean Sensitivity (TPR)",
    color = "Cell Count Range"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text.y = element_text(angle = 0, face = "bold", size = 8), # Flips right labels to read horizontally
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

ggplot(data = ROC_data_prepared,
       mapping = aes(d = Truth,
                     m = Score,
                     colour = Cell_range,               # 🎨 Color represents Cell Count range
                     linetype = Method,                 # ── Linetype represents the Model type
                     group = interaction(Method, Cell_range))) + 
  geom_roc(n.cuts = 0, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'grey50') +
  
  facet_grid(Res_prop ~ Effect, labeller = label_both) + 
  
  # 🌟 high-contrast color palette 🌟
  scale_colour_manual(
    values = c(
      "Wide_High"  = "#2c7bb6",  # Clean Deep Blue
      "Medium_Low" = "#abd9e9",  # Light Ice Blue
      "Sparse"     = "#fdae61",  # Vibrant Orange
      "V_Sparse"   = "#d7191c"   # High-visibility Crimson Red
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
    legend.position = 'bottom',
    legend.box      = 'vertical', 
    plot.title      = element_text(face = 'bold', hjust = 0.5),
    
    # 🌟 FIXES FOR THE RIGHT STRIP LABELS 🌟
    strip.text.x    = element_text(size = 9, face = "bold"), # Keeps top labels clean
    strip.text.y    = element_text(size = 9, face = "bold", angle = 0, hjust = 0), # Un-rotates the right labels
    strip.background = element_rect(fill = "grey95"),
    
    # Adds a small cushion on the right margin of the entire plot canvas so nothing clips
    plot.margin     = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt") 
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

print(head(AUROC))

print(ROC_plot)
