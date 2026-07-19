library(MIMOSA2)
library(ggplot2)

setwd("~/GitHub/comp-2026-mimosa2/")
load('_simulations/Simulation_2.0.Rdata')

library(ggplot2)

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
                                        levels = c("Prop_0.00", "Prop_0.25", "Prop_0.50", "Prop_0.75", "Prop_1.00"),
                                        labels = c("0% Resp", "25% Resp", "50% Resp", "75% Resp", "100% Resp"))

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
  theme_bw(base_size = 11) +
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
