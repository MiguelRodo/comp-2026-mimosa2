# MIMOSA2 presentation plots
# Isabella Lethbridge and Tayyeb Abrahams 
# 28 August 2026

#-------------------------------------------------------------------------------
# Load require libraries: 
library(ggplot2)
library(dplyr)

#-------------------------------------------------------------------------------
# MIMOSA2 Model Framework 
#-------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(ggpattern) # Install via install.packages("ggpattern") if needed

# 1. Update params to include pattern types
params = data.frame(
  Condition = factor(c("u0","s0","u1","s1"), levels = c("u0","s0","u1","s1")),
  Label     = c("unstimulated t0", "stimulated t0", "unstimulated t1", "stimulated t1"),
  mean      = c(0.0010, 0.0015, 0.0055, 0.0060), 
  sd        = c(0.0006, 0.0006, 0.0006, 0.0006)
)

x_grid = seq(-0.001, 0.009, length.out = 1000)

plot_data = params |>
  group_by(Condition, Label) |>
  reframe(
    x = x_grid,
    density = dnorm(x_grid, mean = mean, sd = sd)
  )

# 2. Plot with ggpattern
vaccine_effect_plot = ggplot(data = plot_data, 
                             mapping = aes(x = x, 
                                           y = density, 
                                           color = Condition, 
                                           fill = Condition,
                                           alpha = Condition,
                                           pattern = Condition,
                                           pattern_fill = Condition,
                                           pattern_color = Condition)) +
  geom_area_pattern(
    position = "identity",
    pattern_density = 0.03,    # Adjust density of the hatching lines
    pattern_spacing = 0.05,   # Adjust distance between lines
    pattern_scale = 0.5
  ) +
  geom_line(linewidth = 0.9) +
  geom_vline(
    data = params, 
    aes(xintercept = mean, color = Condition), 
    linetype = "dashed", 
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Annotations
  annotate("segment", x = 0.0010, xend = 0.0015, y = 720, yend = 720,
           arrow = arrow(ends = "both", length = unit(0.12, "cm")), color = "grey10") +
  annotate("text", x = 0.00125, y = 750, label = "(s0 - u0)", size = 3.0, color = "grey10", fontface = "bold") +
  annotate("segment", x = 0.0055, xend = 0.0060, y = 720, yend = 720,
           arrow = arrow(ends = "both", length = unit(0.12, "cm")), color = "grey10") +
  annotate("text", x = 0.00575, y = 750, label = "(s1 - u1)", size = 3.0, color = "grey10", fontface = "bold") +
  annotate("curve", x = 0.00125, xend = 0.00575, y = 820, yend = 820, curvature = -0.15,
           arrow = arrow(ends = "both", length = unit(0.12, "cm")), linetype = "dashed", color = "grey10") +
  annotate("text", x = 0.00350, y = 910, label = "Vaccine effect: (s1 - u1) - (s0 - u0)", 
           size = 4, color = "grey10", fontface = "bold") +
  
  # Cross-hatch for t0 ('crosshatch'), Solid for t1 ('none')
  scale_pattern_manual(
    values = c("u0" = "stripe", "s0" = "stripe", "u1" = "none", "s1" = "none"),
    labels = params$Label
  ) +
  scale_alpha_manual(
    values = c("u0" = 0.1, "s0" = 0.1, "u1" = 0.4, "s1" = 0.4),
    labels = params$Label
  ) +
  # Match hatch line color to the group outlines
  scale_pattern_color_manual(
    values = c("u0" = "orange", "s0" = "deeppink", "u1" = "orange", "s1" = "deeppink"),
    labels = params$Label
  ) +
  scale_color_manual(
    values = c("u0" = "orange", "s0" = "deeppink", "u1" = "orange", "s1" = "deeppink"),
    labels = params$Label
  ) +
  scale_fill_manual(
    values = c("u0" = "orange", "s0" = "deeppink", "u1" = "orange", "s1" = "deeppink"),
    labels = params$Label
  ) +
  scale_pattern_fill_manual(
    values = c("u0" = "orange", "s0" = "deeppink", "u1" = "orange", "s1" = "deeppink"),
    labels = params$Label
  ) +
  scale_x_continuous(limits = c(-0.001,0.008),
                     breaks = seq(0, 0.008, 0.002),
                     labels = sprintf("%.3f", seq(0, 0.008, 0.002))) +
  scale_y_continuous(limits = c(0, 980), expand = c(0, 0)) + 
  theme_minimal(base_size = 14) +
  labs(
    title = "MIMOSA2 model framework.",
    x = "Functional cell proportions",
    y = ""
  ) +
  theme_bw(base_size=14)+
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8.5),
    panel.background = element_rect(fill='grey98'),
    panel.grid.minor = element_line(linewidth=0.5),
    panel.grid.major = element_line(linewidth = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    # Remove axis text (numbers/labels)
    axis.text = element_blank(),
    axis.text.x.bottom = element_blank(),
    
    # Remove axis tick marks
    axis.ticks = element_blank(),
    # Remove the full bounding box around the panel
    panel.border = element_blank(),
    
    # Add back only the bottom x-axis edge line
    axis.line.x = element_line(color = "black", linewidth = 0.5)
)
print(vaccine_effect_plot)

ggsave('MIMOSA2_fw.pdf',plot = vaccine_effect_plot,width=8,height=6)
