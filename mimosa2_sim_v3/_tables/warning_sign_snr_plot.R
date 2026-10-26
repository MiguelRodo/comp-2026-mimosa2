# =============================================================================
# Section 4.6 figure: sensitivity against extra positive cells (A) and against
# extra positive cells / sqrt(background count) (B), coloured by cell count.
# Thesis palette v2 and global style (analysis_functions.R).
# Run after warning_signs_check.R, so that `ws` exists.
# =============================================================================
library(ggplot2)
library(tidyr)
library(dplyr)

# --- Palette v2 (same hex codes as analysis_functions.R) ----------------------
# Names match the Cell_range labels in your CSVs ("Very Low", "Extremely Low").
CELL_COLOURS_WS <- c("High"          = "#081D58",
                     "Medium"        = "#253494",
                     "Low"           = "#225EA8",
                     "Very Low"      = "#1D91C0",
                     "Extremely Low" = "#41B6C4")
cell_labels <- c("High" = "High", "Medium" = "Medium", "Low" = "Low",
                 "Very Low" = "Very low", "Extremely Low" = "Extremely low")

panel_labels <- c(
  excess = "Extra positive cells",
  snr    = "Extra positive cells / √(background count)"
)

plot_df <- ws |>
  select(Cell_range, TPR, excess, snr) |>
  pivot_longer(c(excess, snr), names_to = "measure", values_to = "x") |>
  mutate(measure    = factor(measure, levels = names(panel_labels),
                             labels = panel_labels),
         Cell_range = factor(Cell_range, levels = names(CELL_COLOURS_WS)))

p <- ggplot(plot_df, aes(x = x, y = TPR, colour = Cell_range)) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~ measure, nrow = 1, scales = "free_x",
             strip.position = "bottom") +          # strips act as x-axis titles
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100, 200),
                labels = c("0.5", "1", "2", "5", "10", "20", "50", "100", "200")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = CELL_COLOURS_WS, labels = cell_labels,
                      name = "Cell count") +
  labs(x = NULL, y = "MIMOSA2 sensitivity") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 1))) +
  theme_bw(base_size = 10) +
  theme(axis.title        = element_text(size = 14),
        axis.text         = element_text(size = 10),
        strip.text        = element_text(size = 11),
        strip.background  = element_blank(),
        strip.placement   = "outside",
        legend.position   = "bottom",
        legend.title      = element_text(size = 10),
        legend.text       = element_text(size = 10),
        panel.border      = element_rect(colour = "#4D4D4D", fill = NA),
        panel.grid.major  = element_line(colour = "#EBEBEB"),
        panel.grid.minor  = element_blank())

p

# Same size and device as the other thesis figures (8 in wide, fonts embedded).
ggsave("warning_sign_snr.pdf", p, width = 8, height = 4.2, device = cairo_pdf)
