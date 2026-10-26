## ch5_figures.R — Chapter 5 figures from ch5_results.rds (no refitting)
library(tidyverse)
res <- readRDS("ch5_results.rds")

PAL_ALT     <- "#D55E00"   # vermilion: called by MIMOSA2
PAL_BLUE    <- "#225EA8"   # blue: stimulated samples
PAL_REFLINE <- "#333333"   # dashed reference line
FIG_W       <- 8

theme_mimosa <- function() {
  theme_bw(base_size = 11) +
    theme(axis.title = element_text(size = 14),
          axis.title.y = element_text(margin = margin(r = 8)),
          axis.text = element_text(size = 10),
          strip.text = element_text(face = "bold", size = 11),
          strip.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.3),
          panel.spacing = unit(0.6, "lines"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey80", linewidth = 0.3),
          panel.grid.major.x = element_blank(),
          legend.position = "bottom", legend.box = "horizontal",
          legend.title = element_text(size = 10), legend.text = element_text(size = 10),
          legend.key.width = unit(2, "lines"),
          legend.background = element_blank(), legend.key = element_blank())
}

## The nine populations, in the order of Table 5.1 (plotmath labels)
LABELS <- c("cd4" = 'CD4^"+"', "ccr7~p~" = 'CCR7^"+"', "cd45ra~p~" = 'CD45RA^"+"',
            "klrg1~p~" = 'KLRG1^"+"', "ccr6~p~" = 'CCR6^"+"', "cxcr3~p~" = 'CXCR3^"+"',
            "hladr~p~" = '"HLA-DR"^"+"', "cd103~p~" = 'CD103^"+"', "cd153~p~" = 'CD153^"+"')

subj <- res$subjects |>
  filter(parent_pop %in% names(LABELS),
         pop == ifelse(parent_pop == "cd4", "any_cyt", paste0(parent_pop, "any_cyt"))) |>
  mutate(label    = factor(LABELS[parent_pop], levels = LABELS),
         mim_resp = coalesce(mim_resp_01, FALSE),
         did_resp = did_resp_01,
         call = case_when(mim_resp & did_resp ~ "Both",
                          mim_resp            ~ "MIMOSA2 only",
                          did_resp            ~ "DiD-BH only",
                          TRUE                ~ "Neither"),
         call = factor(call, levels = c("Both", "MIMOSA2 only", "DiD-BH only", "Neither")))

## Panel titles: population and median number of cells per sample
strips <- subj |>
  group_by(label) |>
  summarise(N_med = median(c(Ns1, Nu1, Ns0, Nu0)), .groups = "drop") |>
  mutate(strip = sprintf('%s~"(%s cells)"', as.character(label),
                         format(round(N_med), big.mark = " ", trim = TRUE)))
subj <- subj |>
  left_join(select(strips, label, strip), by = "label") |>
  mutate(strip = factor(strip, levels = strips$strip))

# Check: should match the "Responders called" columns of Table 5.1
print(subj |> group_by(label) |> summarise(P = n(), mimosa2 = sum(mim_resp), did_bh = sum(did_resp)))

## ---- Figure 5.1: parent counts in the four samples ----------------------------
counts_long <- subj |>
  select(pid, strip, Ns1, Nu1, Ns0, Nu0) |>
  pivot_longer(c(Ns1, Nu1, Ns0, Nu0), names_to = "sample", values_to = "N") |>
  mutate(stim  = factor(ifelse(startsWith(sample, "Ns"), "Stimulated", "Unstimulated"),
                        levels = c("Stimulated", "Unstimulated")),
         visit = factor(ifelse(endsWith(sample, "0"), "Day 0", "Day 224"),
                        levels = c("Day 0", "Day 224")))

p_counts <- ggplot(counts_long, aes(visit, N, fill = stim)) +
  geom_boxplot(position = position_dodge(0.75), width = 0.65, outlier.shape = NA,
               colour = "#1A1A1A", linewidth = 0.4, alpha = 0.6) +
  geom_point(aes(group = stim),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
             size = 0.6, colour = "grey30", show.legend = FALSE) +
  scale_fill_manual(values = c(Stimulated = PAL_BLUE, Unstimulated = "white"), name = NULL) +
  scale_y_log10(breaks = as.vector(outer(c(1, 2, 5), 10^(0:6))),   # 1, 2, 5, 10, 20, 50, ...
                labels = scales::label_comma()) +
  facet_wrap(~ strip, scales = "free_y", ncol = 3, labeller = label_parsed) +
  labs(x = NULL, y = "Parent count per sample") +
  theme_mimosa()

print(p_counts)
ggsave("ch5_parent_counts.pdf", p_counts, width = FIG_W, height = 8, device = cairo_pdf)

## ---- Figure 5.2: background-subtracted frequency, d0 -> d224 ------------------
bs_long <- subj |>
  transmute(pid, strip, call,
            bs_0 = 100 * (s0 / Ns0 - u0 / Nu0),
            bs_1 = 100 * (s1 / Ns1 - u1 / Nu1)) |>
  pivot_longer(c(bs_0, bs_1), names_to = "t", values_to = "bs") |>
  mutate(visit = factor(t, levels = c("bs_0", "bs_1"), labels = c("Day 0", "Day 224")))

CALL_COLOURS <- c("Both" = PAL_ALT, "MIMOSA2 only" = PAL_ALT,
                  "DiD-BH only" = "grey20", "Neither" = "grey75")
CALL_LINES   <- c("Both" = "solid", "MIMOSA2 only" = "22",
                  "DiD-BH only" = "solid", "Neither" = "22")
CALL_WIDTHS  <- c("Both" = 0.7, "MIMOSA2 only" = 0.7, "DiD-BH only" = 0.7, "Neither" = 0.4)
DRAW_ORDER   <- c("Neither", "DiD-BH only", "MIMOSA2 only", "Both")   # bottom to top

layers <- lapply(DRAW_ORDER, function(k) list(
  geom_line(data = ~ filter(.x, call == k)),
  geom_point(data = ~ filter(.x, call == k), size = if (k == "Neither") 0.9 else 1.6)))

p_bs <- ggplot(bs_long, aes(visit, bs, group = pid, colour = call,
                            linetype = call, linewidth = call)) +
  geom_hline(yintercept = 0, colour = PAL_REFLINE, linetype = "dashed", linewidth = 0.3) +
  layers +
  scale_colour_manual(values = CALL_COLOURS, drop = FALSE, name = "Responder call") +
  scale_linetype_manual(values = CALL_LINES, drop = FALSE, name = "Responder call") +
  scale_discrete_manual("linewidth", values = CALL_WIDTHS, drop = FALSE, name = "Responder call") +
  scale_x_discrete(expand = expansion(add = 0.2)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4)) +   # about 4 round ticks per panel
  facet_wrap(~ strip, scales = "free_y", ncol = 3, labeller = label_parsed) +
  labs(x = NULL, y = "Background-subtracted frequency (%)") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_mimosa()

print(p_bs)
ggsave("ch5_bs_spaghetti_all.pdf", p_bs, width = FIG_W, height = 9, device = cairo_pdf)
