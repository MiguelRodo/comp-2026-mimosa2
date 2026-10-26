# =============================================================================
# Section 4.6 warning signs: reproduce Table 4.x and test Miguel's p. 36 question
# Run from the folder that holds the four CSVs (your _tables/ folder).
# =============================================================================
library(dplyr)

# Cells per assay for each cell-count label.
# Check these labels match the Cell_range column in your CSVs exactly.
N <- c(High = 250000, Medium = 1e5, Low = 15000,
       "Very Low" = 7000, "Extremely Low" = 3000)
mu0 <- 1e-4   # background proportion used in every scenario

load_study <- function(perf, dgm, study) {
  p <- read.csv(perf) |> filter(Res_prop != "Prop_0.00")   # drop null scenarios
  g <- read.csv(dgm)
  m <- p |> filter(Method == "MIMOSA2", Alpha == 0.01) |>
       select(Scenario_ID, TPR, AUC, Cell_range)
  d <- p |> filter(Method == "DiD", Rule == "BH", Alpha == 0.01) |>
       select(Scenario_ID, TPR_d = TPR)
  m |> inner_join(d, by = "Scenario_ID") |>
       inner_join(select(g, Scenario_ID, mean_true_delta_resp, share_zero_ns1),
                  by = "Scenario_ID") |>
       mutate(study  = study,
              N      = N[Cell_range],
              excess = mean_true_delta_resp * N,   # extra positive cells
              snr    = excess / sqrt(mu0 * N))     # extra cells / sd of background
}

ws <- bind_rows(
  load_study("baseline_performance.csv", "baseline_realised_dgm.csv", "Study 1"),
  load_study("prior_performance.csv",    "prior_realised_dgm.csv",    "Study 2")
)

# --- Sanity checks -----------------------------------------------------------
nrow(ws)            # should be 750 (450 + 300)
sum(is.na(ws$N))    # should be 0; if not, fix the labels in N above

# --- 1. Reproduce the warning-sign table -------------------------------------
ws |>
  mutate(band = cut(excess, c(0, 2, 5, 10, 30, 100, Inf), right = FALSE,
                    labels = c("< 2", "2-5", "5-10", "10-30", "30-100", ">= 100"))) |>
  group_by(band) |>
  summarise(n            = n(),
            sens         = median(TPR),
            q1           = quantile(TPR, .25),
            q3           = quantile(TPR, .75),
            sens_did     = median(TPR_d),
            auc          = median(AUC),
            pct_auc_lt07 = 100 * mean(AUC < 0.7),
            zero_s1      = median(share_zero_ns1)) |>
  print()

# --- 2. Miguel's three candidates: which tracks performance best? ------------
sp <- function(x, y) round(cor(x, y, method = "spearman"), 2)
data.frame(
  quantity    = c("Parent cells (N)",
                  "Extra positive cells (delta x N)",
                  "Extra cells / sqrt(background count)",
                  "Share with zero positive cells in s1"),
  sensitivity = c(sp(ws$N, ws$TPR), sp(ws$excess, ws$TPR),
                  sp(ws$snr, ws$TPR), sp(ws$share_zero_ns1, ws$TPR)),
  AUC         = c(sp(ws$N, ws$AUC), sp(ws$excess, ws$AUC),
                  sp(ws$snr, ws$AUC), sp(ws$share_zero_ns1, ws$AUC))
)
# Earlier session found 0.95 / 0.91 for extra cells and -0.60 for zero counts.

# --- 3. Optional plot: do the cell counts fall on one curve? -----------------
# If the points from different cell counts line up along one curve, the
# x-axis quantity explains performance on its own (Miguel's question).
par(mfrow = c(1, 2))
plot(ws$excess, ws$TPR, log = "x", col = factor(ws$Cell_range), pch = 19,
     xlab = "Extra positive cells (log scale)", ylab = "MIMOSA2 sensitivity")
plot(ws$snr, ws$TPR, log = "x", col = factor(ws$Cell_range), pch = 19,
     xlab = "Extra cells / sqrt(background) (log scale)", ylab = "MIMOSA2 sensitivity")
legend("topleft", legend = levels(factor(ws$Cell_range)),
       col = seq_along(levels(factor(ws$Cell_range))), pch = 19, cex = 0.8)

# --- Table by signal-to-noise bands -------------------------------------------
ws |>
  mutate(band = cut(snr, c(0, 1, 2, 5, 10, 20, Inf), right = FALSE,
                    labels = c("< 1", "1-2", "2-5", "5-10", "10-20", ">= 20"))) |>
  group_by(band) |>
  summarise(n            = n(),
            sens         = round(median(TPR), 2),
            q1           = round(quantile(TPR, .25), 2),
            q3           = round(quantile(TPR, .75), 2),
            sens_did     = round(median(TPR_d), 2),
            auc          = round(median(AUC), 2),
            pct_auc_lt07 = round(100 * mean(AUC < 0.7)),
            zero_s1      = round(median(share_zero_ns1), 2)) |>
  print()
