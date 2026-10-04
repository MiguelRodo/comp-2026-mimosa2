# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# scenarios.R: DGM (scenario) of every study 
# ==============================================================================
# All grids are now defined once, here, and used by the simulation, check and
# analysis scripts. 
# Column In_standard marks the scenarios run in the standard (~20 core) profile
# The extended profile runs every row.
# In_smoke marks a handful of scenarios for the smoke test.
# The factor levels are those of version 2 ("fix errors, keep designs").
# Labels: the 15,000-cell level is called "Low" everywhere 
# =============================================================================
# ------------------------------ Shared lists ----------------------------------
resp_components <- function(rho) c(rep(rho / 4, 4), rep((1 - rho) / 4, 4))

CELL_LEVELS <- c("High", "Medium", "Low", "Very Low", "Extremely Low")
CELL_COUNTS <- c("High" = 250000, "Medium" = 100000, "Low" = 15000,
                 "Very Low" = 7000, "Extremely Low" = 3000)
cell_rng <- function(level) rep(CELL_COUNTS[[level]], 2)

RES_LEVELS <- c("Prop_0.00", "Prop_0.10", "Prop_0.20", "Prop_0.25", "Prop_0.50",
                "Prop_0.75", "Prop_0.80", "Prop_0.90")
res_rho <- function(lbl) as.numeric(sub("Prop_", "", lbl))

# =============================================================================
# Study 1: BASELINE (Beta DGM, MIMOSA2's own assumptions)
# =============================================================================
# As Simulation.R: 
#    5 responder proportions x 6 P x 5 effects x 3 cell counts
#    phi = 2000. 
# Added a NULL scenario with no responders (rho = 0) for every P x cell count 
# (every call there is a false positive, so it checks FDR control of MIMOSA2 and the type I error of the DiD test directly. 
# (delta = 0 would NOT give a null: the rejection sampling still forces Delta > 0 for profiles R1-R4.)
scenarios_baseline <- function() {
  main <- expand.grid(Res_prop = c("Prop_0.10", "Prop_0.25", "Prop_0.50", "Prop_0.75", "Prop_0.90"),
                      P = c(10, 20, 30, 50, 75, 100),
                      Effect = c(1e-3, 5e-4, 2.5e-4, 1.25e-4, 6.25e-5),
                      Cell_range = c("High", "Medium", "Low"),
                      stringsAsFactors = FALSE)
  null <- expand.grid(Res_prop = "Prop_0.00",
                      P = c(10, 20, 30, 50, 75, 100),
                      Effect = 0,
                      Cell_range = c("High", "Medium", "Low"),
                      stringsAsFactors = FALSE)
  s <- rbind(main, null)
  s$Distribution <- "Beta"; s$Phi <- 2000
  s$Scenario_ID <- sprintf("BL%03d", seq_len(nrow(s)))
  s$In_standard <- TRUE                 # baseline is fully factorial in both profiles
  s$In_smoke <- s$P == 20 & s$Cell_range == "Medium" &
    ((s$Res_prop == "Prop_0.50" & s$Effect %in% c(1e-3, 2.5e-4)) | s$Res_prop == "Prop_0.00")
  s
}
simulate_baseline <- function(row) {
  simulate_MIMOSA2_alt_prior(effect = row$Effect, phi = row$Phi, P = row$P,
                             prior = row$Distribution,
                             components = resp_components(res_rho(row$Res_prop)),
                             rng = cell_rng(row$Cell_range))
}

# =============================================================================
# Study 2: PRIOR FAMILY (robustness to the Beta assumption)
# =============================================================================
# As Prior_simulations.R: 
#    5 families x 3 P x 4 effects x 5 cell counts
#    rho = 0.5
#    phi = 10000 is now the precision of the Beta REFERENCE, and every
# family is calibrated to the Beta mean and variance at each condition mean
# (dgm.R, CHANGE 2). 
# Standard profile: P = 50 only
scenarios_prior <- function() {
  s <- expand.grid(Distribution = c("Beta", "EG", "LN", "SX", "BB"),
                   P = c(20, 50, 100),
                   Effect = c(1e-3, 2.5e-4, 1.25e-4, 6.25e-5),
                   Cell_range = CELL_LEVELS,
                   stringsAsFactors = FALSE)
  s$Res_prop <- "Prop_0.50"; s$Phi <- 10000
  s$Scenario_ID <- sprintf("PF%03d", seq_len(nrow(s)))
  s$In_standard <- s$P == 50
  s$In_smoke <- s$P == 50 & s$Effect == 1e-3 & s$Cell_range == "Medium" & s$Distribution %in% c("Beta", "LN", "BB")
  s
}
simulate_prior <- function(row) {
  simulate_MIMOSA2_alt_prior(effect = row$Effect, phi = row$Phi, P = row$P,
                             prior = row$Distribution,
                             components = resp_components(res_rho(row$Res_prop)),
                             rng = cell_rng(row$Cell_range))
}

# =============================================================================
# Study 3: EFFECT-SIZE HETEROGENEITY
# =============================================================================
# Each dataset = a FOCAL subgroup (P_focal subjects, small effect delta_focal)
# + an OTHER subgroup (P_other subjects, larger effect delta_other).
# Performance is measured on the focal (weak-effect) subjects, for MIMOSA2
# fitted to all subjects together (pooled) vs to the focal subjects alone
# (separate), and for DiD.
#
# In Heterogeneous_effect_sim.R the effect pairs are c(smaller, larger), but the 
# code set eff_large <- pair[1] (the SMALLER one) and eff_small <- pair[2] (the 
# LARGER one). 
# The subgroup that was evaluated (rows 1:p_small, labelled "Small_effect") 
# therefore had the LARGER effect, the opposite of the design in Chapter 3. 
# Here the names are explicit:
# Effect_focal < Effect_other, and the focal group is evaluated.
# Grid: 
#    rho in {0.2, 0.5, 0.8}
#    (P_focal, P_other) in {(20,80), (50,50), (80,20)}
#    delta_focal in {5e-4, 6.25e-4}
#    delta_other in {6.25e-4, 8e-4, 1e-3, 1e-2, 5e-2} with delta_other > delta_focal (9 pairs)
#    5 cell counts
#    phi = 5000.
# Standard profile: rho = 0.5, delta_focal = 5e-4, cell counts Medium and Low.
scenarios_heterogeneity <- function() {
  pairs <- expand.grid(Effect_focal = c(5e-4, 6.25e-4),
                       Effect_other = c(6.25e-4, 8e-4, 1e-3, 1e-2, 5e-2))
  pairs <- pairs[pairs$Effect_other > pairs$Effect_focal, ]
  splits <- data.frame(P_focal = c(20, 50, 80), P_other = c(80, 50, 20))
  s <- merge(merge(pairs, splits, by = NULL),
             expand.grid(Res_prop = c("Prop_0.20", "Prop_0.50", "Prop_0.80"),
                         Cell_range = CELL_LEVELS, stringsAsFactors = FALSE), by = NULL)
  s <- s[order(s$Res_prop, s$Cell_range, s$P_focal, s$Effect_focal, s$Effect_other), ]
  rownames(s) <- NULL
  s$Distribution <- "Beta"; s$Phi <- 5000
  s$Scenario_ID <- sprintf("EH%03d", seq_len(nrow(s)))
  s$In_standard <- s$Res_prop == "Prop_0.50" & s$Effect_focal == 5e-4 & s$Cell_range %in% c("Medium", "Low")
  s$In_smoke <- s$In_standard & s$P_focal == 20 & s$Cell_range == "Medium" & s$Effect_other %in% c(1e-3, 5e-2)
  s
}
simulate_heterogeneity <- function(row) {
  comp <- resp_components(res_rho(row$Res_prop))
  focal <- simulate_MIMOSA2_alt_prior(effect = row$Effect_focal, phi = row$Phi, P = row$P_focal,
                                      prior = row$Distribution, components = comp,
                                      rng = cell_rng(row$Cell_range))
  other <- simulate_MIMOSA2_alt_prior(effect = row$Effect_other, phi = row$Phi, P = row$P_other,
                                      prior = row$Distribution, components = comp,
                                      rng = cell_rng(row$Cell_range))
  bind_sims(focal, other)
}

# =============================================================================
# Study 4: CELL-COUNT IMBALANCE
# =============================================================================
# As Count_Imbalanced_Sim.R: 
#    Reference 150,000 cells in every assay
#    one, two or three assays depleted to a fraction f of the reference 
#    ("Imbalanced") and "Balanced" controls with the same TOTAL spread equally 
#    over the four assays, (4 - k + k f) / 4 * 150,000 per assay for k depleted 
#    assays.
#    delta = 5e-4
#    phi = 5000
#    P = 100
#    rho = 0.25
# The 120 scenarios are now generated from (pattern, f) instead of a
# hand-typed list, and the depleted assays are stored explicitly.
# Standard profile: reference, Balanced k = 1, 2, the four single-assay
# patterns and the two "same time point" pairs (s1+u1, s0+u0).
IMB_FRACTIONS <- c(0.90, 0.75, 0.50, 0.25, 0.10, 0.05, 0.025)
IMB_PATTERNS <- c("s1", "u1", "s0", "u0",
                  "s1_s0", "u1_u0", "s1_u1", "s0_u0", "s1_u0", "s0_u1",
                  "s1_u1_u0", "s1_u1_s0", "s1_u0_s0", "u1_u0_s0")
frac_label <- function(f) c("0.9" = ".90", "0.75" = ".75", "0.5" = ".50", "0.25" = ".25",
                            "0.1" = ".10", "0.05" = ".050", "0.025" = ".025", "1" = "1.00")[as.character(f)]
scenarios_imbalance <- function() {
  ref <- data.frame(Type = "Reference", Fraction = 1, Depleted = "none", N_depleted = 0,
                    Cell_range = "1.00_Count_Ref", stringsAsFactors = FALSE)
  bal <- expand.grid(N_depleted = 1:3, Fraction = IMB_FRACTIONS)
  bal <- data.frame(Type = "Balanced", Fraction = bal$Fraction, Depleted = "none",
                    N_depleted = bal$N_depleted,
                    Cell_range = paste0(frac_label(bal$Fraction), "_Count_Bal_", bal$N_depleted),
                    stringsAsFactors = FALSE)
  imb <- expand.grid(Depleted = IMB_PATTERNS, Fraction = IMB_FRACTIONS, stringsAsFactors = FALSE)
  imb <- data.frame(Type = "Imbalanced", Fraction = imb$Fraction, Depleted = imb$Depleted,
                    N_depleted = lengths(strsplit(imb$Depleted, "_")),
                    Cell_range = paste0(frac_label(imb$Fraction), "_Count_Imb_", imb$Depleted),
                    stringsAsFactors = FALSE)
  s <- rbind(ref, bal, imb)
  s$Distribution <- "Beta"; s$Phi <- 5000; s$P <- 100; s$Effect <- 5e-4; s$Res_prop <- "Prop_0.25"
  s$Scenario_ID <- sprintf("CI%03d", seq_len(nrow(s)))
  s$In_standard <- s$Type == "Reference" |
    (s$Type == "Balanced" & s$N_depleted %in% 1:2) |
    (s$Type == "Imbalanced" & s$Depleted %in% c("s1", "u1", "s0", "u0", "s1_u1", "s0_u0"))
  s$In_smoke <- s$Type == "Reference" | (s$Fraction == 0.10 & (s$Depleted == "s1" | (s$Type == "Balanced" & s$N_depleted == 1)))
  s
}
# per-assay totals in the column order used by the simulator: nu1, ns1, nu0, ns0
imbalance_counts <- function(row, N = 150000) {
  cnt <- c(u1 = N, s1 = N, u0 = N, s0 = N)
  if (row$Type == "Balanced") {
    cnt[] <- (4 - row$N_depleted + row$N_depleted * row$Fraction) / 4 * N
  } else if (row$Type == "Imbalanced") {
    dep <- strsplit(row$Depleted, "_")[[1]]
    stopifnot(all(dep %in% names(cnt)))
    cnt[dep] <- row$Fraction * N
  }
  cnt
}
simulate_imbalance <- function(row) {
  cnt <- imbalance_counts(row)
  rng <- as.vector(rbind(cnt, cnt))      # (min, max) pairs: u1,u1, s1,s1, u0,u0, s0,s0
  simulate_MIMOSA2_alt_prior(effect = row$Effect, phi = row$Phi, P = row$P,
                             prior = row$Distribution,
                             components = resp_components(res_rho(row$Res_prop)),
                             rng = rng)
}

# =============================================================================
# Study 5: is maxit = 30 enough?
# =============================================================================
# Re-generates EXACTLY the same datasets as the baseline study (same streams)
# for 12 baseline scenarios and fits MIMOSA2 with maxit = 30 and maxit = 100.
# Justifies 'maxit=30' and shows that the stored random-number states reproduce 
# datasets exactly (the maxit = 30 results must equal the baseline results for the 
# same Task_ID).
maxit_check_subset <- function(scen_baseline) {
  with(scen_baseline, Res_prop == "Prop_0.50" & P %in% c(20, 50, 100) &
         Effect %in% c(2.5e-4, 6.25e-5) & Cell_range %in% c("Medium", "Low"))
}
