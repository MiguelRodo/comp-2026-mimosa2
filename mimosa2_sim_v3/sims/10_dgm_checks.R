# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# 10_dgm_checks.R: checks of the data-generating mechanisms
# ==============================================================================
# Run BEFORE the simulations (takes a few minutes, does not need MIMOSA2):
#    cd mimosa2_sim_v3 
#    Rscript sims/10_dgm_checks.R
# Writes CSV tables to _tables/ and stops with an error if a check fails.
# Checks:
#  A. Every generative family has the intended MEAN and VARIANCE at every
#     condition mean (1e6 draws per family x mean). 
#     Also reports what the version-2 parameters gave, for the write-up.
#  B. Realised number of responders for every (P, rho) used, version 2 vs
#     version 3 allocation.
#  C. One very large dataset per family: profile counts, every constraint of
#     every profile holds, truth == (true Delta > 0), and the REALISED effect
#     size among responders vs the nominal delta.
#  D. Reproducibility: the same stream gives the identical dataset; different
#     streams give different datasets.
#  E. (if MIMOSA2 is installed) MIMOSA2 does not use or reset the random-number
#     generator.
# =============================================================================
for (f in c("R/config.R", "R/dgm.R", "R/methods.R", "R/runner.R", "R/scenarios.R")) source(f)
RNGkind("L'Ecuyer-CMRG"); set.seed(1)
fails <- character(0)
check <- function(ok, what) {
  if (!isTRUE(all(ok))) { fails <<- c(fails, what); message("FAIL: ", what) } else message("ok:   ", what)
}

# ----------------- A. marginal mean / variance of each family -----------------
message("\n== A. Marginal moments of each family (1e6 draws) ==")
phi <- 10000
mus <- sort(unique(unlist(lapply(c(1e-3, 2.5e-4, 1.25e-4, 6.25e-5), function(e) condition_means(e)))))
nd <- 1e6
rowsA <- list()
for (fam in c("Beta", "EG", "LN", "SX", "BB")) for (mu in mus) {
  for (cond in c("u", "s")) {
    if (fam != "BB" && cond == "u") next          # only BB differs by condition
    x <- draw_prop(nd, mu, if (cond == "s") "s1" else "u1", fam, phi)
    rowsA[[length(rowsA) + 1]] <- data.frame(
      Family = fam, Draw = if (fam == "BB") cond else "any", mu = mu,
      target_sd = sqrt(beta_var(mu, phi)), mean_ratio = mean(x) / mu,
      sd_ratio = sd(x) / sqrt(beta_var(mu, phi)), skewness = mean((x - mean(x))^3) / sd(x)^3,
      p_below_mu_over_10 = mean(x < mu / 10))
  }
}
# Version-2 parameters for comparison (reported, not used):
for (mu in mus) {
  x <- exp(-rgamma(nd, 113, 113 / (113 * (mu^(-1/113) - 1))))
  rowsA[[length(rowsA) + 1]] <- data.frame(Family = "EG (v2: k = 113)", Draw = "any", mu = mu,
    target_sd = sqrt(beta_var(mu, phi)), mean_ratio = mean(x) / mu, sd_ratio = sd(x) / sqrt(beta_var(mu, phi)),
    skewness = mean((x - mean(x))^3) / sd(x)^3, p_below_mu_over_10 = mean(x < mu / 10))
  x <- plogis(rnorm(nd, qlogis(mu), 1 / sqrt(2.05)))
  rowsA[[length(rowsA) + 1]] <- data.frame(Family = "LN (v2: phi = 2.05, median = mu)", Draw = "any", mu = mu,
    target_sd = sqrt(beta_var(mu, phi)), mean_ratio = mean(x) / mu, sd_ratio = sd(x) / sqrt(beta_var(mu, phi)),
    skewness = mean((x - mean(x))^3) / sd(x)^3, p_below_mu_over_10 = mean(x < mu / 10))
}
tabA <- do.call(rbind, rowsA)
write.csv(tabA, file.path(TAB_DIR, "dgm_check_A_marginal_moments.csv"), row.names = FALSE)
print(format(tabA, digits = 3), row.names = FALSE)
newA <- tabA[!grepl("v2", tabA$Family), ]
check(abs(newA$mean_ratio - 1) < 0.01, "A: every family has mean = mu (within 1%)")
check(abs(newA$sd_ratio - 1) < 0.03, "A: every family has SD = Beta SD (within 3%)")
calA <- data.frame(t(sapply(mus, function(mu) c(mu = mu, EG_k = calib_eg(mu, phi)$k,
                                                   LN_sd = calib_ln(mu, phi)$sd, LN_loc = calib_ln(mu, phi)$loc,
                                                   logit_mu = qlogis(mu)))))
write.csv(calA, file.path(TAB_DIR, "dgm_check_A_calibrated_parameters.csv"), row.names = FALSE)
message("Calibrated parameters (for the write-up):"); print(signif(calA, 4))

# ------------------------- B. realised responders------------------------------
message("\n== B. Realised number of responders (version 2 vs version 3) ==")
allocate_v2 <- function(P, pis) {
  n <- floor(P * pis)
  while (sum(n) > P) { i <- which.max(n); n[i] <- n[i] - 1 }
  while (sum(n) < P) { i <- which.min(n); n[i] <- n[i] + 1 }
  n
}
gridB <- expand.grid(P = c(10, 20, 30, 50, 75, 80, 100), rho = c(0.1, 0.2, 0.25, 0.5, 0.75, 0.8, 0.9))
gridB$rho_v2 <- mapply(function(P, r) sum(allocate_v2(P, resp_components(r))[1:4]) / P, gridB$P, gridB$rho)
gridB$rho_v3 <- mapply(function(P, r) sum(allocate_profiles(P, resp_components(r))[1:4]) / P, gridB$P, gridB$rho)
write.csv(gridB, file.path(TAB_DIR, "dgm_check_B_realised_rho.csv"), row.names = FALSE)
message("Realised rho, version 2 (rows P, columns nominal rho):")
print(xtabs(rho_v2 ~ P + rho, gridB))
message("Realised rho, version 3:")
print(xtabs(rho_v3 ~ P + rho, gridB))
check(abs(gridB$rho_v3 - gridB$rho) <= 0.5 / gridB$P + 1e-9, "B: realised rho within rounding of nominal rho")

# ------------------ C. one very large dataset per family ----------------------
message("\n== C. Large datasets: constraints and realised effect size ==")
constraint_ok <- function(sim) {
  p <- as.data.frame(sim$p); d <- sim$true_delta; tr <- sim$truth
  ok <- rep(NA, length(tr))
  ok[tr == "R1"]  <- with(p[tr == "R1", ], d[tr == "R1"] > 0 & ps1 > pu1)
  ok[tr == "R2"]  <- with(p[tr == "R2", ], ps0 == pu0 & ps1 > pu1)
  ok[tr == "R3"]  <- with(p[tr == "R3", ], ps1 == ps0 & d[tr == "R3"] > 0 & ps1 > pu1 & pu0 > pu1)
  ok[tr == "R4"]  <- with(p[tr == "R4", ], pu1 == pu0 & d[tr == "R4"] > 0 & ps1 > pu1 & ps1 > ps0)
  ok[tr == "NR1"] <- with(p[tr == "NR1", ], ps1 == pu1 & ps0 == pu0)
  ok[tr == "NR2"] <- with(p[tr == "NR2", ], ps1 == pu1 & ps0 >= pu0)
  ok[tr == "NR3"] <- with(p[tr == "NR3", ], ps1 == ps0 & ps0 == pu1 & pu1 == pu0)
  ok[tr == "NSR"] <- with(p[tr == "NSR", ], ps1 == ps0 & pu1 == pu0)
  ok
}
rowsC <- list()
cases <- rbind(
  data.frame(Family = "Beta", Phi = 2000, Effect = c(1e-3, 6.25e-5), Source = "baseline"),
  data.frame(Family = rep(c("Beta", "EG", "LN", "SX", "BB"), each = 2), Phi = 10000,
             Effect = rep(c(1e-3, 6.25e-5), 5), Source = "prior"))
prewarm_calibration(cases$Family, cases$Effect, cases$Phi)
for (i in seq_len(nrow(cases))) {
  cs <- cases[i, ]
  sim <- simulate_MIMOSA2_alt_prior(effect = cs$Effect, phi = cs$Phi, P = 40000, prior = cs$Family,
                                    components = resp_components(0.5), rng = c(1e5, 1e5))
  ok <- constraint_ok(sim)
  resp <- is_responder(sim$truth) == 1
  check(all(ok), sprintf("C: all profile constraints hold (%s, %s, delta = %g)", cs$Source, cs$Family, cs$Effect))
  check(all((sim$true_delta > 0) == resp), sprintf("C: responder <=> true Delta > 0 (%s, %s, delta = %g)", cs$Source, cs$Family, cs$Effect))
  rowsC[[i]] <- data.frame(cs, n_profile = paste(sim$n_profile, collapse = "/"),
                           rho_realised = mean(resp),
                           mean_true_delta_resp = mean(sim$true_delta[resp]),
                           ratio_to_nominal = mean(sim$true_delta[resp]) / cs$Effect,
                           median_true_delta_resp = median(sim$true_delta[resp]),
                           mean_ps1_resp = mean(sim$p[resp, "ps1"]),
                           # realised mean Delta / nominal delta in each responder profile:
                           # R3 has p(s,1) = p(s,0), so its Delta = p(u,0) - p(u,1) does not
                           # depend on delta at all (it comes from the MIMOSA2 simulator design)
                           ratio_R1 = mean(sim$true_delta[sim$truth == "R1"]) / cs$Effect,
                           ratio_R2 = mean(sim$true_delta[sim$truth == "R2"]) / cs$Effect,
                           ratio_R3 = mean(sim$true_delta[sim$truth == "R3"]) / cs$Effect,
                           ratio_R4 = mean(sim$true_delta[sim$truth == "R4"]) / cs$Effect,
                           expected_extra_cells_at_15k = mean(sim$true_delta[resp]) * 15000)
}
tabC <- do.call(rbind, rowsC)
write.csv(tabC, file.path(TAB_DIR, "dgm_check_C_large_datasets.csv"), row.names = FALSE)
print(format(tabC, digits = 3), row.names = FALSE)
message("NOTE for the write-up: ratio_to_nominal != 1 shows that rejection sampling and the profile definitions make the")
message("realised mean Delta among responders differ from the nominal delta (Morris 3.2, pitfall 4):")
message("larger than delta for small delta (truncation), smaller for large delta (profile R3 does not depend on delta).")

# --------------------------- D. reproducibilit---------------------------------
message("\n== D. Reproducibility of streams ==")
des <- build_design("baseline", scenarios_baseline())
row <- des[des$Scenario_ID == des$Scenario_ID[1] & des$Rep == 1, ]
a <- regenerate_dataset(row, simulate_baseline); b <- regenerate_dataset(row, simulate_baseline)
row2 <- des[des$Scenario_ID == des$Scenario_ID[1] & des$Rep == 2, ]
c2 <- regenerate_dataset(row2, simulate_baseline)
check(identical(a, b), "D: same stream -> identical dataset")
check(!identical(a$ns1, c2$ns1), "D: different streams -> different datasets")
des2 <- build_design("baseline", scenarios_baseline())
check(identical(des$Seed, des2$Seed), "D: design (and all streams) rebuilt identically")

# ------------------------- E. MIMOSA2 and the RNG -----------------------------
if (requireNamespace("MIMOSA2", quietly = TRUE)) {
  message("\n== E. Does MIMOSA2 touch the random-number generator? ==")
  sim <- simulate_baseline(data.frame(Effect = 1e-3, Phi = 2000, P = 30, Distribution = "Beta",
                                      Res_prop = "Prop_0.50", Cell_range = "Medium"))
  s0 <- .Random.seed
  invisible(MIMOSA2::MIMOSA2(Ntot = sim$Ntot, ns1 = sim$ns1, nu1 = sim$nu1, ns0 = sim$ns0, nu0 = sim$nu0, maxit = 5, verbose = FALSE))
  s1 <- .Random.seed
  invisible(MIMOSA2::MIMOSA2(Ntot = sim$Ntot, ns1 = sim$ns1, nu1 = sim$nu1, ns0 = sim$ns0, nu0 = sim$nu0, maxit = 5, verbose = FALSE))
  s2 <- .Random.seed
  check(identical(s0, s1) && identical(s1, s2), "E: MIMOSA2 does not use or reset the RNG")
} else message("\n(E skipped: MIMOSA2 not installed here)")

message("\n==============================================================")
if (length(fails)) {
  message("DGM CHECKS FAILED:\n  ", paste(fails, collapse = "\n  "))
  quit(status = 1)
} else message("ALL DGM CHECKS PASSED. Tables written to ", TAB_DIR)
