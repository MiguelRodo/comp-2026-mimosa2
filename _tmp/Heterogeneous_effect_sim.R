
library(R.utils)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. PARAMETERS & SIMULATION MATRIX SETUP
# -----------------------------------------------------------------------------
responders50 = c(rep(0.50 / 4, 4), rep(0.50 / 4, 4))

component_list = list(
  "Prop_0.50" = responders50
)

rng_list = list(
  "Medium"   = c(100000, 100000)
)

stresstest_mat = expand.grid(
  Distribution   = "Beta",
  Comp_Name      = names(component_list),
  P_small        = 20,
  P_big          = 80,
  Large_effect   = c(5e-2,1e-3,8e-4,6.25e-4),
  Small_effect   = c(5e-4),
  Rng_Name       = names(rng_list),
  Phi            = 5000,
  Replication    = 1:30, 
  stringsAsFactors = FALSE
)

# Containers for observation-level and summary results
obs_results_list     <- vector("list", nrow(stresstest_mat))
summary_results_list <- vector("list", nrow(stresstest_mat))

# -----------------------------------------------------------------------------
# 2. MAIN FOR LOOP
# -----------------------------------------------------------------------------
message("Starting sequential simulation loop...")

for (i in 1:nrow(stresstest_mat)) {
  
  # Extract parameters for row i
  dist       <- stresstest_mat$Distribution[i]
  comp_nm    <- stresstest_mat$Comp_Name[i]
  p_small    <- stresstest_mat$P_small[i]
  p_big      <- stresstest_mat$P_big[i]
  eff_large  <- stresstest_mat$Large_effect[i]
  eff_small  <- stresstest_mat$Small_effect[i]
  rng_nm     <- stresstest_mat$Rng_Name[i]
  phi        <- stresstest_mat$Phi[i]
  rep_id     <- stresstest_mat$Replication[i]
  
  active_components <- component_list[[comp_nm]]
  active_rng        <- rng_list[[rng_nm]]
  
  cat(sprintf("Run %d/%d: Rep=%d, CellRange=%s, SmallEff=%.2e, LargeEff=%.2e\n",
              i, nrow(stresstest_mat), rep_id, rng_nm, eff_small, eff_large))
  
  # --- Step A: Generate Small and Big Simulations ---
  sim_small <- simulate_MIMOSA2_alt_prior(
    effect     = eff_small,
    phi        = phi,
    P          = p_small,
    prior      = dist,
    components = active_components,
    rng        = active_rng
  )
  
  sim_big <- simulate_MIMOSA2_alt_prior(
    effect     = eff_large,
    phi        = phi,
    P          = p_big,
    prior      = dist,
    components = active_components,
    rng        = active_rng
  )
  
  # --- Step B: Combine Dataset (sim_all) ---
  sim_all <- list(
    Ntot  = rbind(sim_small$Ntot, sim_big$Ntot),
    ns1   = c(sim_small$ns1, sim_big$ns1),
    nu1   = c(sim_small$nu1, sim_big$nu1),
    ns0   = c(sim_small$ns0, sim_big$ns0),
    nu0   = c(sim_small$nu0, sim_big$nu0),
    truth = c(sim_small$truth, sim_big$truth)
  )
  
  truth_small <- as.numeric(sim_small$truth %in% c("R1", "R2", "R3", "R4"))
  
  # --- Step D: Fit MIMOSA2 on Combined Data (sim_all) ---
  fit_all <- MIMOSA2(
        Ntot    = sim_all$Ntot,
        ns1     = sim_all$ns1,
        nu1     = sim_all$nu1,
        ns0     = sim_all$ns0,
        nu0     = sim_all$nu0,
        maxit   = 30,
        verbose = FALSE
      )
  
  # Extract predicted probabilities for the sim_small subset (rows 1:p_small)
  prob_mimosa_combined <- if (!is.null(fit_all)) {
    rowSums(fit_all$z[1:p_small, 1:4, drop = FALSE])
  } else {
    rep(NA_real_, p_small)
  }
  
  # --- Step E: Fit MIMOSA2 Independently on sim_small ---
  fit_small <- MIMOSA2(
        Ntot    = sim_small$Ntot,
        ns1     = sim_small$ns1,
        nu1     = sim_small$nu1,
        ns0     = sim_small$ns0,
        nu0     = sim_small$nu0,
        maxit   = 30,
        verbose = FALSE
      )
  prob_mimosa_independent <- if (!is.null(fit_small)) {
    rowSums(fit_small$z[, 1:4, drop = FALSE])
  } else {
    rep(NA_real_, p_small)
  }
  
  # --- Step F: Fit DiD GLM on sim_small ---
  prob_did_glm <- DiD_GLM(
    sim_small$Ntot, 
    sim_small$ns1, 
    sim_small$nu1, 
    sim_small$ns0, 
    sim_small$nu0
  )
  
  # --- Step G: Store Observation-Level Results (results_continuous format) ---
  obs_results_list[[i]] <- data.frame(
    Run_ID                   = i,
    Distribution             = dist,
    Res_prop                 = comp_nm,
    P_small                  = p_small,
    P_big                    = p_big,
    Small_effect             = eff_small,
    Large_effect             = eff_large,
    Cell_range               = rng_nm,
    Phi                      = phi,
    Replication              = rep_id,
    Subject_id               = 1:p_small,
    Truth                    = truth_small,
    MIMOSA2_prob_independent = prob_mimosa_independent,
    MIMOSA2_prob_combined    = prob_mimosa_combined,
    DiD_GLM_prob             = prob_did_glm,
    stringsAsFactors         = FALSE
  )
  
  # --- Step H: Store Metadata Summary per Run ---
  summary_results_list[[i]] <- data.frame(
    Run_ID          = i,
    Distribution    = dist,
    Res_prop        = comp_nm,
    P_small         = p_small,
    P_big           = p_big,
    Small_effect    = eff_small,
    Large_effect    = eff_large,
    Cell_range      = rng_nm,
    Phi             = phi,
    Replication     = rep_id,
    Status_Combined = if (!is.null(fit_all)) "Success" else "Crash/Timeout",
    Status_Indep    = if (!is.null(fit_small)) "Success" else "Crash/Timeout",
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 3. MERGE AND SAVE RESULTS
# -----------------------------------------------------------------------------
results_continuous <- do.call(rbind, obs_results_list)
results_summary    <- do.call(rbind, summary_results_list)

if (!dir.exists("_simulations")) dir.create("_simulations", recursive = TRUE)

save(results_summary, 
     results_continuous, 
     file = "_simulations/Simulation_Small_Vs_Combined_Continuous.Rdata")

message("Simulations completed successfully. Continuous results saved.")