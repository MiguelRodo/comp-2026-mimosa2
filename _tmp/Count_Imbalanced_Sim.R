
my_libs <- "/scratch/abrmoe030/R_libs"
.libPaths(c(my_libs, .libPaths()))

library(R.utils)
library(dplyr)
library(future)
library(future.apply)
plan(multicore, workers = as.numeric(Sys.getenv("SLURM_NTASKS", 2)))

# -----------------------------------------------------------------------------
# 1. PARAMETERS & SIMULATION MATRIX SETUP
# -----------------------------------------------------------------------------
responders25 = c(rep(0.25 / 4, 4), rep(0.75 / 4, 4))

component_list = list(
  "Prop_0.25" = responders25
)

rng_list = list(
  "1.00_Count"  = c(150000,15000),
  ".90_Count"   = c(rep(150000,2),rep(0.9*150000,2),rep(150000,4)),
  ".75_Count"   = c(rep(150000,2),rep(0.75*150000,2),rep(150000,4)),
  ".50_Count"   = c(rep(150000,2),rep(0.5*150000,2),rep(150000,4)),
  ".25_Count"   = c(rep(150000,2),rep(0.25*150000,2),rep(150000,4)),
  ".10_Count"   = c(rep(150000,2),rep(0.10*150000,2),rep(150000,4)),
  ".050_Count"   = c(rep(150000,2),rep(0.05*150000,2),rep(150000,4)),
  ".025_Count"   = c(rep(150000,2),rep(0.025*150000,2),rep(150000,4))
)

stresstest_mat = expand.grid(
  Distribution   = "Beta",
  Comp_Name      = names(component_list),
  P              = 100,
  Effect         = c(5e-4),
  Rng_Name       = names(rng_list),
  Phi            = 5000,
  Replication    = 1:50, 
  stringsAsFactors = FALSE
)

# Containers for observation-level and summary results
obs_results_list     <- vector("list", nrow(stresstest_mat))

# -----------------------------------------------------------------------------
# 2. MAIN FOR LOOP
# -----------------------------------------------------------------------------


log_dir  <- "/scratch/abrmoe030/projects/mimosa2/_tmp"
log_file <- file.path(log_dir, "sim_progress.log")

run_simulation <- function(i) {
  
  # Extract parameters for row i
  dist       <- stresstest_mat$Distribution[i]
  comp_nm    <- stresstest_mat$Comp_Name[i]
  P          <- stresstest_mat$P[i]
  Effect     <- stresstest_mat$Effect[i]
  rng_nm     <- stresstest_mat$Rng_Name[i]
  phi        <- stresstest_mat$Phi[i]
  rep_id     <- stresstest_mat$Replication[i]
  
  active_components <- component_list[[comp_nm]]
  active_rng        <- rng_list[[rng_nm]]
  
  cat(sprintf("Run %d/%d: Rep=%d, CellRange=%s \n",
              i, nrow(stresstest_mat), rep_id, rng_nm),
      file = log_file, append = TRUE)
  
  # --- Step A: Generate Small and Big Simulations ---
  sim <- simulate_MIMOSA2_alt_prior(
    effect     = Effect,
    phi        = phi,
    P          = P,
    prior      = dist,
    components = active_components,
    rng        = active_rng
  )
  
  truth <- as.numeric(sim$truth %in% c("R1", "R2", "R3", "R4"))
  
  # --- Step D: Fit MIMOSA2 on Combined Data (sim_all) ---
  fit_all <- MIMOSA2(
    Ntot    = sim$Ntot,
    ns1     = sim$ns1,
    nu1     = sim$nu1,
    ns0     = sim$ns0,
    nu0     = sim$nu0,
    maxit   = 30,
    verbose = FALSE
  )
  
  # Extract predicted probabilities for the sim_small subset (rows 1:p_small)
  prob_mimosa<- if (!is.null(fit_all)) {
    rowSums(fit_all$z[, 1:4, drop = FALSE])
  } else {
    rep(NA_real_, P)
  }
  
  
  # --- Step F: Fit DiD GLM on sim_small ---
  prob_did_glm <- DiD_GLM(
    sim$Ntot, 
    sim$ns1, 
    sim$nu1, 
    sim$ns0, 
    sim$nu0
  )
  
  # --- Step G: Store Observation-Level Results (results_continuous format) ---
  obs_results_list <- data.frame(
    Run_ID                   = i,
    Distribution             = dist,
    Res_prop                 = comp_nm,
    P                        = P,
    Effect                   = Effect,
    Cell_range               = rng_nm,
    Phi                      = phi,
    Replication              = rep_id,
    Subject_id               = 1:P,
    Truth                    = truth,
    MIMOSA2_prob             = prob_mimosa,
    DiD_GLM_prob             = prob_did_glm,
    stringsAsFactors         = FALSE
  )
  
  return(list(continuous=obs_results_list))
}

# -----------------------------------------------------------------------------
# 3. MERGE AND SAVE RESULTS
# -----------------------------------------------------------------------------

message("Starting parallel simulations...")

master_obs_list <- future_lapply(1:nrow(stresstest_mat), run_simulation, future.seed = TRUE,
                                 future.scheduling = Inf)

results_continuous <- do.call(rbind, lapply(master_obs_list, function(x) x$continuous))


if (!dir.exists("_simulations")) dir.create("_simulations", recursive = TRUE)

save(results_continuous, 
     file = "_simulations/Count_Imbalance.Rdata")

message("Simulations completed successfully. Continuous results saved.")