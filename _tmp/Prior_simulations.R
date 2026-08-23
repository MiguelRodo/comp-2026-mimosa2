# MIMOSA2 Simulation 
# Isabella and Tayyeb
# 26 June 2026

#-------------------------------------------------------------------------------
my_libs <- "/scratch/abrmoe030/R_libs"
.libPaths(c(my_libs, .libPaths()))

library(future)
library(future.apply)
library(R.utils)
plan(multicore, workers = as.numeric(Sys.getenv("SLURM_NTASKS", 2)))

cat("SLURM_NTASKS as seen by R:", Sys.getenv("SLURM_NTASKS", "NOT SET"), "\n")

responders50  = c(rep(0.50 / 4, 4), rep(0.50 / 4, 4))

# Combine into a named list for tracking:
component_list = list(
  "Prop_0.50" = responders50
)

# Cell count scenarios: 
# Combine your rng cell counts into a named list for tracking:
rng_list = list(
  "High"  = c(250000, 250000),
  "Medium" = c(100000, 100000),
  "Low"     = c(15000, 15000),
  "Very Low"= c(7000,7000)
)

# Build simulation grid:
stresstest_mat = expand.grid(
  Distribution_Phi   = list(c("Beta",10000),
                            c("EG",113),
                            c("LN",2.05),
                            c("SX",10000),
                            c("BB",10000)),
  Comp_Name      = names(component_list),
  P              = c(50),
  Effect         = c(1e-3, 2.5e-4, 1.25e-4, 6.25e-5),
  Rng_Name       = names(rng_list),
  Replication    = 1:30, 
  stringsAsFactors = FALSE
)

#-------------------------------------------------------------------------------
# Main Simulation Worker
#-------------------------------------------------------------------------------
log_dir  <- "/scratch/abrmoe030/projects/mimosa2/_tmp"
log_file <- file.path(log_dir, "sim_progress.log")

if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
file.create(log_file)

run_single_simulation <- function(i) {
  dist     = stresstest_mat$Distribution_Phi[[i]][1]
  comp_nm  = stresstest_mat$Comp_Name[i]
  p        = stresstest_mat$P[i]
  eff      = stresstest_mat$Effect[i]
  rng_nm   = stresstest_mat$Rng_Name[i]
  phi      = as.numeric(stresstest_mat$Distribution_Phi[[i]][2])
  rep_id   = stresstest_mat$Replication[i]
  
  active_components = component_list[[comp_nm]]
  active_rng        = rng_list[[rng_nm]]
  
  cat(paste0("Running simulation ", i, "/", nrow(stresstest_mat),
             ": Distribution = ", dist,
             ", Responder Prop = ", comp_nm,
             ", P = ", p,
             ", Effect = ", format(eff, scientific=FALSE, drop0trailing=TRUE),
             ", Rep = ", rep_id,
             ", Cell Range = ", rng_nm,
             ", phi = ", phi, "\n"),
              file = log_file, append = TRUE)

  
  # Simulate:
  sim = simulate_MIMOSA2_alt_prior(
    effect     = eff,
    phi        = phi,
    P          = p,
    prior      = dist,
    components = active_components,
    rng        = active_rng
  )
  
  # Define truth responder logic:
  true_responder = as.numeric(sim$truth %in% c("R1","R2","R3","R4"))
  
  # Call the DiD_GLM function sourced from Non_beta_simulations
  did_glm_prob = DiD_GLM(sim$Ntot, sim$ns1, sim$nu1, sim$ns0, sim$nu0)
  
  #Call the DiD function with shrinkage
  
  # Log-fold change (LFC) using matrix column names
  prop_s = sim$ns1 / sim$Ntot[, "ns1"]
  prop_u = sim$nu1 / sim$Ntot[, "nu1"]
  log_fold_change = log2((prop_s + 1e-5) / (prop_u + 1e-5))
  
  # Fit MIMOSA2 Model:
  fit_error = FALSE        
fit <- tryCatch({
  R.utils::withTimeout({
    MIMOSA2(
      Ntot    = sim$Ntot,
      ns1     = sim$ns1,
      nu1     = sim$nu1,
      ns0     = sim$ns0,
      nu0     = sim$nu0,
      maxit   = 30,
      verbose = FALSE
    )
  }, timeout = 420, onTimeout = "error")
},
error = function(e) {
  fit_error <<- TRUE
  return(NULL)
})
  
  # Initialize point metrics:
  status     = "Success"
  iterations = NA_real_
  TPR_001    = NA_real_
  tFDR_001   = NA_real_
  tFDR_005   = NA_real_
  
  # Evaluate MIMOSA metrics:
  if (fit_error == FALSE && is.null(fit) == FALSE) {
    iterations = length(fit$inds)
    mimosa_prob = rowSums(fit$z[, 1:4, drop = FALSE])
    
    rescall_001 = getResponse(fit, threshold=0.01) 
    rescall_005 = getResponse(fit, threshold=0.05) 
    
    if (any(true_responder == 1)) {
      TPR_001 = sum(rescall_001 & true_responder == 1) / sum(true_responder == 1)
    } else {
      TPR_001 = NA
    } 
    
    tFDR_001 = if (sum(rescall_001) > 0) {
      sum(rescall_001 & true_responder != 1) / sum(rescall_001)
    } else {
      0
    }
    
    tFDR_005 = if (sum(rescall_005) > 0) {
      sum(rescall_005 & true_responder != 1) / sum(rescall_005)
    } else {
      0
    }
    
  } else {
    status = "Optimisation crash"
    mimosa_prob = rep(NA_real_, p)
  }
  
  # Store summary results:
  row_res = data.frame(
    Distribution = dist,
    Res_prop     = comp_nm,
    P            = p,
    Effect       = eff,
    Cell_range   = rng_nm,
    Phi          = phi,
    Replication  = rep_id,
    Status       = status,
    Iterations   = iterations,
    TPR_001      = TPR_001,
    tFDR_001     = tFDR_001,
    tFDR_005     = tFDR_005,
    stringsAsFactors = FALSE
  )
  
  # Store observation-level results:
  scenario_obs = data.frame(
    Distribution = dist,
    Res_prop     = comp_nm,
    P            = p,
    Effect       = eff,
    Cell_range   = rng_nm,
    Phi          = phi,
    Replication  = rep_id,
    Subject_id   = 1:p,
    Truth        = true_responder,
    MIMOSA2_prob = mimosa_prob,
    DiD_GLM_prob = did_glm_prob,
    Log2_FC      = log_fold_change,
    stringsAsFactors = FALSE
  )
  
  return(list(summary = row_res, continuous = scenario_obs))
}

message("Starting parallel simulations...")

master_obs_list <- future_lapply(1:nrow(stresstest_mat), run_single_simulation, future.seed = TRUE)

results_summary    <- do.call(rbind, lapply(master_obs_list, function(x) x$summary))
results_continuous <- do.call(rbind, lapply(master_obs_list, function(x) x$continuous))

if (!dir.exists("_simulations")) dir.create("_simulations", recursive = TRUE)

# Save results: 
save(results_summary, 
     results_continuous,
     file = '_simulations/Simulation_3.0.Rdata')