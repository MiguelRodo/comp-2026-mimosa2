
my_libs <- "/scratch/abrmoe030/R_libs"
.libPaths(c(my_libs, .libPaths()))

# CRITICAL: Prevent OpenMP/BLAS thread deadlocks during process forking
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

library(R.utils)
library(dplyr)
library(parallel)

# -----------------------------------------------------------------------------
# 1. PARAMETERS & SIMULATION MATRIX SETUP
# -----------------------------------------------------------------------------
responders25 = c(rep(0.25 / 4, 4), rep(0.75 / 4, 4))

component_list = list(
  "Prop_0.25" = responders25
)

rng_list = list(
  "1.00_Count_Ref"  = c(150000, 150000),
  ".90_Count_Imb"   = c(rep(150000, 2), rep(0.90 * 150000, 2), rep(150000, 4)),
  ".90_Count_Bal"   = c(rep((3 + 0.90) / 4 * 150000, 2)), 
  ".75_Count_Imb"   = c(rep(150000, 2), rep(0.75 * 150000, 2), rep(150000, 4)),
  ".75_Count_Bal"   = c(rep((3 + 0.75) / 4 * 150000, 2)), 
  ".50_Count_Imb"   = c(rep(150000, 2), rep(0.50 * 150000, 2), rep(150000, 4)),
  ".50_Count_Bal"   = c(rep((3 + 0.50) / 4 * 150000, 2)), 
  ".25_Count_Imb"   = c(rep(150000, 2), rep(0.25 * 150000, 2), rep(150000, 4)),
  ".25_Count_Bal"   = c(rep((3 + 0.25) / 4 * 150000, 2)), 
  ".10_Count_Imb"   = c(rep(150000, 2), rep(0.10 * 150000, 2), rep(150000, 4)),
  ".10_Count_Bal"   = c(rep((3 + 0.10) / 4 * 150000, 2)), 
  ".050_Count_Imb"  = c(rep(150000, 2), rep(0.05 * 150000, 2), rep(150000, 4)),
  ".050_Count_Bal"  = c(rep((3 + 0.05) / 4 * 150000, 2)), 
  ".025_Count_Imb"  = c(rep(150000, 2), rep(0.025 * 150000, 2), rep(150000, 4)),
  ".025_Count_Bal"  = c(rep((3 + 0.025) / 4 * 150000, 2))  
)

stresstest_mat = expand.grid(
  Distribution   = "Beta",
  Comp_Name      = names(component_list),
  P              = 100,
  Effect         = c(5e-4),
  Rng_Name       = names(rng_list),
  Phi            = 5000,
  Replication    = 1:30, 
  stringsAsFactors = FALSE
)

# Set directories
log_dir  <- "/scratch/abrmoe030/projects/mimosa2/_tmp"
out_tmp_dir <- file.path(log_dir, "rds_chunks")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
if (!dir.exists(out_tmp_dir)) dir.create(out_tmp_dir, recursive = TRUE)

log_file <- file.path(log_dir, "sim_progress.log")

# -----------------------------------------------------------------------------
# 2. WORKER FUNCTION WITH DISK-BACKED SAVING
# -----------------------------------------------------------------------------
run_simulation <- function(i) {
  
  # Define per-task output RDS file
  chunk_file <- file.path(out_tmp_dir, sprintf("res_run_%04d.rds", i))
  
  # Skip if already computed (enables easy job resumption if interrupted)
  if (file.exists(chunk_file)) {
    return(TRUE)
  }

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
  
  cat(sprintf("[%s] Run %d/%d: Rep=%d, CellRange=%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), i, nrow(stresstest_mat), rep_id, rng_nm),
      file = log_file, append = TRUE)
  
  # --- Step A: Generate Simulation ---
  sim <- simulate_MIMOSA2_alt_prior(
    effect     = Effect,
    phi        = phi,
    P          = P,
    prior      = dist,
    components = active_components,
    rng        = active_rng
  )
  
  truth <- as.numeric(sim$truth %in% c("R1", "R2", "R3", "R4"))
  
  # --- Step B: Fit MIMOSA2 ---
  fit_all <- tryCatch({
    MIMOSA2(
      Ntot    = sim$Ntot,
      ns1     = sim$ns1,
      nu1     = sim$nu1,
      ns0     = sim$ns0,
      nu0     = sim$nu0,
      maxit   = 30,
      verbose = FALSE
    )
  }, error = function(e) NULL)
  
  prob_mimosa <- if (!is.null(fit_all)) {
    rowSums(fit_all$z[, 1:4, drop = FALSE])
  } else {
    rep(NA_real_, P)
  }
  
  # --- Step C: Fit DiD GLM ---
  prob_did_glm <- DiD_GLM(
    sim$Ntot, 
    sim$ns1, 
    sim$nu1, 
    sim$ns0, 
    sim$nu0
  )
  
  # --- Step D: Store Observation-Level Results ---
  obs_df <- data.frame(
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
  
  # Save chunk immediately to disk to free worker RAM
  saveRDS(obs_df, file = chunk_file)
  
  return(TRUE)
}

# -----------------------------------------------------------------------------
# 3. PARALLEL EXECUTION VIA MCLAPPLY
# -----------------------------------------------------------------------------
n_cores <- as.numeric(Sys.getenv("SLURM_NTASKS", 2))
message(sprintf("Starting parallel execution using %d cores...", n_cores))

# Set parallel seed for reproducible simulations
set.seed(123)

# Execute via base R mclapply
res_status <- mclapply(
  X        = 1:nrow(stresstest_mat),
  FUN      = run_simulation,
  mc.cores = n_cores,
  mc.preschedule = TRUE
)

# -----------------------------------------------------------------------------
# 4. AGGREGATE RESULTS FROM DISK
# -----------------------------------------------------------------------------
message("Parallel tasks completed. Aggregating output files from disk...")

all_rds_files <- list.files(out_tmp_dir, pattern = "^res_run_.*\\.rds$", full.names = TRUE)

if (length(all_rds_files) == 0) {
  stop("No results found! Check logs for errors.")
}

# Read and combine all RDS files
results_continuous <- do.call(rbind, lapply(all_rds_files, readRDS))

# Save final dataset
if (!dir.exists("_simulations")) dir.create("_simulations", recursive = TRUE)
save(results_continuous, file = "_simulations/Count_Imbalance.Rdata")

# Clean up temporary RDS chunks
unlink(out_tmp_dir, recursive = TRUE)

message(sprintf("Successfully aggregated %d runs into _simulations/Count_Imbalance.Rdata", length(all_rds_files)))