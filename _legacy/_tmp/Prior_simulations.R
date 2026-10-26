my_libs <- "/scratch/abrmoe030/R_libs"
.libPaths(c(my_libs, .libPaths()))

library(R.utils)
library(dplyr)
library(parallel)

n_workers   <- as.numeric(Sys.getenv("SLURM_NTASKS", 2))
timeout_sec <- 180   # hard per-task kill threshold
set.seed(9854)

# -----------------------------------------------------------------------------
# 1. PARAMETERS & SIMULATION MATRIX SETUP
# -----------------------------------------------------------------------------
responders50 = c(rep(0.50 / 4, 4), rep(0.50 / 4, 4))
component_list = list("Prop_0.50" = responders50)

rng_list = list(
  "High"          = c(250000, 250000),
  "Medium"        = c(100000, 100000),
  "Low"           = c(15000, 15000),
  "Very Low"      = c(7000, 7000),
  "Extremely Low" = c(3000, 3000)
)

stresstest_mat = expand.grid(
  Distribution_Phi = list(c("Beta", 10000),
                          c("EG", 113),
                          c("LN", 2.05),
                          c("SX", 10000),
                          c("BB", 10000)),
  Comp_Name   = names(component_list),
  P           = c(20, 50, 100),
  Effect      = c(1e-3, 2.5e-4, 1.25e-4, 6.25e-5),
  Rng_Name    = names(rng_list),
  Replication = 1:200,
  stringsAsFactors = FALSE
)

log_dir   <- "/scratch/abrmoe030/projects/mimosa2/_tmp"
start_log <- file.path(log_dir, "task_start.log")
end_log   <- file.path(log_dir, "task_end.log")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
file.create(start_log)
file.create(end_log)

# -----------------------------------------------------------------------------
# 2. PER-TASK WORKER (returns BOTH summary and continuous data frames)
# -----------------------------------------------------------------------------
run_simulation <- function(i) {
  dist    <- stresstest_mat$Distribution_Phi[[i]][1]
  comp_nm <- stresstest_mat$Comp_Name[i]
  p       <- stresstest_mat$P[i]
  eff     <- stresstest_mat$Effect[i]
  rng_nm  <- stresstest_mat$Rng_Name[i]
  phi     <- as.numeric(stresstest_mat$Distribution_Phi[[i]][2])
  rep_id  <- stresstest_mat$Replication[i]
  
  active_components <- component_list[[comp_nm]]
  active_rng        <- rng_list[[rng_nm]]
  
  sim <- simulate_MIMOSA2_alt_prior(effect = eff, phi = phi, P = p,
                                    prior = dist, components = active_components,
                                    rng = active_rng)
  
  true_responder <- as.numeric(sim$truth %in% c("R1", "R2", "R3", "R4"))
  
  # Direct call - no internal timeout needed, the outer scheduler hard-kills
  # the whole task (simulate + fit) if it runs past timeout_sec.
  fit <- MIMOSA2(Ntot = sim$Ntot, ns1 = sim$ns1, nu1 = sim$nu1,
                 ns0 = sim$ns0, nu0 = sim$nu0, maxit = 30, verbose = FALSE)
  
  did_glm_prob <- DiD_GLM(sim$Ntot, sim$ns1, sim$nu1, sim$ns0, sim$nu0)
  
  prop_s <- sim$ns1 / sim$Ntot[, "ns1"]
  prop_u <- sim$nu1 / sim$Ntot[, "nu1"]
  log_fold_change <- log2((prop_s + 1e-5) / (prop_u + 1e-5))
  
  status     <- "Success"
  iterations <- NA_real_
  TPR_001    <- NA_real_
  tFDR_001   <- NA_real_
  tFDR_005   <- NA_real_
  
  if (!is.null(fit)) {
    iterations  <- length(fit$inds)
    mimosa_prob <- rowSums(fit$z[, 1:4, drop = FALSE])
    
    rescall_001 <- getResponse(fit, threshold = 0.01)
    rescall_005 <- getResponse(fit, threshold = 0.05)
    
    TPR_001 <- if (any(true_responder == 1)) {
      sum(rescall_001 & true_responder == 1) / sum(true_responder == 1)
    } else NA
    
    tFDR_001 <- if (sum(rescall_001) > 0) sum(rescall_001 & true_responder != 1) / sum(rescall_001) else 0
    tFDR_005 <- if (sum(rescall_005) > 0) sum(rescall_005 & true_responder != 1) / sum(rescall_005) else 0
  } else {
    status      <- "Optimisation crash"
    mimosa_prob <- rep(NA_real_, p)
  }
  
  row_res <- data.frame(
    Distribution = dist, Res_prop = comp_nm, P = p, Effect = eff,
    Cell_range = rng_nm, Phi = phi, Replication = rep_id,
    Status = status, Iterations = iterations,
    TPR_001 = TPR_001, tFDR_001 = tFDR_001, tFDR_005 = tFDR_005,
    stringsAsFactors = FALSE
  )
  
  scenario_obs <- data.frame(
    Distribution = dist, Res_prop = comp_nm, P = p, Effect = eff,
    Cell_range = rng_nm, Phi = phi, Replication = rep_id,
    Subject_id = 1:p, Truth = true_responder,
    MIMOSA2_prob = mimosa_prob, DiD_GLM_prob = did_glm_prob,
    Log2_FC = log_fold_change, stringsAsFactors = FALSE
  )
  
  list(summary = row_res, continuous = scenario_obs)
}

# -----------------------------------------------------------------------------
# 3. HARD-TIMEOUT PARALLEL SCHEDULER (same mcparallel/mccollect pattern)
# -----------------------------------------------------------------------------
task_ids <- 1:nrow(stresstest_mat)
pending  <- task_ids
running  <- list()
results  <- vector("list", length(task_ids))

message("Starting parallel simulations (base parallel, no future)...")

while (length(pending) > 0 || length(running) > 0) {
  
  while (length(running) < n_workers && length(pending) > 0) {
    tid <- pending[1]
    pending <- pending[-1]
    
    dist   <- stresstest_mat$Distribution_Phi[[tid]][1]
    rng_nm <- stresstest_mat$Rng_Name[tid]
    rep_id <- stresstest_mat$Replication[tid]
    
    job <- mcparallel(run_simulation(tid), silent = TRUE)
    key <- as.character(job$pid)
    running[[key]] <- list(job = job, task_id = tid, start = Sys.time())
    
    cat(sprintf("[%s] Launched task %d/%d (pid=%s) Dist=%s Rep=%d CellRange=%s\n",
                format(Sys.time(), "%H:%M:%S"), tid, length(task_ids), key,
                dist, rep_id, rng_nm),
        file = start_log, append = TRUE)
  }
  
  Sys.sleep(1)
  
  for (key in names(running)) {
    entry <- running[[key]]
    res <- mccollect(entry$job, wait = FALSE)
    
    if (!is.null(res)) {
      results[[entry$task_id]] <- res[[key]]
      elapsed <- as.numeric(difftime(Sys.time(), entry$start, units = "secs"))
      cat(sprintf("[%s] Completed task %d (pid=%s) in %.1fs\n",
                  format(Sys.time(), "%H:%M:%S"), entry$task_id, key, elapsed),
          file = end_log, append = TRUE)
      running[[key]] <- NULL
      next
    }
    
    elapsed <- as.numeric(difftime(Sys.time(), entry$start, units = "secs"))
    if (elapsed > timeout_sec) {
      tools::pskill(as.integer(key), tools::SIGKILL)
      mccollect(entry$job, wait = FALSE)  # reap
      cat(sprintf("[%s] KILLED task %d (pid=%s) after %.0fs\n",
                  format(Sys.time(), "%H:%M:%S"), entry$task_id, key, elapsed),
          file = end_log, append = TRUE)
      running[[key]] <- NULL
    }
  }
  
  # Periodic checkpoint save every ~50 completed tasks
  n_done <- sum(!sapply(results, is.null))
  if (n_done %% 50 == 0 && n_done > 0) {
    done_idx <- which(!sapply(results, is.null))
    partial_summary    <- do.call(rbind, lapply(results[done_idx], function(x) x$summary))
    partial_continuous <- do.call(rbind, lapply(results[done_idx], function(x) x$continuous))
    if (!dir.exists("_simulations")) dir.create("_simulations", recursive = TRUE)
    save(partial_summary, partial_continuous,
         file = "_simulations/Simulation_3.0_partial.Rdata")
  }
}

# -----------------------------------------------------------------------------
# 4. FINAL MERGE (NA placeholder rows for any killed/failed tasks)
# -----------------------------------------------------------------------------
master_list <- lapply(seq_along(results), function(tid) {
  r <- results[[tid]]
  if (!is.null(r) && !is.null(r$summary) && !is.null(r$continuous)) return(r)
  
  dist   <- stresstest_mat$Distribution_Phi[[tid]][1]
  comp_nm<- stresstest_mat$Comp_Name[tid]
  p      <- stresstest_mat$P[tid]
  eff    <- stresstest_mat$Effect[tid]
  rng_nm <- stresstest_mat$Rng_Name[tid]
  phi    <- as.numeric(stresstest_mat$Distribution_Phi[[tid]][2])
  rep_id <- stresstest_mat$Replication[tid]
  
  list(
    summary = data.frame(
      Distribution = dist, Res_prop = comp_nm, P = p, Effect = eff,
      Cell_range = rng_nm, Phi = phi, Replication = rep_id,
      Status = "Crash/Timeout", Iterations = NA, TPR_001 = NA,
      tFDR_001 = NA, tFDR_005 = NA, stringsAsFactors = FALSE
    ),
    continuous = data.frame(
      Distribution = dist, Res_prop = comp_nm, P = p, Effect = eff,
      Cell_range = rng_nm, Phi = phi, Replication = rep_id,
      Subject_id = 1:p, Truth = NA, MIMOSA2_prob = NA, DiD_GLM_prob = NA,
      Log2_FC = NA, stringsAsFactors = FALSE
    )
  )
})

results_summary    <- do.call(rbind, lapply(master_list, function(x) x$summary))
results_continuous <- do.call(rbind, lapply(master_list, function(x) x$continuous))

if (!dir.exists("_simulations")) dir.create("_simulations", recursive = TRUE)
save(results_summary, results_continuous, file = "_simulations/Simulation_3.0.Rdata")

message("All simulations completed successfully.")