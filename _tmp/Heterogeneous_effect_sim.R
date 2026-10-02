my_libs <- "/scratch/abrmoe030/R_libs"
.libPaths(c(my_libs, .libPaths()))

library(R.utils)
library(dplyr)
library(parallel)

n_workers   <- as.numeric(Sys.getenv("SLURM_NTASKS", 2))
timeout_sec <- 180   # hard per-task kill threshold
set.seed(4471)

# -----------------------------------------------------------------------------
# 1. PARAMETERS & SIMULATION MATRIX SETUP
# -----------------------------------------------------------------------------
responders20 = c(rep(0.20 / 4, 4), rep(0.80 / 4, 4))
responders50 = c(rep(0.50 / 4, 4), rep(0.50 / 4, 4))
responders80 = c(rep(0.80 / 4, 4), rep(0.20 / 4, 4))

component_list = list(
  "Prop_0.20" = responders20,
  "Prop_0.50" = responders50,
  "Prop_0.80" = responders80
)

rng_list = list(
  "High"          = c(250000, 250000),
  "Medium"        = c(100000, 100000),
  "Low"           = c(15000, 15000),
  "Very Low"      = c(7000, 7000),
  "Extremely Low" = c(3000, 3000)
)

stresstest_mat = expand.grid(
  Distribution = "Beta",
  Comp_Name    = names(component_list),
  P            = list(c(20, 80), c(50, 50), c(80, 20)),
  effect       = c(5e-2, 1e-2, 1e-3, 8e-4, 6.25e-4, 5e-4) |>
    rev() |>
    combn(2, simplify = FALSE) |>
    (\(x) Filter(\(z) z[1] < z[2] && z[1] %in% c(6.25e-4, 5e-4), x))(),
  Rng_Name     = names(rng_list),
  Phi          = 5000,
  Replication  = 1:50,
  stringsAsFactors = FALSE
)

log_dir   <- "/scratch/abrmoe030/projects/mimosa2/_tmp"
start_log <- file.path(log_dir, "task_start.log")
end_log   <- file.path(log_dir, "task_end.log")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
file.create(start_log)
file.create(end_log)

# -----------------------------------------------------------------------------
# 2. PER-TASK WORKER
# -----------------------------------------------------------------------------
run_simulation <- function(i) {
  dist      <- stresstest_mat$Distribution[i]
  comp_nm   <- stresstest_mat$Comp_Name[i]
  p_small   <- stresstest_mat$P[[i]][1]
  p_big     <- stresstest_mat$P[[i]][2]
  eff_large <- stresstest_mat$effect[[i]][1]
  eff_small <- stresstest_mat$effect[[i]][2]
  rng_nm    <- stresstest_mat$Rng_Name[i]
  phi       <- stresstest_mat$Phi[i]
  rep_id    <- stresstest_mat$Replication[i]
  
  active_components <- component_list[[comp_nm]]
  active_rng        <- rng_list[[rng_nm]]
  
  sim_small <- simulate_MIMOSA2_alt_prior(effect = eff_small, phi = phi, P = p_small,
                                          prior = dist, components = active_components,
                                          rng = active_rng)
  
  sim_big <- simulate_MIMOSA2_alt_prior(effect = eff_large, phi = phi, P = p_big,
                                        prior = dist, components = active_components,
                                        rng = active_rng)
  
  sim_all <- list(
    Ntot  = rbind(sim_small$Ntot, sim_big$Ntot),
    ns1   = c(sim_small$ns1, sim_big$ns1),
    nu1   = c(sim_small$nu1, sim_big$nu1),
    ns0   = c(sim_small$ns0, sim_big$ns0),
    nu0   = c(sim_small$nu0, sim_big$nu0),
    truth = c(sim_small$truth, sim_big$truth)
  )
  
  truth_small <- as.numeric(sim_small$truth %in% c("R1", "R2", "R3", "R4"))
  
  fit_all <- MIMOSA2(Ntot = sim_all$Ntot, ns1 = sim_all$ns1, nu1 = sim_all$nu1,
                     ns0 = sim_all$ns0, nu0 = sim_all$nu0, maxit = 30, verbose = FALSE)
  
  prob_mimosa_combined <- if (!is.null(fit_all)) {
    rowSums(fit_all$z[1:p_small, 1:4, drop = FALSE])
  } else rep(NA_real_, p_small)
  
  fit_small <- MIMOSA2(Ntot = sim_small$Ntot, ns1 = sim_small$ns1, nu1 = sim_small$nu1,
                       ns0 = sim_small$ns0, nu0 = sim_small$nu0, maxit = 30, verbose = FALSE)
  
  prob_mimosa_independent <- if (!is.null(fit_small)) {
    rowSums(fit_small$z[, 1:4, drop = FALSE])
  } else rep(NA_real_, p_small)
  
  prob_did_glm <- DiD_GLM(sim_small$Ntot, sim_small$ns1, sim_small$nu1,
                          sim_small$ns0, sim_small$nu0)
  
  list(continuous = data.frame(
    Run_ID = i, Distribution = dist, Res_prop = comp_nm,
    P_small = p_small, P_big = p_big,
    Small_effect = eff_small, Large_effect = eff_large,
    Cell_range = rng_nm, Phi = phi, Replication = rep_id,
    Subject_id = 1:p_small, Truth = truth_small,
    MIMOSA2_prob_independent = prob_mimosa_independent,
    MIMOSA2_prob_combined    = prob_mimosa_combined,
    DiD_GLM_prob = prob_did_glm,
    stringsAsFactors = FALSE
  ))
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
    
    rng_nm <- stresstest_mat$Rng_Name[tid]
    rep_id <- stresstest_mat$Replication[tid]
    eff_sm <- stresstest_mat$effect[[tid]][2]
    eff_lg <- stresstest_mat$effect[[tid]][1]
    
    job <- mcparallel(run_simulation(tid), silent = TRUE)
    key <- as.character(job$pid)
    running[[key]] <- list(job = job, task_id = tid, start = Sys.time())
    
    cat(sprintf("[%s] Launched task %d/%d (pid=%s) Rep=%d CellRange=%s SmallEff=%.2e LargeEff=%.2e\n",
                format(Sys.time(), "%H:%M:%S"), tid, length(task_ids), key,
                rep_id, rng_nm, eff_sm, eff_lg),
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
      mccollect(entry$job, wait = FALSE)
      cat(sprintf("[%s] KILLED task %d (pid=%s) after %.0fs\n",
                  format(Sys.time(), "%H:%M:%S"), entry$task_id, key, elapsed),
          file = end_log, append = TRUE)
      running[[key]] <- NULL
    }
  }
  
  n_done <- sum(!sapply(results, is.null))
  if (n_done %% 50 == 0 && n_done > 0) {
    done_idx <- which(!sapply(results, is.null))
    partial_continuous <- do.call(rbind, lapply(results[done_idx], function(x) x$continuous))
    if (!dir.exists("_simulations")) dir.create("_simulations", recursive = TRUE)
    save(partial_continuous, file = "_simulations/EffectSize_Heterogeneity_partial.Rdata")
  }
}

# -----------------------------------------------------------------------------
# 4. FINAL MERGE
# -----------------------------------------------------------------------------
master_list <- lapply(seq_along(results), function(tid) {
  r <- results[[tid]]
  if (!is.null(r) && !is.null(r$continuous)) return(r)
  
  comp_nm <- stresstest_mat$Comp_Name[tid]
  p_small <- stresstest_mat$P[[tid]][1]
  p_big   <- stresstest_mat$P[[tid]][2]
  eff_lg  <- stresstest_mat$effect[[tid]][1]
  eff_sm  <- stresstest_mat$effect[[tid]][2]
  rng_nm  <- stresstest_mat$Rng_Name[tid]
  phi     <- stresstest_mat$Phi[tid]
  rep_id  <- stresstest_mat$Replication[tid]
  
  list(continuous = data.frame(
    Run_ID = tid, Distribution = "Beta", Res_prop = comp_nm,
    P_small = p_small, P_big = p_big,
    Small_effect = eff_sm, Large_effect = eff_lg,
    Cell_range = rng_nm, Phi = phi, Replication = rep_id,
    Subject_id = 1:p_small, Truth = NA,
    MIMOSA2_prob_independent = NA, MIMOSA2_prob_combined = NA,
    DiD_GLM_prob = NA, stringsAsFactors = FALSE
  ))
})

results_continuous <- do.call(rbind, lapply(master_list, function(x) x$continuous))

if (!dir.exists("_simulations")) dir.create("_simulations", recursive = TRUE)
save(results_continuous, file = "_simulations/EffectSize_Heterogeneity.Rdata")

message("All simulations completed successfully.")