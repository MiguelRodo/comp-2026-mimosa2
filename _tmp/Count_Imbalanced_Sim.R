
my_libs <- "/scratch/abrmoe030/R_libs"
.libPaths(c(my_libs, .libPaths()))

library(R.utils)
library(dplyr)
library(parallel)

n_workers   <- as.numeric(Sys.getenv("SLURM_NTASKS", 2))
timeout_sec <- 180   # hard per-task kill threshold

# -----------------------------------------------------------------------------
# 1. PARAMETERS & SIMULATION MATRIX SETUP  (unchanged from your script)
# -----------------------------------------------------------------------------
responders25 = c(rep(0.25 / 4, 4), rep(0.75 / 4, 4))
component_list = list("Prop_0.25" = responders25)

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
  Replication    = 1:50,
  stringsAsFactors = FALSE
)

log_dir    <- "/scratch/abrmoe030/projects/mimosa2/_tmp"
start_log  <- file.path(log_dir, "task_start.log")
end_log    <- file.path(log_dir, "task_end.log")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
file.create(start_log)
file.create(end_log)

# -----------------------------------------------------------------------------
# 2. PER-TASK WORKER  (unchanged from your script)
# -----------------------------------------------------------------------------
run_simulation <- function(i) {
  dist    <- stresstest_mat$Distribution[i]
  comp_nm <- stresstest_mat$Comp_Name[i]
  P       <- stresstest_mat$P[i]
  Effect  <- stresstest_mat$Effect[i]
  rng_nm  <- stresstest_mat$Rng_Name[i]
  phi     <- stresstest_mat$Phi[i]
  rep_id  <- stresstest_mat$Replication[i]

  active_components <- component_list[[comp_nm]]
  active_rng        <- rng_list[[rng_nm]]

  sim <- simulate_MIMOSA2_alt_prior(effect = Effect, phi = phi, P = P,
                                     prior = dist, components = active_components,
                                     rng = active_rng)

  truth <- as.numeric(sim$truth %in% c("R1", "R2", "R3", "R4"))

  fit_all <- MIMOSA2(Ntot = sim$Ntot, ns1 = sim$ns1, nu1 = sim$nu1,
                      ns0 = sim$ns0, nu0 = sim$nu0, maxit = 30, verbose = FALSE)

  prob_mimosa <- if (!is.null(fit_all)) rowSums(fit_all$z[, 1:4, drop = FALSE]) else rep(NA_real_, P)
  prob_did_glm <- DiD_GLM(sim$Ntot, sim$ns1, sim$nu1, sim$ns0, sim$nu0)

  list(continuous = data.frame(
    Run_ID = i, Distribution = dist, Res_prop = comp_nm, P = P, Effect = Effect,
    Cell_range = rng_nm, Phi = phi, Replication = rep_id, Subject_id = 1:P,
    Truth = truth, MIMOSA2_prob = prob_mimosa, DiD_GLM_prob = prob_did_glm,
    stringsAsFactors = FALSE
  ))
}

# -----------------------------------------------------------------------------
# 3. HARD-TIMEOUT PARALLEL SCHEDULER (parallel::mcparallel / mccollect)
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
    job <- mcparallel(run_simulation(tid), silent = TRUE)
    key <- as.character(job$pid)
    running[[key]] <- list(job = job, task_id = tid, start = Sys.time())
    cat(sprintf("[%s] Launched task %d/%d (pid=%s)\n",
                format(Sys.time(), "%H:%M:%S"), tid, length(task_ids), key),
        file = start_log, append = TRUE)
  }

  Sys.sleep(1)

  for (key in names(running)) {
    entry <- running[[key]]
    res <- mccollect(entry$job, wait = FALSE)

    if (!is.null(res)) {
      results[[entry$task_id]] <- res[[key]]
      cat(sprintf("[%s] Completed task %d (pid=%s)\n",
                  format(Sys.time(), "%H:%M:%S"), entry$task_id, key),
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
    partial <- do.call(rbind, Filter(Negate(is.null), lapply(results, function(x) x$continuous)))
    if (!dir.exists("_simulations")) dir.create("_simulations", recursive = TRUE)
    save(partial, file = "_simulations/Count_Imbalance_partial.Rdata")
  }
}

# -----------------------------------------------------------------------------
# 4. FINAL MERGE (fills NA placeholder rows for any killed/failed tasks)
# -----------------------------------------------------------------------------
master_obs_list <- lapply(seq_along(results), function(tid) {
  r <- results[[tid]]
  if (is.null(r) || is.null(r$continuous)) {
    Pi <- stresstest_mat$P[tid]
    data.frame(
      Run_ID = tid, Distribution = stresstest_mat$Distribution[tid],
      Res_prop = stresstest_mat$Comp_Name[tid], P = Pi,
      Effect = stresstest_mat$Effect[tid], Cell_range = stresstest_mat$Rng_Name[tid],
      Phi = stresstest_mat$Phi[tid], Replication = stresstest_mat$Replication[tid],
      Subject_id = 1:Pi, Truth = NA, MIMOSA2_prob = NA, DiD_GLM_prob = NA,
      stringsAsFactors = FALSE
    )
  } else r$continuous
})

results_continuous <- do.call(rbind, master_obs_list)

if (!dir.exists("_simulations")) dir.create("_simulations", recursive = TRUE)
save(results_continuous, file = "_simulations/Count_Imbalance.Rdata")
message("Simulations completed successfully.")