# =============================================================================
# runner.R : reproducible, resumable parallel execution of a simulation study
# =============================================================================
#
# [CHANGE] Replaces the scheduler loop copied into each version-2 script.
# What it does differently (and why):
#
# 1. RANDOM NUMBERS (Morris 4.1, 4.1.1; audit B1). The master seed is set ONCE
#    per study with the L'Ecuyer-CMRG generator, and every dataset (task) gets
#    its own random-number stream (parallel::nextRNGStream). The stream seed
#    of every task is stored in the design table (the "states" dataset of
#    Morris Table 5). Any single dataset, e.g. one where MIMOSA2 failed, can
#    be regenerated exactly with regenerate_dataset(). Results no longer
#    depend on how many workers were used or in which order tasks finished.
#
# 2. ONE FILE PER TASK, RESUMABLE. Each task writes its own small .rds file.
#    If the job hits the wall-time or crashes, re-submitting the same command
#    skips completed tasks. Version 2 kept everything in memory and its
#    "checkpoint every 50 tasks" condition (n_done %% 50 == 0) re-saved the
#    whole partial result on every pass of the loop while n_done stayed at a
#    multiple of 50.
#
# 3. REPETITION-MAJOR ORDER. Tasks run in the order rep 1 of every scenario,
#    then rep 2 of every scenario, and so on. If the run stops early, every
#    scenario has (almost) the same number of repetitions, rather than some
#    scenarios being complete and others empty.
#
# 4. NESTED PROFILES. The design is built for the largest nsim of any profile
#    and the full (extended) grid; the standard profile is a subset of rows.
#    Task IDs and streams are therefore identical across profiles.
#
# 5. FAILURES ARE RECORDED, NOT DROPPED (Morris 5.1). A MIMOSA2 fit that errors
#    or exceeds FIT_TIMEOUT is stored with status "error"/"timeout". A task
#    whose whole process dies is recorded as "task_failed" when results are
#    combined, so every planned dataset appears in the results.
# =============================================================================

suppressPackageStartupMessages(library(parallel))

# ---- Design table with one RNG stream per task -----------------------------
# scenarios: data.frame, one row per scenario (data-generating mechanism),
#            must contain Scenario_ID and a logical column In_standard.
# Returns one row per (scenario, repetition) with Task_ID and Seed.
build_design <- function(study, scenarios, nsim_max = NSIM_MAX[[study]]) {
  stopifnot("Scenario_ID" %in% names(scenarios), "In_standard" %in% names(scenarios))
  stopifnot(!anyDuplicated(scenarios$Scenario_ID))
  grid <- expand.grid(s = seq_len(nrow(scenarios)), Rep = seq_len(nsim_max))
  grid <- grid[order(grid$Rep, grid$s), ]                 # repetition-major order
  design <- cbind(scenarios[grid$s, , drop = FALSE], Rep = grid$Rep)
  rownames(design) <- NULL
  design$Task_ID <- seq_len(nrow(design))
  # set the seed ONCE, then one independent stream per task
  old_kind <- RNGkind()[1]
  RNGkind("L'Ecuyer-CMRG")
  set.seed(MASTER_SEED[[study]])
  s <- .Random.seed
  seeds <- character(nrow(design))
  for (i in seq_len(nrow(design))) {
    seeds[i] <- paste(s, collapse = ",")
    s <- nextRNGStream(s)
  }
  design$Seed <- seeds
  attr(design, "end_state") <- paste(s, collapse = ",")   # Morris: store final state too
  attr(design, "study") <- study
  design
}

# Rows to run in the current profile
select_profile <- function(design, study, profile = PROFILE) {
  nsim <- NSIM[[profile]][[study]]
  keep <- design$Rep <= nsim
  if (profile %in% c("standard", "smoke")) keep <- keep & design$In_standard
  if (profile == "smoke") keep <- keep & design$In_smoke
  design[keep, , drop = FALSE]
}

set_task_seed <- function(seed_string) {
  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", as.integer(strsplit(seed_string, ",")[[1]]), envir = .GlobalEnv)
}

task_dir  <- function(study) file.path(OUT_DIR, study, "tasks")
task_file <- function(study, task_id) file.path(task_dir(study), sprintf("task_%06d.rds", task_id))

# ---- Run one task in the current process ------------------------------------
run_task <- function(study, row, task_fun) {
  set_task_seed(row$Seed)
  t0 <- Sys.time()
  res <- task_fun(row)
  res$task_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res$Task_ID <- row$Task_ID
  f <- task_file(study, row$Task_ID)
  tmp <- paste0(f, ".tmp", Sys.getpid())
  saveRDS(res, tmp)
  file.rename(tmp, f)        # atomic: a half-written file is never seen as done
  invisible(TRUE)
}

# ---- Parallel scheduler with an outer safety time limit ---------------------
# Each MIMOSA2 fit has its own FIT_TIMEOUT (methods.R). The outer limit only
# catches a task whose process hangs for some other reason.
run_study <- function(study, design, task_fun, n_fits = 1,
                      n_workers = N_WORKERS, chunk_id = CHUNK_ID, n_chunks = N_CHUNKS) {
  dir.create(task_dir(study), recursive = TRUE, showWarnings = FALSE)
  stopifnot(chunk_id >= 0, chunk_id < n_chunks)
  # [CHANGE v3.1, 6 Oct] SAFETY GUARD. On 6 Oct the baseline was re-submitted
  # after its task files had gone, so it silently restarted all 46,800 tasks
  # (and, before doing so, deleted the 'finished' marker and overwrote
  # run_info/design_full). If a combined results file already exists but no
  # task file of this chunk is found, stop instead: the study is already done.
  # Set MIMOSA2_FORCE_RERUN=1 to re-run it anyway.
  res_file <- file.path(OUT_DIR, study, sprintf("%s_results_%s.rds", study, PROFILE))
  todo0 <- select_profile(design, study)
  todo0 <- todo0[(seq_len(nrow(todo0)) - 1) %% n_chunks == chunk_id, , drop = FALSE]
  if (file.exists(res_file) && !any(file.exists(task_file(study, todo0$Task_ID))) &&
      !nzchar(Sys.getenv("MIMOSA2_FORCE_RERUN")))
    stop(sprintf(paste0("[%s] %s already exists but the tasks/ folder has none of this chunk's task files. ",
                        "The study has already been run and combined; nothing was changed. ",
                        "To re-run from scratch anyway: export MIMOSA2_FORCE_RERUN=1"), study, res_file))
  # this chunk is (re)starting: remove its old 'finished' marker
  unlink(file.path(OUT_DIR, study, sprintf("finished_%s_chunk%02d.txt", PROFILE, chunk_id)))
  dfile <- file.path(OUT_DIR, study, "design_full.rds")             # states dataset (Morris Table 5)
  tmpd <- paste0(dfile, ".tmp", Sys.getpid())                       # atomic write: array jobs may
  saveRDS(design, tmpd); file.rename(tmpd, dfile)                   # start at the same moment
  todo <- select_profile(design, study)
  todo <- todo[(seq_len(nrow(todo)) - 1) %% n_chunks == chunk_id, , drop = FALSE]
  done <- file.exists(task_file(study, todo$Task_ID))
  log_file <- file.path(OUT_DIR, study, sprintf("progress_chunk%02d.log", chunk_id))
  logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                            sprintf(...), "\n", sep = "", file = log_file, append = TRUE)
  logf("study=%s profile=%s chunk=%d/%d tasks=%d already_done=%d workers=%d",
       study, PROFILE, chunk_id + 1, n_chunks, nrow(todo), sum(done), n_workers)
  message(sprintf("[%s] %d tasks in this chunk, %d already done, %d to run on %d workers",
                  study, nrow(todo), sum(done), sum(!done), n_workers))
  saveRDS(list(profile = PROFILE, maxit = MAXIT, fit_timeout = FIT_TIMEOUT, alphas = ALPHAS,
               did_cc = if (exists("DID_CC")) DID_CC else 0,   # [v3.1] DiD variance correction used in this run
               master_seed = MASTER_SEED[[study]], nsim = NSIM[[PROFILE]][[study]],
               n_workers = n_workers, chunk = c(chunk_id, n_chunks),
               session = utils::sessionInfo(), started = Sys.time(),
               mimosa2_version = tryCatch(as.character(utils::packageVersion("MIMOSA2")), error = function(e) NA)),
          file.path(OUT_DIR, study, sprintf("run_info_chunk%02d.rds", chunk_id)))
  todo <- todo[!done, , drop = FALSE]
  if (nrow(todo) == 0) { message("Nothing to do."); mark_finished(study, chunk_id, n_chunks); return(invisible(TRUE)) }

  outer_timeout <- n_fits * FIT_TIMEOUT + 600
  # Windows (local testing): run sequentially
  if (.Platform$OS.type == "windows" || n_workers <= 1) {
    for (i in seq_len(nrow(todo))) {
      ok <- tryCatch({ run_task(study, todo[i, ], task_fun); TRUE },
                     error = function(e) { logf("ERROR task %d: %s", todo$Task_ID[i], conditionMessage(e)); FALSE })
      if (ok) logf("done task %d", todo$Task_ID[i])
    }
    mark_finished(study, chunk_id, n_chunks)
    return(invisible(TRUE))
  }

  pending <- seq_len(nrow(todo)); running <- list(); n_ok <- 0; n_bad <- 0
  t_start <- Sys.time(); last_report <- 0
  while (length(pending) > 0 || length(running) > 0) {
    while (length(running) < n_workers && length(pending) > 0) {
      i <- pending[1]; pending <- pending[-1]
      row <- todo[i, ]
      job <- mcparallel(run_task(study, row, task_fun), silent = TRUE, mc.set.seed = FALSE)
      running[[as.character(job$pid)]] <- list(job = job, id = row$Task_ID, start = Sys.time())
    }
    Sys.sleep(0.25)
    for (key in names(running)) {
      e <- running[[key]]
      res <- mccollect(e$job, wait = FALSE)
      el <- as.numeric(difftime(Sys.time(), e$start, units = "secs"))
      if (!is.null(res)) {
        r <- res[[1]]
        if (inherits(r, "try-error")) { n_bad <- n_bad + 1; logf("ERROR task %d after %.1fs: %s", e$id, el, trimws(as.character(r))) }
        else { n_ok <- n_ok + 1; logf("done task %d in %.1fs", e$id, el) }
        running[[key]] <- NULL
      } else if (el > outer_timeout) {
        # kill the task AND any MIMOSA2 fit it started (its child process)
        try(system(sprintf("pkill -KILL -P %s", key), ignore.stdout = TRUE, ignore.stderr = TRUE), silent = TRUE)
        tools::pskill(as.integer(key), tools::SIGKILL)
        suppressWarnings(mccollect(e$job, wait = TRUE, timeout = 5))
        n_bad <- n_bad + 1
        logf("KILLED task %d after %.0fs (outer limit)", e$id, el)
        running[[key]] <- NULL
      }
    }
    if ((n_ok + n_bad) - last_report >= 500) {
      rate <- (n_ok + n_bad) / as.numeric(difftime(Sys.time(), t_start, units = "hours"))
      message(sprintf("[%s] %d/%d finished (%d failed). %.0f tasks/hour; ~%.1f h left",
                      format(Sys.time(), "%H:%M"), n_ok + n_bad, nrow(todo), n_bad, rate,
                      (nrow(todo) - n_ok - n_bad) / rate))
      last_report <- n_ok + n_bad
    }
  }
  logf("finished: %d ok, %d failed", n_ok, n_bad)
  mark_finished(study, chunk_id, n_chunks)
  message(sprintf("[%s] finished: %d ok, %d failed (see %s)", study, n_ok, n_bad, log_file))
  invisible(TRUE)
}

# A chunk that reaches the end of its task list writes a marker file. If every
# chunk of this profile has one, the run is complete and combine_study() treats
# any missing task as a genuine failure; otherwise it treats the unfinished
# last repetition as "not run yet".
mark_finished <- function(study, chunk_id, n_chunks, profile = PROFILE) {
  writeLines(as.character(n_chunks),
             file.path(OUT_DIR, study, sprintf("finished_%s_chunk%02d.txt", profile, chunk_id)))
}
run_is_finished <- function(study, profile = PROFILE) {
  f <- list.files(file.path(OUT_DIR, study), pattern = sprintf("^finished_%s_chunk[0-9]+\\.txt$", profile), full.names = TRUE)
  if (length(f) == 0) return(FALSE)
  n_chunks <- max(as.integer(vapply(f, function(x) readLines(x, n = 1), "")))
  length(f) >= n_chunks
}

# ---- Combine task files into four tables ------------------------------------
# Writes _simulations/<study>/<study>_results.rds containing:
#   design    : one row per planned dataset in this profile (+ Seed)
#   datasets  : realised DGM quantities per dataset (n responders, true Delta...)
#   fits      : one row per MIMOSA2 fit (status, time, iterations, ...)
#   estimates : one row per dataset x method x rule x alpha (TP, FP, TN, FN, AUC)
#   subjects  : one row per subject (truth, scores, q-values)
# Planned datasets with no task file are added with Status = "task_failed".
combine_study <- function(study, profile = PROFILE, keep_subjects = TRUE) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("combine_study() needs dplyr")
  design <- readRDS(file.path(OUT_DIR, study, "design_full.rds"))
  plan <- select_profile(design, study, profile)
  have <- file.exists(task_file(study, plan$Task_ID))
  # If the run is unfinished (time limit), tasks run in repetition order, so
  # the last repetition reached is partly done and later ones not started.
  # Those are NOT failures: drop them, so every scenario keeps the same nsim.
  # Missing tasks in earlier repetitions are genuine task failures.
  if (any(have) && !all(have) && !run_is_finished(study, profile)) {
    r_last <- max(plan$Rep[have])
    r_keep <- if (all(have[plan$Rep == r_last])) r_last else r_last - 1
    if (r_keep < max(plan$Rep)) {
      message(sprintf("[%s] run not finished: using repetitions 1-%d (planned %d)", study, r_keep, max(plan$Rep)))
      plan <- plan[plan$Rep <= r_keep, , drop = FALSE]
      have <- file.exists(task_file(study, plan$Task_ID))
    }
  }
  files <- task_file(study, plan$Task_ID)
  out_file <- file.path(OUT_DIR, study, sprintf("%s_results_%s.rds", study, profile))
  # [CHANGE v3.1] never replace a results file by an (almost) empty one
  if (sum(have) == 0)
    stop(sprintf("[%s] no task files found in %s; refusing to (over)write %s", study, task_dir(study), out_file))
  message(sprintf("[%s] combining %d of %d planned tasks", study, sum(have), nrow(plan)))
  res <- lapply(files[have], readRDS)
  pick <- function(nm) dplyr::bind_rows(lapply(res, `[[`, nm))
  out <- list(design = plan,
              datasets = pick("datasets"),
              fits = pick("fits"),
              estimates = pick("estimates"),
              subjects = if (keep_subjects) pick("subjects") else NULL,
              task_time = data.frame(Task_ID = vapply(res, function(r) as.numeric(r$Task_ID), 1),
                                     Task_time = vapply(res, function(r) as.numeric(r$task_time), 1)),
              run_info = lapply(list.files(file.path(OUT_DIR, study), "^run_info_chunk", full.names = TRUE), readRDS),
              profile = profile, combined = Sys.time())
  missing <- plan$Task_ID[!have]
  if (length(missing) > 0) {
    warning(sprintf("%d planned tasks have no result file; recorded as task_failed", length(missing)))
    out$fits <- dplyr::bind_rows(out$fits, data.frame(Task_ID = missing, Method = "ALL", Status = "task_failed",
                                                      stringsAsFactors = FALSE))
  }
  # [CHANGE v3.1] keep the previous results file instead of overwriting it
  if (file.exists(out_file)) {
    prev <- sub("\\.rds$", "_previous.rds", out_file)
    file.rename(out_file, prev)
    message(sprintf("[%s] previous results kept as %s", study, prev))
  }
  saveRDS(out, out_file)
  invisible(out)
}

# ---- Re-create one dataset exactly (Morris 4.1) ------------------------------
# e.g. bad <- subset(res$fits, Status == "timeout")[1, ]
#      row <- subset(res$design, Task_ID == bad$Task_ID)
#      sim <- regenerate_dataset(row, simulate_fun)   # then inspect / refit
regenerate_dataset <- function(row, simulate_fun) {
  set_task_seed(row$Seed)
  simulate_fun(row)
}
