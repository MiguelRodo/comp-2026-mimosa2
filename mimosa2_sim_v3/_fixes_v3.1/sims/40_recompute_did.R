# =============================================================================
# 40_recompute_did.R : re-compute the DiD comparator with the corrected
#                      (continuity-corrected) variance, WITHOUT refitting MIMOSA2.
# NEW in v3.1 (6 Oct 2026).
# =============================================================================
# Why: in v3.0 DiD_wald() used p_hat (1 - p_hat) / N, so an assay with a zero
# count added no variance. See R/methods.R and R/config.R (DID_CC).
#
# How: the DiD needs only the cell counts, which are not stored, but every
# dataset can be re-created exactly from its stored random-number stream
# (runner.R, regenerate_dataset). For every dataset this script
#   1. regenerates the counts,
#   2. CHECKS that the uncorrected DiD estimate equals the one stored in the
#      results file (this is also a full reproducibility check, Morris 4.1),
#   3. recomputes DiD with cc = DID_CC and rebuilds the DiD rows of
#      res$estimates and the DiD columns of res$subjects.
# The MIMOSA2 results are not touched. The old DiD estimates are kept in
# res$estimates_did_v30, and the original file is copied to
# <study>_results_<profile>_before_didfix.rds before anything is written.
#
# Usage (from the mimosa2_sim_v3 folder; ~5-15 min on 20 cores, all studies):
#   SIM_PROFILE=standard Rscript sims/40_recompute_did.R baseline prior heterogeneity imbalance maxit_check
# On the cluster:  apptainer-rscript -f mimosa2 -- 'source("sims/40_recompute_did.R")'
# (no arguments = all five studies). Then re-run analysis/30-35.
# =============================================================================
source("sims/_header.R")
studies <- commandArgs(trailingOnly = TRUE)
if (length(studies) == 0 && nzchar(Sys.getenv("STUDIES"))) studies <- strsplit(Sys.getenv("STUDIES"), "[ ,]+")[[1]]
if (length(studies) == 0) studies <- c("baseline", "prior", "heterogeneity", "imbalance", "maxit_check")

SIMULATOR <- list(baseline = simulate_baseline, prior = simulate_prior,
                  heterogeneity = simulate_heterogeneity, imbalance = simulate_imbalance,
                  maxit_check = simulate_baseline)          # maxit_check re-uses the baseline datasets
cores <- if (.Platform$OS.type == "windows") 1L else as.integer(N_WORKERS)

one_task <- function(row, study) {
  sim <- regenerate_dataset(row, SIMULATOR[[study]])
  P <- length(sim$ns1)
  eval_idx <- if (study == "heterogeneity") seq_len(row$P_focal) else seq_len(P)
  old <- DiD_wald(sim$Ntot, sim$ns1, sim$nu1, sim$ns0, sim$nu0, cc = 0)        # v3.0 statistic
  new <- DiD_wald(sim$Ntot, sim$ns1, sim$nu1, sim$ns0, sim$nu0, cc = DID_CC)   # v3.1 statistic
  truth <- is_responder(sim$truth)[eval_idx]
  p_eval <- new$DiD_p[eval_idx]
  calls_raw <- setNames(lapply(ALPHAS, function(a) p_eval <= a), paste0("a", ALPHAS))
  p_bh <- p.adjust(p_eval, method = "BH")
  calls_bh  <- setNames(lapply(ALPHAS, function(a) p_bh <= a), paste0("a", ALPHAS))
  est <- rbind(estimate_rows(row$Task_ID, "DiD", "unadjusted", truth, new$DiD_z[eval_idx], calls_raw, "ok"),
               estimate_rows(row$Task_ID, "DiD", "BH",         truth, new$DiD_z[eval_idx], calls_bh,  "ok"))
  list(est = est,
       subj = data.frame(Task_ID = row$Task_ID, Subject_id = seq_len(P), old_est = old$DiD_est,
                         new, DiD_z_v30 = old$DiD_z, DiD_p_v30 = old$DiD_p))
}

for (st in studies) {
  f <- file.path(OUT_DIR, st, sprintf("%s_results_%s.rds", st, PROFILE))
  if (!file.exists(f)) { message("skip ", st, ": ", f, " not found"); next }
  message(sprintf("\n==== %s ====", st))
  res <- readRDS(f)
  if (!is.null(res$did_fix)) { message("already corrected (cc = ", res$did_fix$cc, "); skipped"); next }
  # a study run with the v3.1 code already has the corrected DiD (run_info records did_cc)
  cc_run <- vapply(res$run_info, function(r) if (is.null(r$did_cc)) 0 else as.numeric(r$did_cc), 0)
  if (length(cc_run) && all(cc_run > 0)) { message("run with the v3.1 code (did_cc = ", cc_run[1], "): DiD already corrected; skipped"); next }
  bak <- sub("\\.rds$", "_before_didfix.rds", f)
  if (!file.exists(bak)) { stopifnot(file.copy(f, bak)); message("backup: ", bak) }
  if (st == "prior") prewarm_calibration(res$design$Distribution, res$design$Effect, res$design$Phi)

  rows <- split(res$design, seq_len(nrow(res$design)))
  t0 <- Sys.time()
  out <- parallel::mclapply(rows, one_task, study = st, mc.cores = cores, mc.set.seed = FALSE, mc.preschedule = TRUE)
  bad <- vapply(out, function(x) inherits(x, "try-error") || is.null(x$est), TRUE)
  if (any(bad)) stop(sprintf("[%s] %d datasets could not be regenerated, e.g. %s", st, sum(bad), as.character(out[[which(bad)[1]]])))
  message(sprintf("[%s] regenerated %d datasets in %.1f min", st, length(out), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  new_subj <- as.data.frame(dplyr::bind_rows(lapply(out, `[[`, "subj")))
  new_est  <- as.data.frame(dplyr::bind_rows(lapply(out, `[[`, "est")))

  # ---- reproducibility check: regenerated data must give the stored DiD ----
  key_old <- paste(res$subjects$Task_ID, res$subjects$Subject_id)
  m <- match(key_old, paste(new_subj$Task_ID, new_subj$Subject_id))
  if (anyNA(m)) stop(sprintf("[%s] %d stored subjects were not regenerated", st, sum(is.na(m))))
  d <- abs(new_subj$old_est[m] - res$subjects$DiD_est)
  message(sprintf("[%s] reproducibility: max |DiD_est regenerated - stored| = %.3g over %d subjects", st, max(d), length(d)))
  if (max(d) > 1e-12) stop(sprintf("[%s] regenerated datasets differ from the stored ones: NOT written", st))

  # ---- write the corrected DiD into the results ----
  res$estimates_did_v30 <- subset(res$estimates, Method == "DiD")
  res$estimates <- as.data.frame(dplyr::bind_rows(subset(res$estimates, Method != "DiD"), new_est[, names(res$estimates)]))
  res$estimates <- res$estimates[order(res$estimates$Task_ID), ]
  for (cl in c("DiD_se", "DiD_z", "DiD_p", "DiD_score", "DiD_z_v30", "DiD_p_v30")) res$subjects[[cl]] <- new_subj[[cl]][m]
  res$subjects$DiD_GLM_prob <- res$subjects$DiD_score            # legacy column follows the new score
  res$did_fix <- list(cc = DID_CC, date = Sys.time(), max_repro_diff = max(d),
                      note = "DiD variance continuity-corrected; MIMOSA2 untouched; old DiD rows in estimates_did_v30")
  tmp <- paste0(f, ".tmp"); saveRDS(res, tmp); file.rename(tmp, f)
  # quick summary: false-positive rate of unadjusted DiD at 5% among non-responders, v3.0 vs v3.1
  s <- res$subjects[res$subjects$Eval & res$subjects$Truth == 0, ]
  message(sprintf("[%s] non-responders called by unadjusted DiD at 5%%: v3.0 %.3f -> v3.1 %.3f",
                  st, mean(s$DiD_p_v30 <= 0.05), mean(s$DiD_p <= 0.05)))
  message("written: ", f)
}
message("\nDone. Now re-run analysis/30_performance_tables.R and analysis/31-35.")
