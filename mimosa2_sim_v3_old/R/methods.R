# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# methods.R: MIMOSA2 vs DiD 
# ==============================================================================
# 1. Difference-in-Differences (DiD)
# DiD_wald() replaces DiD_GLM(): 
# The statistic is the same
#    The per-subject binomial GLM with identity link and the Time x Stimulation
#    interaction is saturated (4 parameters, 4 observations), so its MLE is the
#    observed difference-in-differences and its Wald variance is sum over the 4
#    assays of p_hat (1 - p_hat) / N.
# Computing this directly gives the identical estimate and z statistic when
# the GLM converges (checked in sims/11_smoke_test.R) and fixes two problems
# in DiD_GLM():
#   (a) every subject with a NEGATIVE estimate, and every subject whose GLM did
#       not converge, was given the same score 0.05. These ties distort the
#       ROC curve and AUC, and non-convergence (common when a count is zero,
#       i.e. at low cell counts) was silently treated as "confident
#       non-responder".
#   (b) it ran one glm() per subject (slow).
# Scores: score = Phi(z) = 1 - one-sided p-value (as before for est >= 0, and
# now also defined for est < 0). 
# Ranking for ROC/AUC uses z itself, which has no ties even when Phi(z) rounds to 1
# If all four observed proportions are 0 (SE = 0 and estimate 0), z = 0.
# ==============================================================================
# Function to compute wald intervals for DiD:
DiD_wald <- function(Ntot, ns1, nu1, ns0, nu0) {
  N_s1 <- Ntot[, "ns1"]; N_u1 <- Ntot[, "nu1"]; N_s0 <- Ntot[, "ns0"]; N_u0 <- Ntot[, "nu0"]
  p_s1 <- ns1 / N_s1; p_u1 <- nu1 / N_u1; p_s0 <- ns0 / N_s0; p_u0 <- nu0 / N_u0
  est <- (p_s1 - p_u1) - (p_s0 - p_u0)
  se  <- sqrt(p_s1 * (1 - p_s1) / N_s1 + p_u1 * (1 - p_u1) / N_u1 +
              p_s0 * (1 - p_s0) / N_s0 + p_u0 * (1 - p_u0) / N_u0)
  z <- est / se
  z[se == 0 & est == 0] <- 0
  data.frame(DiD_est = est, DiD_se = se, DiD_z = z,
             DiD_p = pnorm(z, lower.tail = FALSE),   # one-sided H1: Delta > 0
             DiD_score = pnorm(z))                   
}

# The version-2 function, kept ONLY so the smoke test can show that DiD_wald()
# reproduces it (do not use for the study).
DiD_GLM_legacy <- function(Ntot, ns1, nu1, ns0, nu0) {
  P <- nrow(Ntot)
  out <- numeric(P)
  for (i in 1:P) {
    df <- data.frame(Time = factor(c("Active", "Active", "Baseline", "Baseline"), levels = c("Baseline", "Active")),
                     Stim = factor(c("Stimulated", "Unstimulated", "Stimulated", "Unstimulated"), levels = c("Unstimulated", "Stimulated")),
                     Pos  = c(ns1[i], nu1[i], ns0[i], nu0[i]),
                     Tot  = c(Ntot[i, "ns1"], Ntot[i, "nu1"], Ntot[i, "ns0"], Ntot[i, "nu0"]))
    df$Neg <- df$Tot - df$Pos
    pr <- pmax(pmin(df$Pos / df$Tot, 1 - 1e-5), 1e-5)
    st <- c(pr[4], pr[2] - pr[4], pr[3] - pr[4], (pr[1] - pr[2]) - (pr[3] - pr[4]))
    fit <- suppressWarnings(tryCatch(glm(cbind(Pos, Neg) ~ Time * Stim, family = binomial(link = "identity"),
                                         data = df, start = st, control = glm.control(maxit = 200, epsilon = 1e-8)),
                                     error = function(e) NULL))
    if (is.null(fit) || !fit$converged) { out[i] <- 0.95; next }
    cm <- summary(fit)$coefficients
    nm <- "TimeActive:StimStimulated"
    if (nm %in% rownames(cm) && !is.na(cm[nm, 1]) && cm[nm, 1] >= 0) out[i] <- cm[nm, 4] / 2 else out[i] <- 0.95
  }
  1 - out
}

# ==============================================================================
# 2. Hard time limit for one function call (Linux / macOS: fork-based)
# ==============================================================================
# R.utils::withTimeout() cannot interrupt MIMOSA2's compiled optimiser, which is 
# why version 2 moved to mcparallel(). 
# Version 2 killed the WHOLE task on timeout, so the DiD results for that dataset 
# were also lost and the dataset could not be identified. 
# Now each MIMOSA2 fit runs in its own child process; if it runs past the limit 
# the child is killed and the fit is recorded as "timeout", and the rest of the 
# task (DiD, bookkeeping) still completes. 
# On Windows (no fork) the call simply runs without a limit
# ==============================================================================
run_with_timeout <- function(fun, timeout, poll = 0.1) {
  t0 <- Sys.time()
  if (.Platform$OS.type == "windows") {
    v <- tryCatch(list(status = "ok", value = fun(), msg = NA_character_),
                  error = function(e) list(status = "error", value = NULL, msg = conditionMessage(e)))
    v$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    return(v)
  }
  job <- parallel::mcparallel(list(value = fun()), silent = TRUE, mc.set.seed = FALSE)
  repeat {
    res <- parallel::mccollect(job, wait = FALSE, timeout = poll)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(res)) {
      r <- res[[1]]
      if (inherits(r, "try-error"))
        return(list(status = "error", value = NULL, msg = trimws(as.character(r)), elapsed = elapsed))
      if (is.null(r) || !is.list(r))
        return(list(status = "error", value = NULL, msg = "child process returned no result (crashed?)", elapsed = elapsed))
      return(list(status = "ok", value = r$value, msg = NA_character_, elapsed = elapsed))
    }
    if (elapsed > timeout) {
      tools::pskill(job$pid, tools::SIGKILL)
      suppressWarnings(parallel::mccollect(job, wait = TRUE, timeout = 5))
      return(list(status = "timeout", value = NULL,
                  msg = sprintf("killed after %.0f s", elapsed), elapsed = elapsed))
    }
  }
}

# ==============================================================================
# 3. MIMOSA2 fit + everything we need from it
# ==============================================================================
# Bayesian FDR q-values (direct posterior probability approach, Newton et al.
# 2004): order subjects by decreasing P(responder)
# The q-value of the subject in position m is the mean of (1 - P(responder)) 
# over the top m subjects.
# A subject is called a responder at level alpha if q < alpha (getResponse).
bayes_fdr_q <- function(prob) {
  o <- order(prob, decreasing = TRUE)
  q <- numeric(length(prob))
  q[o] <- cumsum(1 - prob[o]) / seq_along(o)
  q
}

# NOTE: 
# Version 2 stored Iterations = length(fit$inds), which is P x 11
# (the size of an indicator matrix), not the number of EM iterations.
# MIMOSA2 (0.99.x) does not return its iteration count or a convergence flag
# (the fit contains z, inds, pi_est, thetahat, ps1_hat, ps0_hat, pu1_hat,
# pu0_hat, Ntot and the counts), so Iter / Hit_maxit will normally be NA.
# This function is kept in case a later package version adds it; whether
# maxit = 30 is enough is assessed directly by Study 5 (maxit 30 vs 100).
# sims/11_smoke_test.R prints names(fit) so you can see what is available.
fit_iterations <- function(fit) {
  if (!is.list(fit)) return(NA_real_)
  for (nm in c("iter", "iterations", "niter", "n_iter", "iteration")) {
    v <- fit[[nm]]
    if (is.numeric(v) && length(v) == 1) return(as.numeric(v))
  }
  for (nm in c("ll", "loglik", "logLik", "LL", "l", "lls", "loglikelihood")) {
    v <- fit[[nm]]
    if (is.numeric(v) && length(v) > 1) return(as.numeric(length(v)))
  }
  NA_real_
}

fit_mimosa2 <- function(sim, idx = seq_along(sim$ns1), maxit = MAXIT, timeout = FIT_TIMEOUT) {
  P <- length(idx)
  out <- list(prob = rep(NA_real_, P), q = rep(NA_real_, P),
              calls = setNames(vector("list", length(ALPHAS)), paste0("a", ALPHAS)),
              status = NA_character_, msg = NA_character_, time = NA_real_,
              iter = NA_real_, hit_maxit = NA, rho_hat = NA_real_,
              calls_source = NA_character_, getresponse_agree = NA)
  Ntot <- sim$Ntot[idx, , drop = FALSE]
  r <- run_with_timeout(function() {
    MIMOSA2::MIMOSA2(Ntot = Ntot, ns1 = sim$ns1[idx], nu1 = sim$nu1[idx],
                     ns0 = sim$ns0[idx], nu0 = sim$nu0[idx],
                     maxit = maxit, verbose = FALSE)
  }, timeout = timeout)
  out$status <- r$status; out$msg <- r$msg; out$time <- r$elapsed
  if (r$status != "ok") return(out)
  fit <- r$value
  z <- tryCatch(as.matrix(fit$z), error = function(e) NULL)
  if (is.null(z) || nrow(z) != P || ncol(z) < 4 || any(!is.finite(z[, 1:4]))) {
    out$status <- "error"; out$msg <- "fit$z missing, wrong size or not finite"
    return(out)
  }
  prob <- rowSums(z[, 1:4, drop = FALSE])     # P(responder): components 1-4 (as version 2)
  prob <- pmin(pmax(prob, 0), 1)
  out$prob <- prob
  out$q    <- bayes_fdr_q(prob)
  out$iter <- fit_iterations(fit)
  out$hit_maxit <- if (is.na(out$iter)) NA else out$iter >= maxit
  out$rho_hat <- mean(prob)
  # Calls: the package's getResponse()
  # If it is not available or fails, fall back to q-values (same rule).
  agree <- logical(0)
  for (a in ALPHAS) {
    own <- out$q < a  # MIMOSA2::getResponse is getFDR(fit) < threshold
    pk <- tryCatch(as.logical(MIMOSA2::getResponse(fit, threshold = a)), error = function(e) NULL)
    if (!is.null(pk) && length(pk) == P && !anyNA(pk)) {
      out$calls[[paste0("a", a)]] <- pk; out$calls_source <- "getResponse"
      agree <- c(agree, all(pk == own))
    } else {
      out$calls[[paste0("a", a)]] <- own; out$calls_source <- "own_q"
    }
  }
  out$getresponse_agree <- if (length(agree)) all(agree) else NA
  out
}

# ==============================================================================
# 4. Per-dataset performance 
#===============================================================================
# AUC for ONE dataset (Mann-Whitney form; ties count 1/2). 
# NA if the dataset has no responders or no non-responders.
auc_mw <- function(score, truth) {
  ok <- !is.na(score); score <- score[ok]; truth <- truth[ok]
  n1 <- sum(truth == 1); n0 <- sum(truth == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
confusion <- function(call, truth) {
  c(TP = sum(call & truth == 1), FP = sum(call & truth == 0),
    TN = sum(!call & truth == 0), FN = sum(!call & truth == 1))
}

# Build the per-dataset estimate rows for one method:
#    calls_by_alpha: named list (a0.01, a0.05) of logical vectors, NULL if failed
#    score: ranking score (NA if failed)
estimate_rows <- function(task_id, method, rule, truth, score, calls_by_alpha, status) {
  rows <- lapply(ALPHAS, function(a) {
    cl <- if (is.null(calls_by_alpha)) NULL else calls_by_alpha[[paste0("a", a)]]
    cf <- if (is.null(cl) || status != "ok") c(TP = NA, FP = NA, TN = NA, FN = NA) else confusion(cl, truth)
    data.frame(Task_ID = task_id, Method = method, Rule = rule, Alpha = a,
               Status = status,
               n_resp = sum(truth == 1), n_nonresp = sum(truth == 0),
               TP = cf[["TP"]], FP = cf[["FP"]], TN = cf[["TN"]], FN = cf[["FN"]],
               AUC = if (status == "ok") auc_mw(score, truth) else NA_real_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ==============================================================================
# 5. Analyse one simulated dataset with all methods
# ==============================================================================
#   sim            : output of simulate_MIMOSA2_alt_prior() (or bind_sims())
#   task_id        : global task ID
#   eval_idx       : subjects on which performance is measured
#   mimosa_subsets : named list; one MIMOSA2 fit per element, fitted to those
#                    subjects (e.g. list(MIMOSA2 = all) or, for the
#                    heterogeneity study, pooled and separate fits)
#   group          : optional subject-level label stored in the subject table
analyse_dataset <- function(sim, task_id, eval_idx = seq_along(sim$ns1),
                            mimosa_subsets = list(MIMOSA2 = seq_along(sim$ns1)),
                            group = NULL, maxit = MAXIT, timeout = FIT_TIMEOUT) {
  P <- length(sim$ns1)
  truth_all <- is_responder(sim$truth)
  truth <- truth_all[eval_idx]

  # DiD (always succeeds): 
  did <- DiD_wald(sim$Ntot, sim$ns1, sim$nu1, sim$ns0, sim$nu0)
  p_eval <- did$DiD_p[eval_idx]
  did_calls_raw <- setNames(lapply(ALPHAS, function(a) p_eval <= a), paste0("a", ALPHAS))
  p_bh <- p.adjust(p_eval, method = "BH")
  did_calls_bh  <- setNames(lapply(ALPHAS, function(a) p_bh <= a), paste0("a", ALPHAS))

  est <- list(
    estimate_rows(task_id, "DiD", "unadjusted", truth, did$DiD_z[eval_idx], did_calls_raw, "ok"),
    estimate_rows(task_id, "DiD", "BH",         truth, did$DiD_z[eval_idx], did_calls_bh,  "ok"))

  subj <- data.frame(Task_ID = task_id, Subject_id = seq_len(P),
                     Group = if (is.null(group)) NA_character_ else group,
                     Eval = seq_len(P) %in% eval_idx,
                     Profile = sim$truth, Truth = truth_all,
                     True_delta = sim$true_delta,
                     did, stringsAsFactors = FALSE)
  # Legacy columns kept so old plotting code still runs
  subj$DiD_GLM_prob <- subj$DiD_score
  prop_s <- sim$ns1 / sim$Ntot[, "ns1"]; prop_u <- sim$nu1 / sim$Ntot[, "nu1"]
  subj$Log2_FC <- log2((prop_s + 1e-5) / (prop_u + 1e-5))

  fits <- list()
  for (nm in names(mimosa_subsets)) {
    idx <- mimosa_subsets[[nm]]
    # maxit may be one number, or a named vector with one value per MIMOSA2 fit
    mx <- if (!is.null(names(maxit)) && nm %in% names(maxit)) maxit[[nm]] else maxit[[1]]
    f <- fit_mimosa2(sim, idx, maxit = mx, timeout = timeout)
    prob_full <- q_full <- rep(NA_real_, P)
    prob_full[idx] <- f$prob; q_full[idx] <- f$q
    subj[[paste0(nm, "_prob")]] <- prob_full
    subj[[paste0(nm, "_q")]]    <- q_full
    # evaluate on eval subjects that were part of this fit
    pos <- match(eval_idx, idx)
    if (anyNA(pos)) stop("eval_idx must be a subset of every MIMOSA2 subset")
    calls_eval <- if (f$status == "ok") lapply(f$calls, function(cl) cl[pos]) else NULL
    est[[length(est) + 1]] <- estimate_rows(task_id, nm, "BFDR", truth, f$prob[pos], calls_eval, f$status)
    fits[[nm]] <- data.frame(Task_ID = task_id, Method = nm, Status = f$status, Msg = f$msg,
                             Fit_time = f$time, Iter = f$iter, Hit_maxit = f$hit_maxit,
                             Rho_hat = f$rho_hat, Calls_source = f$calls_source,
                             GetResponse_agree = f$getresponse_agree, P_fit = length(idx), Maxit = mx,
                             stringsAsFactors = FALSE)
  }
  resp <- truth_all == 1
  datasets <- data.frame(Task_ID = task_id, P_total = P,
                         n_resp_total = sum(resp), rho_realised = mean(resp),
                         n_eval = length(eval_idx), n_resp_eval = sum(truth == 1),
                         mean_true_delta_resp = if (any(resp)) mean(sim$true_delta[resp]) else NA_real_,
                         median_true_delta_resp = if (any(resp)) median(sim$true_delta[resp]) else NA_real_,
                         n_zero_ns1 = sum(sim$ns1 == 0),
                         stringsAsFactors = FALSE)
  for (lab in PROFILE_LABELS) datasets[[paste0("n_", lab)]] <- sum(sim$truth == lab)
  list(datasets = datasets, fits = do.call(rbind, fits),
       estimates = do.call(rbind, est), subjects = subj)
}
