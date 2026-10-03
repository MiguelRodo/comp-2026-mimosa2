# MIMOSA2 Simulation Study
# Isabella Lethbridge and Tayyeb Abrahams 
# October 2026
# ==============================================================================
# 11_smoke_test.R: Run first on cluster 
# ==============================================================================
#   cd mimosa2_sim_v3
#   Rscript sims/11_smoke_test.R          (on a compute node, e.g. in an
#                                          interactive job with 4+ cores)
# It checks, in order:
#   1. Required packages are installed
#   2. What the MIMOSA2 fit object contains (names(fit)),
#      that the posterior probabilities are valid, 
#      That getResponse() agrees with our own Bayesian-FDR q-values
#   3. The new DiD_wald() reproduces the version-2 GLM (DiD_GLM) exactly
#      wherever the GLM converged with a non-negative estimate
#   4. The time limit works (a sleeping job is killed; errors are caught)
#   5. How long MIMOSA2 takes at P = 10 ... 100 on THIS machine, and from that
#      the projected wall-clock time of every study for the standard and
#      extended profiles on N_WORKERS cores
#   6. The whole pipeline end to end 
#      (all five studies + all analysis scripts) in the tiny "smoke" profile, 
#      writing to _smoke/ so nothing real is touched.
# A summary of PASS / FAIL lines is printed at the end.
# =============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
for (f in c("R/config.R", "R/dgm.R", "R/methods.R", "R/runner.R", "R/scenarios.R")) source(f)
results <- character(0)
record <- function(ok, what) {
  tag <- if (isTRUE(ok)) "PASS" else "FAIL"
  results <<- c(results, sprintf("%s  %s", tag, what)); message(tag, "  ", what)
}

# ------------------------------ 1. packages -----------------------------------
message("\n== 1. Packages ==")
need <- c("MIMOSA2", "dplyr", "tidyr", "purrr", "ggplot2", "knitr")
nice <- c("plotROC", "scales")
for (p in need) record(requireNamespace(p, quietly = TRUE), paste("package installed:", p))
for (p in nice) if (!requireNamespace(p, quietly = TRUE)) message("note: optional package ", p, " not installed")
if (!requireNamespace("MIMOSA2", quietly = TRUE)) stop("MIMOSA2 is not installed; nothing else can be tested.")
suppressPackageStartupMessages(library(MIMOSA2))
message("MIMOSA2 version: ", as.character(packageVersion("MIMOSA2")), " | R ", R.version.string)

# -------------------------- 2. MIMOSA2 fit object -----------------------------
message("\n== 2. MIMOSA2 fit object ==")
RNGkind("L'Ecuyer-CMRG"); set.seed(11)
row <- data.frame(Effect = 1e-3, Phi = 2000, P = 50, Distribution = "Beta", Res_prop = "Prop_0.50", Cell_range = "Medium")
sim <- simulate_baseline(row)
fit <- MIMOSA2(Ntot = sim$Ntot, ns1 = sim$ns1, nu1 = sim$nu1, ns0 = sim$ns0, nu0 = sim$nu0, maxit = MAXIT, verbose = FALSE)
message("names(fit): ", paste(names(fit), collapse = ", "))
str(fit, max.level = 1, give.attr = FALSE)
z <- as.matrix(fit$z)
record(nrow(z) == 50 && ncol(z) >= 8, sprintf("fit$z is P x K (got %d x %d)", nrow(z), ncol(z)))
record(all(abs(rowSums(z) - 1) < 1e-6), "rows of fit$z sum to 1 (posterior component probabilities)")
prob <- rowSums(z[, 1:4])
message(sprintf("mean P(responder) among true responders %.3f, non-responders %.3f",
                mean(prob[is_responder(sim$truth) == 1]), mean(prob[is_responder(sim$truth) == 0])))
record(mean(prob[is_responder(sim$truth) == 1]) > mean(prob[is_responder(sim$truth) == 0]),
       "components 1-4 of fit$z are the responder components (responders score higher)")
message("fit_iterations(fit) = ", fit_iterations(fit),
        "  (NA means the fit object has no iteration count; Iter will be NA in the results)")
q <- bayes_fdr_q(prob)
for (a in ALPHAS) {
  gr <- tryCatch(as.logical(getResponse(fit, threshold = a)), error = function(e) NULL)
  record(!is.null(gr) && length(gr) == 50, sprintf("getResponse(fit, %g) returns P logical calls", a))
  if (!is.null(gr)) record(all(gr == (q < a)), sprintf("getResponse(fit, %g) == our Bayesian-FDR rule q < %g", a, a))
}
if (exists("getFDR", envir = asNamespace("MIMOSA2"))) {
  gf <- tryCatch(MIMOSA2::getFDR(fit), error = function(e) NULL)
  if (!is.null(gf) && length(gf) == 50) {
    message(sprintf("max |getFDR - our q| = %.2e", max(abs(gf - q))))
    record(max(abs(gf - q)) < 1e-6, "getFDR(fit) equals our q-values")
  } else message("note: getFDR(fit) did not return a vector of length P")
}
f2 <- fit_mimosa2(sim)
record(f2$status == "ok" && isTRUE(f2$getresponse_agree), "fit_mimosa2() wrapper: status ok and calls agree")

# ------------------------------- 3. DiD ---------------------------------------
message("\n== 3. DiD_wald() vs version-2 DiD_GLM() ==")
for (cell in c("Medium", "Low")) {
  row$Cell_range <- cell
  sim <- simulate_baseline(row)
  new <- DiD_wald(sim$Ntot, sim$ns1, sim$nu1, sim$ns0, sim$nu0)
  old <- DiD_GLM_legacy(sim$Ntot, sim$ns1, sim$nu1, sim$ns0, sim$nu0)
  tied <- abs(old - 0.05) < 1e-9                  # 1 - 0.95 is not exactly 0.05 in floating point
  allpos <- sim$ns1 > 0 & sim$nu1 > 0 & sim$ns0 > 0 & sim$nu0 > 0
  same <- !tied & allpos                           # compare where the GLM is not at a boundary
  message(sprintf("%s cells: GLM gave the tie value 0.05 to %d of %d subjects (negative estimate or no convergence); %d subjects have a zero count",
                  cell, sum(tied), length(old), sum(!allpos)))
  record(max(abs(new$DiD_score[same] - old[same])) < 1e-4,
         sprintf("DiD_wald score == DiD_GLM where GLM converged with est >= 0 (%s cells; max diff %.1e)",
                 cell, max(abs(new$DiD_score[same] - old[same]))))
}

# ------------------------------ 4. time limit ---------------------------------
message("\n== 4. Time limit ==")
r1 <- run_with_timeout(function() { Sys.sleep(5); 1 }, timeout = 1)
record(r1$status == "timeout" && r1$elapsed < 4, "a 5-second job is killed after 1 second")
r2 <- run_with_timeout(function() stop("boom"), timeout = 5)
record(r2$status == "error" && grepl("boom", r2$msg), "an error inside the job is caught and recorded")
r3 <- run_with_timeout(function() 42, timeout = 5)
record(r3$status == "ok" && r3$value == 42, "a normal job returns its value")

# --------------------------- 5. timing benchmark ------------------------------
message("\n== 5. Timing benchmark (this takes a few minutes) ==")
bench <- list()
for (P in c(10, 20, 50, 100)) for (cell in c("Medium", "Low")) for (rep in 1:2) {
  rw <- data.frame(Effect = 2.5e-4, Phi = 2000, P = P, Distribution = "Beta", Res_prop = "Prop_0.50", Cell_range = cell)
  s <- simulate_baseline(rw)
  t <- system.time(MIMOSA2(Ntot = s$Ntot, ns1 = s$ns1, nu1 = s$nu1, ns0 = s$ns0, nu0 = s$nu0, maxit = MAXIT, verbose = FALSE))[["elapsed"]]
  bench[[length(bench) + 1]] <- data.frame(P = P, Cell = cell, maxit = MAXIT, sec = t)
}
s <- simulate_baseline(data.frame(Effect = 2.5e-4, Phi = 2000, P = 50, Distribution = "Beta", Res_prop = "Prop_0.50", Cell_range = "Medium"))
t100 <- system.time(MIMOSA2(Ntot = s$Ntot, ns1 = s$ns1, nu1 = s$nu1, ns0 = s$ns0, nu0 = s$nu0, maxit = 100, verbose = FALSE))[["elapsed"]]
bench <- do.call(rbind, bench)
print(aggregate(sec ~ P + Cell, bench, mean))
tm <- lm(sec ~ P, data = bench)
fit_time <- function(P) pmax(predict(tm, newdata = data.frame(P = P)), 0.5) + 0.5   # + DiD / overhead
t30_50 <- fit_time(50)
ratio100 <- max(1, t100 / t30_50)
message(sprintf("time per fit ~ %.2f + %.3f * P seconds; maxit=100 takes %.1fx maxit=30 at P = 50",
                coef(tm)[1], coef(tm)[2], ratio100))
project <- function(study, scen, prof, task_sec) {
  des <- build_design(study, scen)
  sel <- select_profile(des, study, prof)
  sum(task_sec(sel)) / 3600
}
proj <- list()
for (prof in c("standard", "extended")) {
  b <- project("baseline", scenarios_baseline(), prof, function(d) fit_time(d$P))
  p <- project("prior", scenarios_prior(), prof, function(d) fit_time(d$P))
  h <- project("heterogeneity", scenarios_heterogeneity(), prof, function(d) fit_time(d$P_focal + d$P_other) + fit_time(d$P_focal))
  i <- project("imbalance", scenarios_imbalance(), prof, function(d) fit_time(d$P))
  mx_des <- build_design("baseline", scenarios_baseline())
  mx_des <- mx_des[maxit_check_subset(mx_des) & mx_des$Rep <= NSIM[[prof]][["maxit_check"]], ]
  m <- sum(fit_time(mx_des$P) * (1 + ratio100)) / 3600
  proj[[prof]] <- data.frame(Profile = prof, Study = c("baseline", "prior", "heterogeneity", "imbalance", "maxit_check", "TOTAL"),
                             Core_hours = round(c(b, p, h, i, m, b + p + h + i + m), 1))
}
proj <- do.call(rbind, proj)
proj$Hours_on_20_cores <- round(proj$Core_hours / 20, 1)
proj$Hours_on_64_cores <- round(proj$Core_hours / 64, 1)
proj$Hours_on_this_N_WORKERS <- round(proj$Core_hours / N_WORKERS, 1)
message("Projected run time (assumes this node is as fast as the cluster nodes):")
print(proj, row.names = FALSE)
write.csv(proj, file.path(TAB_DIR, "smoke_projected_runtime.csv"), row.names = FALSE)

# ------------------------ 6. end-to-end smoke run -----------------------------
message("\n== 6. End-to-end run of every script in the 'smoke' profile (writes to _smoke/) ==")
env <- c("SIM_PROFILE=smoke", "MIMOSA2_OUT_DIR=_smoke/_simulations", "MIMOSA2_FIG_DIR=_smoke/_fig",
         "MIMOSA2_TAB_DIR=_smoke/_tables", sprintf("N_WORKERS=%d", max(2, min(N_WORKERS, 8))))
scripts <- c("sims/20_sim_baseline.R", "sims/21_sim_prior.R", "sims/22_sim_heterogeneity.R",
             "sims/23_sim_imbalance.R", "sims/24_sim_maxit_check.R",
             "analysis/30_performance_tables.R", "analysis/31_plots_baseline.R",
             "analysis/32_plots_prior.R", "analysis/33_plots_heterogeneity.R",
             "analysis/34_plots_imbalance.R", "analysis/35_maxit_check.R")
rscript <- file.path(R.home("bin"), "Rscript")
for (sc in scripts) {
  logf <- file.path("_smoke", paste0(basename(sc), ".log"))
  dir.create("_smoke", showWarnings = FALSE)
  t <- system.time(st <- system2(rscript, sc, env = env, stdout = logf, stderr = logf))[["elapsed"]]
  record(st == 0, sprintf("%s ran without error in %.0fs (log: %s)", sc, t, logf))
}
for (s in c("baseline", "prior", "heterogeneity", "imbalance", "maxit_check")) {
  f <- file.path("_smoke/_simulations", s, sprintf("%s_results_smoke.rds", s))
  if (file.exists(f)) {
    r <- readRDS(f)
    nok <- sum(r$fits$Status == "ok")
    record(nrow(r$design) > 0 && nrow(r$estimates) > 0 && nok > 0,
           sprintf("%s: %d datasets, %d MIMOSA2 fits ok of %d", s, nrow(r$design), nok, nrow(r$fits)))
  } else record(FALSE, paste(s, ": no combined results file"))
}

message("\n================ SMOKE TEST SUMMARY ================")
cat(results, sep = "\n")
if (any(grepl("^FAIL", results))) {
  message("\nSome checks FAILED. Look at the lines above and the logs in _smoke/ before running the full study.")
  quit(status = 1)
} else message("\nAll checks passed. Delete the _smoke folder and start the full runs.")
