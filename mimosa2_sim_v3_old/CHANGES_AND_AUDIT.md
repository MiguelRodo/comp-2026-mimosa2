# Audit of the version-2 files against Morris et al. (2019), and what changed

Severity:

- **A**: the results were wrong or invalid.
- **B**: Morris process requirement not met.
- **C**: cosmetic or minor.

## 1. Compliance with Morris et al. (2019): before and after

| Morris recommendation (section) | Version 2 | Version 3 |
|---|---|---|
| Plan with ADEMP; report in that order (3, 6.1) | Partly. Estimand not defined; nsim not justified. | `STUDY_DESCRIPTIONS.md` gives ADEMP for every study. |
| DGM produces the intended data (3.2, pitfall 4) | No. Wrong responder proportions; prior families mismatched (A1, A6). | Fixed and checked by `10_dgm_checks.R`. |
| Factorial design, or justify non-factorial (3.2) | Baseline factorial; other studies one-factor-away with different base cases | Unchanged by your choice ("keep designs"). Justify each base case in the write-up (Section 4 below). |
| Define estimands / targets (3.3) | AUC target undefined (pooled across datasets) | Target = each subject's responder status; all measures defined per dataset (Section 0.3 of `STUDY_DESCRIPTIONS.md`). |
| Include methods with known properties (4.2) | None | Null scenarios (ρ = 0); DiD-BH as an FDR-controlling comparator. |
| Same simulated data for all methods (4.3) | Yes | Yes, and comparisons are now **paired**. |
| Set seed once; store RNG state of every repetition (4.1, Table 1) | No seed (baseline); `set.seed` with forked Mersenne-Twister (others), so not reproducible | One master seed per study; one L'Ecuyer stream per dataset; streams stored. |
| Separate streams when parallelising (4.1.1) | Only `future.seed = TRUE` (baseline) | Yes. |
| Check the method does not reset the seed (4.1.1) | Not checked | `10_dgm_checks.R`, section E. |
| Start small, build in checks (4.2) | No checks | `10_dgm_checks.R` + `11_smoke_test.R`. |
| Record failures as missing with reasons; failures are the first performance measure (4.2, 5.1) | Inconsistent time limits (600/180/120/none); killed tasks lost all results; failures not reported | One 600 s limit per fit; status and message stored; failure rate reported per scenario. |
| Choose nsim for a target MCSE (5.3, Table 1) | 30–50, no justification; saved results used fewer (10–24) | 100 (standard) / 200 (extended); MCSE argument + post-hoc precision report. |
| Report MCSE for every performance estimate (5.2, 6.2) | No (bootstrap ribbons only; invalid Hanley–McNeil CIs) | Every table "estimate (MCSE)", every figure ±1.96 MCSE. |
| Do not average over DGMs (5.2) | AUC averaged over P and ρ; ROC and FDR pooled over P | Every measure per scenario. |
| Consider measures jointly (5.2) | Only TPR (sensitivity) plotted | TPR, specificity, FDR, P(any FP), AUC. |
| Separate simulation and analysis scripts (Table 1) | Yes | Yes. |
| Publish code that reproduces the results (Table 1, 8) | Saved .Rdata did not come from the committed scripts | Combined results store the profile, seeds, `sessionInfo()` and MIMOSA2 version. |

## 2. File-by-file assessment of the version-2 files

### `Non_beta_simulations.R` → `R/dgm.R`, `R/methods.R`

**A1. Responder allocation is wrong.**

Version 2 computed `n = floor(P*pis)` and then added the missing subjects one at a time to the smallest component (R1, then R2, …). The realised proportion of responders was therefore not the nominal one:

| nominal ρ | P = 10 | P = 20 | P = 30 | P = 50 | P = 75 | P = 100 |
|---|---|---|---|---|---|---|
| 0.10 | 0.20 | 0.20 | 0.20 | 0.12 | 0.15 | 0.12 |
| 0.25 | **0.60** | 0.40 | 0.33 | 0.28 | 0.25 | 0.28 |
| 0.50 | 0.60 | 0.60 | 0.53 | 0.52 | 0.52 | 0.52 |
| 0.75 | **0.60** | 0.60 | 0.67 | 0.72 | 0.75 | 0.72 |
| 0.90 | 0.80 | 0.80 | 0.80 | 0.88 | 0.85 | 0.88 |

At P = 10 the "25%" and "75%" designs were identical. Version 3 rounds ρP once (0.10 → 0.10, 0.25 → 0.30, 0.75 → 0.80 at P = 10) and stores the realised value per dataset.

**A2. Logit-normal mean.**

The location was logit(μ), which makes μ the median. The mean was 1.28μ. Version 3 solves the location so that the mean is μ.

**A3. Variance matching only at μ₀ = 10⁻⁴.**

The EG (k = 113) and LN (φ = 2.05) parameters were reused at every condition mean. At the stimulated responder mean (1.1×10⁻³ for δ = 10⁻³) the SDs relative to the Beta reference were EG 2.3× and LN 3.3×. They were also too large at μ(s,0) = 3.5×10⁻⁴ (EG 1.55×, LN 1.9×). The prior-family comparison therefore confounded "shape" with "more between-subject variance". Version 3 calibrates at every mean. The `prior_v2_vs_v3` figure shows this, and `_tables/prior_v2_vs_v3_moments_at_MS1.csv` gives the numbers.

**B. Bimodal family applied inconsistently.**

p(s,0) was bimodal in R3 and NSR but Beta in R1, R4 and NR2. Version 3 applies the bimodal family to every stimulated draw.

**B. "Simplex" naming.**

The "simplex" sampler is an inverse-Gaussian distribution on the odds scale, not the Barndorff-Nielsen–Jørgensen simplex. Its parameter (10,000) happens to match the Beta variance at every mean, so it was kept. Rename it in the write-up.

**B. `DiD_GLM()` scoring.**

Every subject with a negative estimate, or whose GLM did not converge, got the same score of 0.05. This produced tied ranks in the ROC/AUC. Non-convergence, which is common when a count is 0 at low cell counts, was treated as a confident "non-responder".

Version 3 computes the same Wald statistic in closed form. The model is saturated, so the estimate is the observed DiD and the SE is √Σ p̂(1 − p̂)/N. The closed form never fails and uses z for ranking. The smoke test shows it equals the GLM wherever the GLM converged with a non-negative estimate and no zero counts (elsewhere the GLM either gave the tie value or sat on a boundary).

**C. Unused code removed.**

Unused families (uniform, odds-exponential, odds-gamma, unit-lognormal) and the abandoned shrinkage and bivariate EM comparators were removed. Keep the old file if you describe those comparators in an appendix.

**C. `sum(components) != 1`.**

This exact check can fail on floating-point sums. It now uses a tolerance.

### `Simulation.R` (baseline) → `sims/20_sim_baseline.R`

**A. The script could not have produced the saved results.**

- It puts `DiD_ASH_prob = did_glm_shrink_prob` into the results, and that object is never defined.
- `withTimeout()` is used without loading `R.utils`, so every fit would fail inside the `tryCatch`.

So `Simulation_2.0.Rdata` came from another version of the code. Re-running from the committed scripts fixes this (Morris Table 1).

**B. Seeds.**

There was no `set.seed()`. `future.seed = TRUE` alone gives independent streams but different ones on every run.

**B. Time limit.**

The comment says "45 seconds" but the code uses 600 s. Version 3 uses 600 s everywhere and records each fit's status.

**B. Iterations column.**

`Iterations = length(fit$inds)` stores P × 11, not EM iterations. MIMOSA2 does not return its iteration count, which is why Study 5 was added.

**B. Number of repetitions.**

nsim = 30 gave an MCSE of about 0.045 for TPR when SD ≈ 0.25. Version 3 uses 100 or 200.

### `Prior_simulations.R` → `sims/21_sim_prior.R`

**A. Results not comparable across families.** See A2 and A3 above.

**A. Saved results did not match the script.**

The count table in `Prior_simulation_plots.R` shows 500 rows per cell, which is 10 repetitions. `simulation_times.csv` has 1,919 fits, about 24 repetitions. The script asks for 200 repetitions × 300 scenarios = 60,000 tasks, about 33 h on 20 cores.

**B. RNG.**

`set.seed(9854)` was used with forked Mersenne-Twister workers. That is not reproducible: each child inherits the same state unless streams are used. The seed 9854 is kept as the master seed.

**B. Checkpoint bug.**

The checkpoint condition `n_done %% 50 == 0` stays true on every polling pass until the next task finishes. The partial results were therefore re-saved roughly every second, with cost growing with n.

### `Heterogeneous_effect_sim.R` → `sims/22_sim_heterogeneity.R`

**A. Effect sizes swapped.**

The pairs are (smaller, larger), but `eff_large <- pair[1]` and `eff_small <- pair[2]`. The evaluated group (rows `1:p_small`, column "Small_effect") therefore had the **larger** effect. That is the opposite of the design in Chapter 3. **Any version-2 heterogeneity result should be discarded.**

**B. No error handling around `MIMOSA2()`.**

One failure lost the task, and there was no time limit per fit (only the outer task kill).

**C. Size of the version-2 grid.**

The grid had 3 ρ × 3 splits × 9 pairs × 5 N × 50 repetitions = 20,250 tasks, each with two fits. That is the biggest computational cost, which is why the standard profile uses a subset.

### `Count_Imbalanced_Sim.R` → `sims/23_sim_imbalance.R`

**C. Scenario construction.**

The hand-typed list of 120 `rng` vectors was checked and every vector depletes the intended assays. It is now generated from (pattern, f) with explicit columns. Generating it from parameters removes the chance of a typing error and also fixes the plot script (below).

**B. RNG.** As for the prior study.

**A. Saved results did not match the script.**

`sim_progress.log` shows 450 tasks (15 scenarios × 30 repetitions), so `Count_Imbalance.Rdata` came from an older 15-scenario version.

### `Baseline_simulation_plots.R`, `Simulation_2.0.R` → `analysis/31_plots_baseline.R`

**A. AUC estimand undefined.**

AUC was computed after pooling subjects **across datasets**, and in the tables across P as well. That mixes posterior probabilities from different fitted models. Version 3 computes AUC per dataset, then the mean and MCSE.

**A. Invalid ΔAUC "95% CIs".**

The Hanley–McNeil SE used n = all pooled subjects (8,550 or 42,750), as if they were independent. It also fixed the positive fraction at 50%, which is wrong for 10% and 90% responders, and assumed a correlation r = 0.5. The resulting CIs were far too narrow. Version 3 uses paired per-dataset differences with MCSE.

**A. FDR plot used the wrong rule.**

`selected <- (1 - MIMOSA2_prob) <= a` thresholds each subject's *local* posterior. It does not use the cumulative Bayesian FDR (q-value) that `getResponse()` uses, and it pooled subjects across datasets and P. This makes FDR control look better than it is. Version 3 uses q < α per dataset, then the mean FDP with MCSE.

**B. Averaging over DGMs.**

`AUROC_matrix` averaged over P (and over ρ in the first table), and the ROC curves pooled over P.

**C. Other errors.**

- The AUROC computation was commented out, but `AUROC` was still used.
- The first AUROC table's method columns were swapped.
- The x-axis label "1-Sensitivity" should be 1 − specificity.
- The y-axis label had the typo "Observed FDF".
- `Simulation_2.0.R` calls `print(ROC_plot)` on an undefined object and `View()`, and has two `setwd()` calls.

### `Prior_simulation_plots.R` → `analysis/32_plots_prior.R`

**A. "Extremely Low" dropped.**

The factor levels omitted "Extremely Low" (3,000 cells), so those rows became NA in every figure and table.

**A/B.** The same AUC, Hanley–McNeil, FDR-rule and pooling issues as the baseline plots. The table marked "# INCORRECT" had its method columns swapped.

### `Count_Imbalanced_Plot.R` → `analysis/34_plots_imbalance.R`

**A. Every scenario mislabelled or dropped.**

The scenario names end in `_Imb_s1`, `_Bal_1`, and so on, but the script tested `grepl("_Imb$")` and `grepl("_Bal$")`. In the ROC figure every scenario was labelled "Balanced". In the AUC odds-ratio figure every non-reference scenario got `Balance = NA` and was filtered out.

**B. DeLong SE on pooled subjects.**

The DeLong SE was computed on subjects pooled across datasets. It treats subjects fitted by the same model as independent, which understates uncertainty.

**C.** `pROC::roc(direction = "<")` was set correctly here.

### `Variance_Solver.R`, `Bimodal_prior.R`, `Wonky_prior.R` → `R/dgm.R` (calibration), `analysis/32_plots_prior.R` (figures)

**B. Bimodal figure did not match the simulated distribution.**

`Bimodal_prior.R` plots a *different* bimodal distribution from the one simulated: a 50:50 mixture with mean 6×10⁻⁴, variance 2×10⁻⁸ and separation ratio 0.8. The thesis figure must use the simulated sampler. `prior_density` does.

**B. Solver output not used.**

`match_simplex_phi()` was never used, and its hard-coded comments disagree with the values in `Prior_simulations.R`. The solvers are now inside `R/dgm.R` and run automatically.

**C. Missing package.** `%>%` was used without loading magrittr or dplyr.

### Logs (`sim_progress.log`, `task_start.log`, `task_end.log`, `simulation_times.csv`)

**Run log.** `task_start.log`/`task_end.log` have 12,000 tasks = 120 scenarios × 100 repetitions, which is exactly the grid of `Count_Imbalanced_Sim.R` (P = 100). It used 20 workers, took 6 h 46 min (about 40 s per task per worker) and had no kills. `sim_progress.log` (450 tasks) is from an older 15-scenario version. Keep a note of which script and which log produced each `.Rdata`.

**Fit times.** `simulation_times.csv` has 1,919 fits of P = 50, with a median of 15.8 s and a maximum of 23.5 s. A 600 s limit therefore only catches genuinely stuck fits.

These timings were used to size the profiles. `11_smoke_test.R` re-measures them on your cluster.

## 3. Things that cannot be fixed by code (state them as limitations)

1. **Convergence is not observable.** MIMOSA2 does not return its iteration count or a convergence flag. Study 5 measures the practical consequence of stopping at 30 iterations instead.
2. **The meaning of the effect size.** δ is the nominal shift in the mean of p(s,1). The realised Δ among responders differs, because of truncation by rejection sampling and because profile R3's Δ does not depend on δ. Report the realised Δ (`<study>_realised_dgm.csv`).
3. **Base cases differ across studies.** φ is 2000, 10,000 or 5000; P is 10–100, 50, 100 or 20 + 80; ρ is 0.10–0.90, 0.50 or 0.25. The four studies therefore cannot be compared directly with each other. You kept the designs. Justify each choice in one sentence, or list it as a limitation (Morris Section 3.2).
4. **Fixed cell counts.** min = max in every range, so there is no subject-to-subject variation in N. That is unrealistic, but it keeps the cell count a clean experimental factor.
5. **Different error rates at the same nominal α.** The DiD unadjusted rule controls per-subject type I error, and MIMOSA2 controls the Bayesian FDR. The like-for-like comparator is DiD-BH; the threshold-free comparison is the AUC.

## 4. Notes for the write-up

### Chapter 3 (methods)

- **Section order.** Follow `STUDY_DESCRIPTIONS.md`: common DGM, then methods, targets, performance measures, nsim and seeds, then one subsection per study.
- **Sentence on random numbers.** "Each study used a single master seed with R's L'Ecuyer-CMRG generator; each simulated dataset was generated from its own random-number stream (`parallel::nextRNGStream`), whose state was stored so that any dataset can be regenerated exactly (Morris et al., 2019, Section 4.1)."
- **Sentence on nsim.** "We used n_sim = 100 repetitions per scenario. For a performance measure that is a mean of per-dataset proportions, the Monte Carlo SE is at most 0.5/√100 = 0.05 and was at most [max MCSE from `precision_report.csv`] in practice." Use 200 if you ran the extended profile.
- **Failures.** "A fit was recorded as failed if MIMOSA2 returned an error or did not finish within 600 seconds; failed fits are excluded from other performance measures and failure rates are reported per scenario."
- **Calls.** "MIMOSA2 calls used `getResponse()`, i.e. q < α with the Bayesian FDR of Newton et al. (2004); DiD calls used the one-sided Wald p-value, unadjusted or Benjamini–Hochberg adjusted."
- **AUC definition.** "AUC was computed within each dataset (Mann–Whitney form) and averaged across repetitions; ROC curves in figures are vertically averaged across datasets (Fawcett, 2006)."
- **Update `chapter3_methodology.tex`.** It was drafted for version 2. Change:
  - the replicate counts in the scenarios table: 30 → 100/200, and imbalance 15 scenarios → 57/120;
  - the timeouts paragraph: now 600 s everywhere;
  - the AUC paragraph: now per dataset;
  - "Uncertainty in performance estimates": now MCSE;
  - add the DiD-BH rule, the null scenarios and Study 5;
  - remove the TODOs that this version resolves.

### Chapter 4 (results)

- Lead with the failure rates (Morris Section 5.1). Then put the methods side by side, with MCSE in brackets or as the maximum MCSE in the caption (Morris Section 6.2).
- Round to what the MCSE supports: 2 decimals when the MCSE ≈ 0.01–0.02, 3 decimals for AUC when the MCSE ≈ 0.003.
- The null scenarios are a natural first result: does each method do what it promises when there is no signal?
- The prior study's key result is the difference from Beta, not the raw AUCs.
- The heterogeneity study's key result is the paired pooled − separate ΔAUC.
- The imbalance study's key result is imbalanced − matched balanced.
