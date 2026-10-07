# MIMOSA2 simulation study: version 3.1

Isabella Lethbridge & Tayyeb Abrahams. Prepared 6 October 2026.

Version 3.1 is a small update to version 3 (`mimosa2_sim_v3`). It fixes the problems found after the first full run (4–6 October). It also adds a way to correct the existing results without refitting MIMOSA2, and safeguards so that finished results cannot be overwritten again.

> **Not run yet.** R was not available where these changes were written. Every changed file was read line by line and checked for balanced brackets, and the key ideas were tested separately. Run the smoke test first (step 3 below). It runs the changed code end to end on a tiny grid.

The detailed evidence for every point (numbers, figures, log lines) is in `FIXES_v3.1.md`.

---

## 1. What happened in the first run

| Study | Datasets | Status |
|---|---|---|
| prior | 10,000 | finished and combined, 0 failures |
| imbalance | 5,700 | finished and combined, 0 failures |
| heterogeneity | 3,000 | finished and combined, 0 failures |
| maxit_check | 600 | finished and combined, 0 failures |
| baseline | 46,800 | finished and combined on 5 Oct, but the results file has since been lost |

**The baseline did not crash.**

- **5 Oct, 06:24:** the baseline finished all 46,800 tasks with 0 failed tasks. It was combined into `baseline_results_standard.rds`, and the performance tables were made from it (`slurm_mimosa2_all_1401889.out`). Within those tasks, 76 MIMOSA2 fits (0.16%) returned errors (section 4).
- **Afterwards:** the task files and the results file disappeared.
- **6 Oct:** the baseline was resubmitted. It found no task files, so it silently restarted from task 1. That job was then cancelled at 12:45 (`slurm_mimosa2_sim_1414662.out`).

Because the results file no longer exists, **the baseline has to be re-run** (about 10 h on 20 cores). The other four studies do not need to be re-run.

---

## 2. Changes made to the version 3 scripts

Each change is marked `[CHANGE v3.1]` or `[v3.1]` in the code.

### 2.1 DiD was overconfident when a cell count was zero (this changes results)

**Files:** `R/methods.R` (`DiD_wald`), `R/config.R` (new setting `DID_CC`)

**What DiD does.** For each subject, DiD estimates whether the stimulated-minus-unstimulated gap grew from baseline to after vaccination. It then divides that change by its standard error. A large z means "probably a responder".

**Where it went wrong.**

- The standard error used the observed proportion in each assay, √Σ p̂(1−p̂)/N.
- When an assay had zero positive cells, p̂ = 0, so that assay counted as having no uncertainty at all. But seeing 0 positives in a few thousand cells does not mean the true proportion is exactly 0. There are just too few cells to tell.
- The standard error was therefore too small, and z was far too large.

**Effect.** In the imbalance study, with one assay depleted to 2.5%:

- Depleting **s0**: 26% of non-responders were called responders at α = 0.05, and the DiD-BH FDR was 0.53.
- Depleting **u1**: 20% false calls, and an FDR of 0.38.
- Depleting **s1 or u0**: the opposite problem. DiD found almost no responders (DiD-BH sensitivity ≈ 0.005).
- It also inflated DiD's false calls at Very Low and Extremely Low cell counts in the prior study.

**How it slipped in.**

- Version 2 used a per-subject binomial GLM. With a zero count, that GLM usually failed or did not converge, and the subject was given the tie score 0.05 ("not a responder"). So version 2 never produced these false positives.
- Version 3 replaced the GLM with a closed-form formula that gives the same answer when there are no zero counts.
- The smoke test compared the two only on subjects without zero counts, so the difference was never seen.

**Fix.** The variance now uses p̃ = (x + 0.5)/(N + 1) instead of p̂. The DiD estimate itself is unchanged.

- This is a standard continuity correction, and the same formula the legacy `DiD_bvn_mixture` already used.
- With normal cell counts it changes the variance by about 3%.
- In a stand-alone check, false calls under s0 depletion dropped from 26–33% to 0.3–7.6%. Balanced scenarios were essentially unchanged (0.050 became 0.048).
- `DiD_wald(cc = 0)` reproduces version 3.0 exactly.

### 2.2 Applying the DiD fix without refitting MIMOSA2 (new script)

**File:** `sims/40_recompute_did.R`

DiD needs only the cell counts, and every dataset can be rebuilt exactly from its stored random-number stream. For each study the script:

1. **Rebuilds every dataset** from its seed.
2. **Checks reproducibility.** The uncorrected DiD must equal the value stored in the results file for every subject. If the maximum difference is not 0, the script stops without writing anything.
3. **Writes the corrected DiD** into the results file. It replaces the DiD rows of `estimates` and the DiD columns of `subjects`, and leaves MIMOSA2 untouched.
4. **Keeps the old results.** It makes a backup, `<study>_results_standard_before_didfix.rds`, and keeps the old DiD rows in `res$estimates_did_v30`.

It skips any study that was already corrected, or that was run with the v3.1 code (such as the new baseline).

### 2.3 ROC curves had a flat line and then a jump at the end (figures only)

**File:** `analysis/analysis_functions.R` (`roc_at_grid`)

**Where it went wrong.**

- About 10% of MIMOSA2 posterior probabilities are exactly 0 (numerical underflow). For NR2 non-responders it is 85–93%.
- Tied subjects have no correct order among themselves, so the fair way to draw them is a straight diagonal line. The code drew them as a flat line followed by a vertical jump at FPR = 1.
- The ROC curves therefore looked worse than the method performed. The AUC values were already correct, because they count ties as ½.

**Fix.** Tied blocks are interpolated linearly. Legacy plotting (`plotROC`, `pROC`) also draws ties as a diagonal. The ROC panels now also show the full 0–1 axes; the y-axis used to start at about 0.25.

**Files:** `analysis/31_plots_baseline.R`, `analysis/32_plots_prior.R`

### 2.4 The prior-distribution figure showed the wrong distribution

**File:** `analysis/32_plots_prior.R`

**Where it went wrong.** Every panel of `prior_density` was drawn as a stimulated draw (`cond = "s1"`). The "unstimulated background" panel therefore showed a two-peaked Bimodal Beta, but the simulation always draws unstimulated proportions from the Beta.

**Fix.** Each panel uses the condition it is labelled with. A third panel was added for the stimulated baseline mean (3.5×10⁻⁴), where the bimodal family is used.

### 2.5 The "Δ AUC" axis labels showed as "AUC"

**Files:** `analysis/31_plots_baseline.R`, `33_plots_heterogeneity.R`, `34_plots_imbalance.R`

The PDF font on the cluster has no Greek Δ, so it was silently dropped. The labels now read, for example, "Difference in AUC (MIMOSA2 − DiD)".

### 2.6 Overlapping x-axis in the heterogeneity figures

**File:** `analysis/33_plots_heterogeneity.R`

On a log scale, the effect sizes 6.25×10⁻⁴, 8×10⁻⁴ and 1×10⁻³ sat on top of each other. The axis is now categorical, with evenly spaced values.

### 2.7 Detection by response profile (new table)

**Files:** `analysis/analysis_functions.R` (`profile_table`), `analysis/30_performance_tables.R`

The new file `_tables/<study>_by_profile.csv` gives, for each scenario and each of the eight profiles (R1–R4, NR1–NR3, NSR), the share of subjects called at α = 0.01 and the median true Δ. It shows which profiles limit sensitivity (R3) and which drive false calls (section 4).

### 2.8 Safeguards: finished results can no longer be overwritten

**Files:** `R/runner.R`, `sims/29_combine.R`

- **`run_study()` stops instead of silently restarting.** If `<study>_results_<profile>.rds` exists but none of the task files do, it stops and explains why. Before, it restarted all tasks, which is what happened to the baseline on 6 Oct. To re-run anyway, set `MIMOSA2_FORCE_RERUN=1`.
- **`combine_study()` refuses to run with 0 task files.** Before, `29_combine.R baseline` would have overwritten the results with an empty file.
- **An existing results file is kept.** It is renamed to `*_previous.rds` instead of being overwritten.
- **One failing study no longer blocks the others.** `29_combine.R` now continues to the next study if one fails to combine.
- **The run records which DiD was used.** `run_info` stores `did_cc`, which lets `40_recompute_did.R` skip runs that already have the corrected DiD.

### 2.9 Smoke test

**File:** `sims/11_smoke_test.R`

The test that compares DiD with the version-2 GLM now uses the uncorrected statistic (`cc = 0`). That is the statistic the GLM computes, so the test still checks the right thing.

### 2.10 A check that had never been run

`analysis/35_maxit_check.R` was not run on the full results. The maxit = 30 justification and the baseline reproducibility check were therefore missing. `slurm/fix_v31.sh` now runs it. The script itself is unchanged.

### 2.11 New job script

**File:** `slurm/fix_v31.sh`

It runs `40_recompute_did.R`, then `analysis/30` to `35`, and stops at the first error. It needs all five results files.

---

## 3. Consistency with the legacy (version 2) code

The changes were checked against `Non_beta_simulations.R`, `Simulation.R`, `Prior_simulations.R`, `Heterogeneous_effect_sim.R`, `Count_Imbalanced_Sim.R` and the legacy plotting scripts.

| v3.1 change | Legacy code | Verdict |
|---|---|---|
| DiD continuity correction | The legacy GLM did not call zero-count subjects. The legacy `DiD_bvn_mixture` uses the same (x + 0.5)/(N + 1) formula. | Restores the legacy behaviour, using a legacy formula |
| ROC tie handling | `plotROC` and `pROC` draw ties as a diagonal | Matches legacy |
| Density figure | Legacy unstimulated draws were always Beta | Matches legacy |
| Labels, axes, safeguards, new scripts and tables | No legacy equivalent | No change to the concepts |

**The data-generating mechanism is identical to legacy:** the eight profiles, the condition means (μ₀ = 1×10⁻⁴, δ₀ = 2.5×10⁻⁴, γ = 0), the rejection constraints, the binomial counts and the definition of R3.

**One earlier version-3 change does depart from legacy: which draws are bimodal** (version 3, CHANGE 4). It affects only the Bimodal Beta scenarios of the prior study. Decide which rule to use, and describe it in Chapter 3.

| Profile | Legacy: bimodal draws | Version 3: bimodal draws |
|---|---|---|
| R1, R4 | p(s,1) only; p(s,0) Beta | p(s,1) and p(s,0) |
| R2, R3 | p(s,1) (= p(s,0) in R3) | same |
| NR1 | shared p(s,1) = p(u,1) | none (shared draws treated as unstimulated) |
| NR2 | shared p(s,1) = p(u,1); p(s,0) Beta | p(s,0) only |
| NR3 | all four (shared) | none |
| NSR | p(s,0) = p(s,1) | same |

---

## 4. Results that look odd but are not bugs (explain them in the write-up)

1. **Sensitivity and AUC barely rise with δ, and sensitivity stays below 1 even at δ = 10⁻³ with high cell counts.** This comes from the original design, not from version 3:
   - **Rejection sampling changes the realised effect.** The mean realised Δ among responders is 4.2× δ at δ = 6.25×10⁻⁵ but 0.79× δ at δ = 10⁻³. The nominal 16-fold range of δ is only a 3-fold range of realised Δ (2.6×10⁻⁴ to 7.9×10⁻⁴; `baseline_realised_dgm.csv`).
   - **R3's effect does not depend on δ.** In profile R3, p(s,1) = p(s,0), so Δ = p(u,0) − p(u,1) (median ≈ 1×10⁻⁴). In the prior study (Beta, High cells, δ = 10⁻³), MIMOSA2 calls 96–100% of R1, R2 and R4 subjects but only 68% of R3. DiD's median z for R3 is 0.9.
   - Describe δ as the "nominal effect" and report the realised Δ and the by-profile table.
2. **MIMOSA2's Bayesian FDR is not controlled when responders are rare.**
   - Null scenarios at α = 0.05: P(any false call) reaches 0.19 (P = 100, Low cells), against at most 0.06 for DiD-BH.
   - At ρ = 10%, the mean FDP reaches 0.07 at P = 100.
   - At α = 0.01 the FDR is controlled everywhere. This is a real finding.
3. **76 MIMOSA2 fit errors in the baseline (0.16%).** All are at P = 10 or 20 with at most 3 responders (ρ ≤ 0.25 or ρ = 0). To see the messages: `table(res$fits$Msg[res$fits$Status == "error"])`.
4. **Posterior underflow costs MIMOSA2 some AUC.** Tied zeros count ½, which partly explains MIMOSA2 − DiD < 0 at small δ in `base_dAUC`.
5. **Some comparisons are not paired.** Prior family vs Beta, and imbalanced vs balanced, use independent datasets. Differences in AUC smaller than about ±0.025 cannot be detected in the prior study. State this as a limitation.
6. **Rounding at P = 10.** The realised ρ is 0.30 for "25%" and 0.80 for "75%". The realised value is stored for every dataset.
7. **Heterogeneity.** Pooling hurts the weak responders most when they are the majority (80 focal : 20 other). This is not a labelling error: the focal subjects are the first `P_focal` subjects.

---

## 5. How to proceed

Run everything from inside the `mimosa2_sim_v3` folder on the cluster.

**Step 0. Make sure the baseline results file is really gone:**

```bash
find /scratch/abrmoe030 -name "baseline_results*" 2>/dev/null
```

If this finds it, copy it to `_simulations/baseline/`, skip step 4, and submit only `sbatch slurm/fix_v31.sh` in step 5.

**Step 1. Copy in the v3.1 files.** Copy every file from `_fixes_v3.1/` into the same subfolder of `mimosa2_sim_v3/`, for example `_fixes_v3.1/R/methods.R` to `mimosa2_sim_v3/R/methods.R`. Do this before anything is run. Keep the existing copies of the files that were not changed: `dgm.R`, `scenarios.R`, sims 10 and 20–24, `35_maxit_check.R` and the other slurm scripts.

**Step 2. Move the leftovers of the old baseline run aside:**

```bash
mkdir -p _old_baseline_run
mv _simulations/baseline/progress_chunk00.log _simulations/baseline/run_info_chunk00.rds _old_baseline_run/
```

**Step 3. Run the smoke test (about 30 min):**

```bash
sbatch slurm/smoke_test.sh
```

The end of `slurm_mimosa2_smoke_<jobid>.out` must show **"ALL DGM CHECKS PASSED"** and **"All checks passed"**. Then delete the `_smoke/` folder.

**Step 4 (baseline re-run, about 10 h) and step 5 (fix job, 1–2 h).** The fix job is set to start automatically when the baseline finishes successfully:

```bash
JOB=$(sbatch --parsable slurm/run_study.sh baseline)
sbatch --dependency=afterok:$JOB slurm/fix_v31.sh
```

The baseline uses exactly the same 46,800 datasets as on 5 October (same seeds), now with the corrected DiD. The fix job corrects DiD in the other four studies (no MIMOSA2 refits) and then remakes every table and figure, including the maxit check.

**Step 6. Check the logs.**

- Baseline `.out`: it should end with `finished: 46800 ok` and `combining 46800 of 46800 planned tasks`.
- Fix job `.out`:
  - The baseline should be reported as "DiD already corrected; skipped".
  - For prior, heterogeneity, imbalance and maxit_check, `max |DiD_est regenerated - stored|` must be **0**.
  - The "non-responders called by unadjusted DiD" line should drop clearly for the imbalance study.
  - `35_maxit_check.R` should report "identical AUC and TP: TRUE".
  - The job should end with "done".
- Figures in `_fig/`:
  - No ROC curve should jump at the right edge.
  - `prior_density` should have three panels.
  - The ΔAUC figures should have axes reading "Difference in AUC".

**`.Rdata` copies.** At the end, the fix job also writes `_simulations/<study>/<study>_results_standard.Rdata` next to each `.rds` file. `load()` it to get `<study>_results` (the full list) and `<study>_design`, `_datasets`, `_fits`, `_estimates` and `_subjects` as separate data frames. The pipeline itself keeps reading the `.rds` files. To make the copies on their own: `SIM_PROFILE=standard Rscript sims/41_export_rdata.R`.

**Step 7. Back up the results.** Copy all five `_simulations/<study>/<study>_results_standard.rds` files somewhere safe. Delete a `tasks/` folder only after its results file exists and has been backed up.

**If something fails:** send the `.out` file. If the baseline job hits its time limit, submit the same command again. It continues where it stopped, provided the `tasks/` folder is left in place.

---

## 6. Files in `_fixes_v3.1`

| File | Status | Change |
|---|---|---|
| `R/config.R` | changed | `DID_CC = 0.5` |
| `R/methods.R` | changed | continuity-corrected DiD variance |
| `R/runner.R` | changed | safeguards; records `did_cc` |
| `sims/11_smoke_test.R` | changed | compares the uncorrected DiD with the GLM |
| `sims/29_combine.R` | changed | one failing study no longer stops the others |
| `sims/40_recompute_did.R` | **new** | applies the DiD fix to existing results |
| `analysis/analysis_functions.R` | changed | ROC ties; `profile_table()` |
| `analysis/30_performance_tables.R` | changed | writes `<study>_by_profile.csv` |
| `analysis/31_plots_baseline.R` | changed | ROC axes; ΔAUC label |
| `analysis/32_plots_prior.R` | changed | density figure; ROC axes |
| `analysis/33_plots_heterogeneity.R` | changed | categorical x-axis; ΔAUC label |
| `analysis/34_plots_imbalance.R` | changed | ΔAUC label |
| `sims/41_export_rdata.R` | **new** | saves a `.Rdata` copy of every results file (the `.rds` files are kept) |
| `slurm/fix_v31.sh` | **new** | runs 40, analysis 30–35, then 41 |
| `FIXES_v3.1.md` | **new** | detailed review with evidence |
| `README_v3.1.md` | **new** | this file |

---

## 7. For Chapter 3 (methods)

- **DiD:** a one-sided Wald test of the difference-in-differences. Its variance uses the continuity-corrected proportions (x + 0.5)/(N + 1). Calls are unadjusted, or Benjamini–Hochberg adjusted.
- **Effect size:** δ is the nominal post-vaccination effect. Report the realised Δ among responders and explain the R3 profile.
- **Bimodal family:** state which draws are bimodal (section 3).
- **Limitations:**
  - Comparisons across scenarios are not paired (section 4, item 5).
  - MIMOSA2 does not return its iteration count; Study 5 addresses this.
