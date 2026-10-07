# MIMOSA2 simulation v3.1: review of the 4–6 October run and fixes

Prepared 6 October 2026. These are the files changed relative to `mimosa2_sim_v3`. Copy them over the originals, keeping the same paths. Every change is marked `[CHANGE v3.1]` or `[v3.1]` in the code.

## 1. The baseline did not crash

| Evidence | What it shows |
|---|---|
| `slurm_mimosa2_all_1401889.out`, lines 851–853 | `[baseline] finished: 46800 ok, 0 failed`, then `combining 46800 of 46800 planned tasks`, then `end sims/20_sim_baseline.R` (5 Oct, 06:24 SAST) |
| same file, lines 859–864 | `30_performance_tables.R` loaded `baseline_results_standard.rds` and wrote all baseline tables |
| `_tables/baseline_dAUC_table.csv`, `baseline_fdr_calibration.csv` | `31_plots_baseline.R` was also run on the baseline results later |
| `slurm_mimosa2_sim_1414662.out` | the re-submission on 6 Oct found 0 task files, restarted from task 1, and was **cancelled** (not crashed) at 12:45 |

So `baseline_results_standard.rds` was written on 5 Oct. It is not in the copy of the folder you shared, and neither are any `tasks/` folders.

**Next step on the cluster:** run `ls -la _simulations/baseline/`.

- If the results file is there, back it up and use it.
- If it is gone, the baseline has to be run again (about 10 h). Submit `sbatch slurm/run_study.sh baseline` with the v3.1 `runner.R`. It cannot silently overwrite anything any more.

The second run did change `design_full.rds` and `run_info_chunk00.rds`, and it removed the finished marker. The new design is identical (same seed), so this is harmless. The `run_info` copy inside the results file still describes the 5 Oct run.

## 2. Problems found, and what was changed

### A. Must fix: these change results

1. **The DiD standard error is zero for a zero count** (`R/methods.R`, `DiD_wald`). The SE used p̂(1−p̂)/N, so an assay with no positive cells contributed no variance.
   - Imbalance study, 2.5% depletion: depleting s0 gave 26% false calls among non-responders at α = 0.05 (DiD-BH FDR 0.53). Depleting u1 gave 20% (FDR 0.38). Depleting s1 or u0 left DiD with almost no power.
   - This also inflates DiD's unadjusted FDR in the prior study at Very Low and Extremely Low cell counts, and it affects DiD's AUC everywhere counts are small.
   - **Fix:** a continuity-corrected variance, p̃ = (x + 0.5)/(N + 1). The estimate itself is unchanged. `DID_CC = 0.5` is set in `config.R`, and `DiD_wald(cc = 0)` reproduces v3.0.
   - **Check:** a stand-alone re-implementation reproduced the problem (s0 depletion: 26–33% false calls). With the correction the rate is 0.3–7.6%, and it barely changes anything with balanced counts (0.050 becomes 0.048).
   - **Applying it without refitting MIMOSA2:** `sims/40_recompute_did.R` regenerates every dataset from its stored stream and recomputes DiD. It first checks that the uncorrected DiD equals the stored value for every subject, which also completes the reproducibility check. The original results files are backed up.

2. **`35_maxit_check.R` was never run on the full results** (there is no `maxit_check_summary.csv` in `_tables`). The maxit = 30 justification and the baseline reproducibility check therefore exist only from the smoke test. `slurm/fix_v31.sh` runs it.

### B. Figures that looked wrong (analysis code)

3. **ROC curves were flat and then jumped at FPR = 1** (`analysis_functions.R`, `roc_at_grid`). About 10% of MIMOSA2 posteriors are exactly 0, a numerical underflow. NR2 non-responders: 85–93%.
   - A block of tied scores is a diagonal segment of the ROC curve, but it was drawn as a step: a long flat line, then a jump at FPR = 1. The curves therefore looked worse than the AUC says.
   - **Fix:** linear interpolation inside tie blocks, the same convention as the Mann–Whitney AUC. The ROC panels now also use 0–1 axes; before, the y-axis started at about 0.25.

4. **`prior_density` showed the wrong distribution.**
   - Both panels were drawn with `cond = "s1"`, so the "unstimulated background" panel showed a bimodal curve for Bimodal Beta. The simulation draws unstimulated proportions from the Beta.
   - **Fix:** each panel uses its own condition, and a panel for the stimulated baseline (3.5×10⁻⁴) was added.

5. **"Δ AUC" axis labels rendered as "AUC"** in `base_dAUC`, `het_dAUC_pooled_vs_separate` and `imbalance_vs_balanced`. The Greek letter is missing from the cluster's PDF fonts. The labels are now plain text.

6. **Heterogeneity figures:** on a log axis, 6.25×10⁻⁴, 8×10⁻⁴ and 1×10⁻³ sat on top of each other, so the axis is now categorical.

### C. Safety (what went wrong on 6 Oct)

7. `run_study()` now **stops** if `<study>_results_<profile>.rds` exists but none of the chunk's task files do. Before, it silently restarted all tasks. Set `MIMOSA2_FORCE_RERUN=1` to override.
8. `combine_study()` refuses to run with 0 task files. Before, `29_combine.R baseline` would have overwritten the results with an empty file. It also renames an existing results file to `*_previous.rds` instead of overwriting it.

### D. Not code bugs: explain these in the write-up

9. **The TPR and AUC curves against δ are flat, and TPR stays below 1 even at δ = 10⁻³ with high cell counts.** There are two reasons:
   - **The realised effect differs from the nominal δ** (`baseline_realised_dgm.csv`). Rejection sampling makes the mean realised Δ among responders 4.2× δ at 6.25×10⁻⁵, but 0.79× δ at 10⁻³. The nominal 16-fold range of δ is only a 3-fold range of realised Δ (2.6×10⁻⁴ to 7.9×10⁻⁴).
   - **Profile R3's Δ does not depend on δ** (median Δ ≈ 1×10⁻⁴). In the prior study (Beta, High, δ = 10⁻³), MIMOSA2 calls R1/R2/R4 96–100% of the time but R3 only 68%. DiD's median z for R3 is 0.9.
   - The new table `<study>_by_profile.csv` (from `30_performance_tables.R`) shows this for every scenario.
10. **MIMOSA2's Bayesian FDR is not controlled when responders are rare.**
    - Null scenarios at α = 0.05: P(any false call) reaches 0.19 at P = 100, Low cells (DiD-BH: at most 0.06).
    - ρ = 10%: the mean FDP reaches 0.07 at P = 100.
    - In `fdr_plot_clean`, MIMOSA2 is above the diagonal at 10% responders and below it at 50% and 90%.
    - At α = 0.01, the FDR is controlled everywhere. This is a real finding; report it.
11. **The 76 baseline MIMOSA2 errors (0.16%)** are all at P = 10 or 20 with ρ ≤ 0.25 or ρ = 0, so with at most 3 responders. To see the messages: `table(res$fits$Msg[res$fits$Status == "error"])`.
12. **The posterior underflow (item 3) also costs MIMOSA2 some AUC.** Tied zeros count 1/2, which partly explains MIMOSA2 − DiD < 0 at small δ in `base_dAUC`. Mention it when interpreting ΔAUC.
13. **Comparisons across scenarios are not paired.** Prior family vs Beta, and imbalanced vs balanced, use independent datasets, so their MCSE is about √2 × the per-scenario MCSE. In the prior study the bars are about ±0.025 AUC, so differences smaller than that cannot be detected. State this as a limitation.
14. **Realised ρ at P = 10 is 0.30 for "25%" and 0.80 for "75%"** (rounding of ρP). This is already stored in the table; mention it.
15. **Heterogeneity:** pooling hurts the weak responders most when they are the majority (80 : 20). That's the opposite of what one might expect, so it needs explaining. It is not a bug: the focal subjects are the first `P_focal` subjects, as the labels say.

## 3. What to run (one job, about 1–2 h, no MIMOSA2 refits)

```bash
# 0. copy the v3.1 files over mimosa2_sim_v3 (same paths)
# 1. make sure all five results files exist (restore the baseline one if needed)
ls -la _simulations/*/*_results_standard.rds
# 2. fix DiD + redo every table and figure
sbatch slurm/fix_v31.sh
```

The job log prints, for each study, the reproducibility check (`max |DiD_est regenerated - stored|` must be 0) and the DiD false-call rate before and after the fix.

## 4. Files in this folder

| File | Change |
|---|---|
| `R/config.R` | `DID_CC = 0.5` |
| `R/methods.R` | `DiD_wald(..., cc = DID_CC)` |
| `R/runner.R` | safety guards (items 7–8) |
| `sims/11_smoke_test.R` | compares the uncorrected DiD with the GLM (`cc = 0`) |
| `sims/29_combine.R` | one failing study no longer stops the others |
| `sims/40_recompute_did.R` | **new**: applies the DiD fix to existing results |
| `analysis/analysis_functions.R` | ROC tie fix; new `profile_table()` |
| `analysis/30_performance_tables.R` | writes `<study>_by_profile.csv` |
| `analysis/31`–`34` | axis and label fixes; corrected density figure |
| `slurm/fix_v31.sh` | **new**: runs steps 40, 30–35 |

R was not available where this was written. The files were checked for balanced brackets and read through line by line, but they have not been run. `fix_v31.sh` stops at the first error.

## 5. Consistency with the legacy (version-2) code

The legacy files compared were `Non_beta_simulations.R` (DGM, `DiD_GLM`, `DiD_bvn_mixture`), `Simulation.R`, `Prior_simulations.R`, `Heterogeneous_effect_sim.R`, `Count_Imbalanced_Sim.R` and the plotting scripts.

| v3.1 change | Legacy concept | Verdict |
|---|---|---|
| DiD variance with p̃ = (x + 0.5)/(N + 1) | `DiD_GLM` was the identity-link binomial GLM. When a count was 0, the fitted proportion sat on the boundary, the GLM usually failed or did not converge, and the subject got the tie score 0.05, i.e. "not a responder". So legacy never produced the false positives seen in v3.0. The v3.0 closed form matched the GLM only for subjects with no zero counts, which is all the smoke test checked. The same correction, `(x + 0.5)/(N + 1)`, is already used for the variance in the legacy `DiD_bvn_mixture` | Restores the legacy behaviour at zero counts without the tie problem. Uses a correction from the legacy code. With no zero counts it differs from the GLM by about 3% in variance at N = 150,000, p = 1e-4 |
| ROC tie interpolation | Legacy ROC/AUC used `plotROC::geom_roc` (baseline, prior) and `pROC` (imbalance, heterogeneity). Both draw tied scores as a diagonal segment | Back in line with legacy |
| Density figure: unstimulated panel Beta | Legacy unstimulated draws (pu0, pu1) were always Beta under the bimodal family | In line with legacy |
| Axis labels, 0–1 ROC axes, categorical heterogeneity axis | Cosmetic | No conceptual change |
| Safety guards in `runner.R`; `40_recompute_did.R`; profile table | No legacy equivalent (legacy kept results in memory) | New infrastructure only. DGM, effect sizes, profiles and constraints are unchanged |

**DGM structure is unchanged from legacy.** The 8 profiles, condition means (μ0 = 1e-4, δ0 = 2.5e-4, γ = 0), rejection constraints, binomial counts and the R3 definition (p(s,1) = p(s,0) at μ(s,1), with p(u,0) > p(u,1)) are identical. The flat δ curves and the R3 cap are therefore features of the original design, not of v3.

**One v3.0 departure from legacy that v3.1 does not undo: which draws are bimodal** (v3.0 CHANGE 4, made after the audit). Decide which you want and describe it in Chapter 3.

| Profile | Legacy bimodal draws | v3 bimodal draws |
|---|---|---|
| R1, R4 | p(s,1) only (p(s,0) Beta) | p(s,1) and p(s,0) |
| R2, R3 | p(s,1) (and p(s,0) = p(s,1) in R3) | same |
| NR1 | shared p(s,1) = p(u,1) at μ(u,1) | none (shared draws treated as unstimulated) |
| NR2 | shared p(s,1) = p(u,1); p(s,0) Beta | p(s,0) only |
| NR3 | all four (shared) | none |
| NSR | p(s,0) = p(s,1) | same |

It affects only the Bimodal Beta scenarios of the prior study. Every other study uses the Beta.
