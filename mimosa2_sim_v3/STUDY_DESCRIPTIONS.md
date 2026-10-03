# What each simulation study does

This file describes the five simulation studies in version 3 of the pipeline. Each one is written in the ADEMP order of Morris, White and Crowther (2019, Section 3): **A**ims, **D**ata-generating mechanisms, **E**stimands/targets, **M**ethods, **P**erformance measures. Morris et al. recommend that order for the methods chapter (Section 6.1), so each block below can be turned into Chapter 3 text almost line by line.

Studies 1 to 4 are the four analyses you already had. Study 5 is new and small. It answers the supervisor's "justify maxit = 30" comment and also shows that the stored random-number states reproduce datasets exactly.

---

## Part 0. What all studies share

### 0.1 Data-generating mechanism (common core)

All studies use `simulate_MIMOSA2_alt_prior()` in `R/dgm.R`. This is your modified version of `MIMOSA2::simulate_MIMOSA2`. Version 3 fixes some bugs in it (see `CHANGES_AND_AUDIT.md`).

**Step 1: assign subjects to response profiles.** Each of the P subjects gets one of eight profiles:

| Profile | Status | How the four true proportions are tied | Constraints (enforced by redrawing) |
|---|---|---|---|
| R1 | responder | all four drawn separately | Δ > 0, p(s,1) > p(u,1) |
| R2 | responder | p(s,0) = p(u,0) | p(s,1) > p(u,1) |
| R3 | responder | p(s,1) = p(s,0), drawn at the mean μ(s,1) | Δ > 0, p(s,1) > p(u,1), p(u,0) > p(u,1) |
| R4 | responder | p(u,1) = p(u,0) | Δ > 0, p(s,1) > p(u,1), p(s,1) > p(s,0) |
| NR1 | non-responder | p(s,1) = p(u,1), p(s,0) = p(u,0) | none |
| NR2 | non-responder | p(s,1) = p(u,1) | p(s,0) ≥ p(u,0) |
| NR3 | non-responder | all four equal | none |
| NSR | non-responder | p(s,1) = p(s,0), drawn at the mean μ(s,0); p(u,1) = p(u,0) | none |

Here Δ = (p(s,1) − p(u,1)) − (p(s,0) − p(u,0)) is the true difference-in-differences. A subject is a true responder exactly when Δ > 0, and `10_dgm_checks.R` confirms this.

With nominal responder proportion ρ, the number of responders is fixed first as n_R = round(ρP), with halves rounded up. The responders are then split as evenly as possible over R1–R4, and the P − n_R non-responders over NR1–NSR. Any remainder goes to randomly chosen profiles.

**Step 2: draw the true proportions.** The four condition means are:

- μ(u,0) = μ₀ = 10⁻⁴ (background)
- μ(s,0) = μ₀ + δ₀, with δ₀ = 2.5×10⁻⁴ (antigen-specific response before vaccination)
- μ(u,1) = μ₀ + γ, with γ = 0 (no change in background)
- μ(s,1) = μ₀ + γ + δ, where δ is the nominal **effect size**

Each proportion is drawn from the study's generative family. The Beta family uses Beta(μφ, (1 − μ)φ). Inequality constraints are imposed by redrawing until they hold.

**Step 3: draw the counts.** Total cells N per assay are fixed by the scenario. In the imbalance study they differ by assay. Positive cells are drawn as n ~ Binomial(N, p).

**Realised effect size (report this).** The nominal δ is not the realised difference-in-differences among responders. There are two reasons:

1. Rejection sampling truncates the distributions, which inflates small effects.
2. In profile R3, p(s,1) = p(s,0), so its Δ = p(u,0) − p(u,1). This does not depend on δ at all.

A Python check of the version 3 DGM with Beta and φ = 2000 gave these ratios of mean realised Δ to nominal δ:

| δ | R1 | R2 | R3 | R4 |
|---|---|---|---|---|
| 1×10⁻³ | 0.98 | 1.05 | 0.16 | 0.97 |
| 2.5×10⁻⁴ | 1.8 | 1.5 | 0.6 | 1.9 |
| 6.25×10⁻⁵ | 5.5 | 3.5 | 2.3 | 5.9 |

`sims/10_dgm_checks.R` writes the exact values to `_tables/dgm_check_C_large_datasets.csv`, and every dataset stores its mean and median realised Δ. In the write-up, describe δ as the "nominal post-vaccination effect" and report the realised Δ (Morris Section 3.2, pitfall 4). This structure comes from MIMOSA2's own simulator. It is not a bug, but it explains why performance changes less across δ than the nominal values suggest.

### 0.2 Methods (every study)

1. **MIMOSA2** (R package `MIMOSA2`). The settings are `maxit = 30`, `verbose = FALSE` and all other arguments at their defaults. Each fit has a hard time limit of 600 seconds; a fit that runs longer is recorded as a "timeout" failure. Outputs:
   - **P(responder) for each subject:** the sum of the posterior probabilities of components 1–4.
   - **Bayesian FDR q-values:** the running mean of (1 − P(responder)) over subjects ranked by P(responder).
   - **Responder calls at level α:** `getResponse(fit, α)`, which is q < α.
2. **DiD comparator.** This is a per-subject Wald test of H₀: Δ ≤ 0 against H₁: Δ > 0. The estimate is the observed difference-in-differences, and its SE is √(Σ p̂(1 − p̂)/N) over the four assays. This is exactly the interaction term of the saturated identity-link binomial GLM from version 2; the smoke test checks the two agree wherever the GLM converged with a non-negative estimate and no zero counts. Two decision rules are applied:
   - **unadjusted:** one-sided p ≤ α, which controls the per-subject type I error;
   - **BH:** Benjamini–Hochberg adjusted p ≤ α, which controls the FDR. This rule is the like-for-like comparison with MIMOSA2's Bayesian FDR.

   The ranking score for ROC/AUC is the z statistic.

Both methods are applied to **the same simulated datasets** (Morris Section 4.3). Each comparison is therefore paired, which reduces the Monte Carlo error of method differences.

α ∈ {0.01, 0.05} throughout.

### 0.3 Targets and performance measures (every study)

The **target** is each subject's binary responder status (R = 1 if Δ > 0). In Morris's terms (Table 3) this is a classification/testing target, not an estimand.

These measures are computed **in each simulated dataset** and then averaged over the n_sim datasets of a scenario. Scenarios are never pooled (Morris Section 5.2).

| Measure | Per-dataset definition | Monte Carlo SE of the mean |
|---|---|---|
| Fit failures | error or timeout (MIMOSA2); task crash | √(p(1 − p)/n_planned) |
| Sensitivity (TPR) | TP/(TP + FN), when the dataset has ≥ 1 responder | SD/√n |
| Specificity (TNR), FPR | TN/(TN + FP), FPR = 1 − TNR | SD/√n |
| FDR | mean of FDP = FP/(TP + FP), with FDP = 0 if nothing is called | SD/√n |
| P(any false positive) | 1 if FP > 0 | SD/√n |
| AUC | Mann–Whitney AUC within the dataset | SD/√n |
| ΔAUC (MIMOSA2 − DiD) | paired difference in the same dataset | SD(differences)/√n |
| FDR calibration curve | FDP at α = 0.01, …, 0.20 using q < α (MIMOSA2) or BH (DiD) | SD/√n |
| ROC curve (figures) | TPR of each dataset at a fixed grid of FPR values, averaged ("vertical averaging", Fawcett 2006) | SD/√n at each grid point |

Failed fits are excluded from the other measures, and the failure rate is reported next to them (Morris Section 5.1). Tables report "estimate (MCSE)". Figures show bands of ±1.96 MCSE (Morris Section 6.2).

### 0.4 Number of repetitions, seeds and profiles

- **n_sim.** The "standard" profile uses 100 repetitions per scenario. The "extended" profile uses 200 (100 for heterogeneity). For a mean of per-dataset proportions, MCSE = SD/√n_sim ≤ 0.5/√n_sim. At n_sim = 100 the worst case is 0.05 and the typical case (SD ≈ 0.2) is 0.02; at 200 these are 0.035 and 0.014. AUC (SD ≈ 0.05–0.15) has an MCSE of 0.005–0.015 at n_sim = 100. `analysis/30_performance_tables.R` writes `_tables/precision_report.csv`. It gives the maximum MCSE actually achieved and the n_sim that a target MCSE would need (Morris Section 5.3, eq. 1), so you can quote both in Chapter 3.
- **Random numbers.** Each study has one master seed (`R/config.R`), set once with the L'Ecuyer-CMRG generator. Every dataset gets its own independent stream (`parallel::nextRNGStream`), stored in `design_full.rds` (Morris Section 4.1, 4.1.1). Datasets are therefore identical whatever the number of cores, chunks or order of completion. Any dataset can be regenerated with `regenerate_dataset()`. MIMOSA2 itself uses no random numbers, and `10_dgm_checks.R` section E verifies this.
- **Profiles.** "standard" is sized for about 20 cores in about a day. "extended" runs the full version 2 grids with more repetitions. The standard run is an exact subset of the extended run (same task IDs and streams), so you can run standard first and extend later without wasting anything.

Dataset counts are given below. Run times are planning estimates. They come from your logs: a P = 50 fit took a median of 15.8 s (`simulation_times.csv`), and a P = 100 imbalance task took about 40 s per worker (the 12,000-task run). The fitted time per fit is about 0.11 × P^1.27 seconds. `11_smoke_test.R` re-measures this on your cluster.

| Study | standard: scenarios × n_sim = datasets | extended: scenarios × n_sim = datasets |
|---|---|---|
| 1 Baseline | 468 × 100 = 46,800 | 468 × 200 = 93,600 |
| 2 Prior family | 100 × 100 = 10,000 | 300 × 200 = 60,000 |
| 3 Heterogeneity | 30 × 100 = 3,000 (2 fits each) | 405 × 100 = 40,500 (2 fits each) |
| 4 Imbalance | 57 × 100 = 5,700 | 120 × 200 = 24,000 |
| 5 maxit check | 12 × 50 = 600 (2 fits each) | 12 × 100 = 1,200 (2 fits each) |
| **Estimated time** | **~375 core-hours: ≈ 19 h on 20 cores, ≈ 9.5 h on 40** | **~1,650 core-hours: ≈ 26 h on 64 cores, ≈ 17 h on 100** |

---

## Study 1. Baseline sensitivity analysis (`sims/20_sim_baseline.R`)

**Aims.**
1. Measure how well MIMOSA2 identifies vaccine responders when its own assumptions hold (Beta-distributed proportions, binomial counts).
2. Map how performance depends on effect size, cell count, number of subjects and the proportion of responders.
3. Compare it with the DiD test on the same data.

This is a "proof-of-concept and limits" study in Morris's terms (Section 3.1).

**Data-generating mechanisms.** Beta family with φ = 2000. The design is fully factorial (Morris Section 3.2):

| Factor | Levels |
|---|---|
| Effect size δ | 1×10⁻³, 5×10⁻⁴, 2.5×10⁻⁴, 1.25×10⁻⁴, 6.25×10⁻⁵ |
| Cells per assay N | High 250,000; Medium 100,000; Low 15,000 |
| Subjects P | 10, 20, 30, 50, 75, 100 |
| Responders ρ | 0.10, 0.25, 0.50, 0.75, 0.90 |

This gives 450 scenarios. **Null scenarios (new): ρ = 0 for each P × N**, adding 18 scenarios, for 468 in total. With no responders, every call is a false positive. These are "methods with known properties" checks (Morris Section 4.2):

- MIMOSA2's Bayesian FDR and DiD-BH should give P(any call) ≤ α. When all subjects are null, the FDR equals the probability of at least one false call.
- The unadjusted DiD test should give a per-subject false-positive rate ≤ α.

Setting δ = 0 would not give a null here, because the profile constraints still force Δ > 0 for R1–R4.

Profiles: standard and extended both run all 468 scenarios, with 100 and 200 repetitions respectively.

**Targets.** Responder status of each subject, as in Section 0.3.

**Methods.** MIMOSA2 fitted to all P subjects, and DiD (unadjusted and BH).

**Performance measures.** All measures in Section 0.3. Main outputs:

- `base_plot` / `base_plot_clean`: TPR at 1% vs δ, MIMOSA2 and DiD-BH, ±1.96 MCSE.
- `base_specificity_*`, `base_fdr_*` (new): specificity and observed FDR at 1% and 5%.
- `base_auc`, `base_dAUC`: per-dataset AUC and paired ΔAUC.
- `ROC_plot_clean`: vertically averaged ROC curves at P = 50.
- `fdr_plot_clean`: calibration of the FDR at P = 50, using the correct q-value rule.
- `baseline_null_scenarios.csv`: type I error and FDR checks.
- `baseline_summary_P50.csv`: TPR, FDR and AUC (MCSE) for both methods side by side.

**What to look for.**

- Where TPR collapses: the low-N, small-δ corner, where the expected extra positive cells δ×N fall below about 2.
- Whether MIMOSA2's FDR stays at or below nominal when P is small or ρ is extreme. Here the mixture weights are poorly estimated, so the posterior can be overconfident.
- Whether ΔAUC favours MIMOSA2. Version 2's tables suggested it does by about 0.03–0.05, but that came from pooled AUCs.
- Whether the null scenarios behave as expected.

---

## Study 2. Robustness to the generative distribution (`sims/21_sim_prior.R`)

**Aims.** MIMOSA2 assumes the subject-level proportions are Beta. This study asks how much its performance changes when the proportions come from other distributions with **the same mean and variance**. Holding mean and variance fixed isolates the effect of distributional shape (skewness, behaviour near 0, bimodality) from the effect of heterogeneity. This is a robustness-to-misspecification study (Morris Section 3.1).

**Data-generating mechanisms.** Five families. Each is calibrated **at every condition mean** to the mean and variance of Beta(μφ, (1 − μ)φ) with φ = 10,000. `10_dgm_checks.R` section A verifies the calibration to within 1% (mean) and 3% (SD) from 10⁶ draws.

| Code | Family | Definition | Calibration |
|---|---|---|---|
| Beta | Beta (reference) | Beta(μφ, (1 − μ)φ) | — |
| EG | Exponential–gamma | p = e^(−G), G ~ Gamma(k, k/m) | m = k(μ^(−1/k) − 1) gives E[p] = μ exactly; k is solved for the variance (k ≈ 113 at μ = 10⁻⁴ and ≈ 527 at μ = 1.1×10⁻³) |
| LN | Logit-normal | logit p ~ N(ℓ, s²) | ℓ solved so that E[p] = μ, s solved for the variance (Gauss–Hermite quadrature); s ≈ 0.83 at 10⁻⁴ and ≈ 0.30 at 1.1×10⁻³ |
| SX | Inverse-Gaussian on the odds ("simplex-type") | odds ~ IG(μ/(1 − μ), φμ²(1 − μ)²) | matched by construction (variance = Beta variance to within 0.02%) |
| BB | Bimodal Beta | for stimulated draws: B ~ Bernoulli(0.25) chooses one of two Beta components; within-component variance = 1/8 of the total; separation set so the overall mean and variance equal the Beta ones. Unstimulated draws are Beta. | exact by construction |

The families are crossed fully factorially with:

| Factor | Levels |
|---|---|
| Effect size δ | 1×10⁻³, 2.5×10⁻⁴, 1.25×10⁻⁴, 6.25×10⁻⁵ |
| Cells per assay | 250,000; 100,000; 15,000; 7,000; 3,000 |
| Subjects P | 20, 50, 100 |

ρ = 0.5 throughout. This gives 300 scenarios. Standard: P = 50 only (100 scenarios × 100 repetitions). Extended: all 300 × 200.

**Targets and methods.** As in Section 0.

**Performance measures.**

- TPR at 1%, FDR at 5% and AUC by family (`prior_plot_P*`, `fdr_prior_plot_P*`, `prior_auc_P*`).
- **Difference from Beta** in the same (P, δ, N), with MCSE = √(MCSE₁² + MCSE₂²) (`prior_dAUC_vs_beta_P*`). This is the key robustness result.
- Vertically averaged ROC curves by family (`prior_ROC_comp_plot`).
- The density figure of the five families (`prior_density`).
- A figure showing why version 2 had to be re-run (`prior_v2_vs_v3`).

**What to look for.**

- Families whose density goes to 0 at p = 0 (LN, SX) put much less mass near zero than the Beta. This changes how often positive counts are zero at low N.
- The bimodal family breaks the "single responder distribution" assumption.
- Whether any of these changes MIMOSA2's FDR calibration (posterior over-confidence) more than its ranking (AUC).

---

## Study 3. Effect-size heterogeneity (`sims/22_sim_heterogeneity.R`)

**Aims.** MIMOSA2 assumes that all responders share one Beta distribution for p(s,1). If some responders respond much more strongly than others, the fitted responder component is pulled towards the strong responders, and weak responders may then look like non-responders. This study measures the cost to weak responders of being analysed together with strong ones.

**Data-generating mechanisms.** Each dataset combines two subgroups, simulated separately and then stacked:

- a **focal** subgroup of P_focal subjects with the smaller effect δ_focal;
- an **other** subgroup of P_other subjects with the larger effect δ_other.

Both subgroups use the same ρ, φ = 5000 and N.

| Factor | Levels |
|---|---|
| (P_focal, P_other) | (20, 80), (50, 50), (80, 20) |
| δ_focal | 5×10⁻⁴, 6.25×10⁻⁴ |
| δ_other (> δ_focal) | 6.25×10⁻⁴, 8×10⁻⁴, 1×10⁻³, 1×10⁻², 5×10⁻² |
| ρ | 0.2, 0.5, 0.8 |
| N | 250k, 100k, 15k, 7k, 3k |

There are 9 effect pairs, giving 405 scenarios. Standard: ρ = 0.5, δ_focal = 5×10⁻⁴, N ∈ {100k, 15k}, giving 30 scenarios × 100 repetitions.

> Version 2 had the two effects swapped in the code. The evaluated "small" group actually received the larger effect. This is fixed.

**Targets.** Responder status of the **focal** subjects.

**Methods.**

- **MIMOSA2, pooled fit:** fitted to all P_focal + P_other subjects. This is what an analyst would do.
- **MIMOSA2, separate fit:** fitted to the focal subjects only. This is the reference.
- **DiD:** per subject. BH adjustment is within the focal subjects.

**Performance measures.**

- AUC, TPR at 1% and FDR at 5% on the focal subjects for the three methods (`het_auc_*`, `het_tpr_*`, `het_fdr_*`).
- The **paired** difference pooled − separate in AUC per dataset, with MCSE (`het_dAUC_pooled_vs_separate`). This is the main answer: below 0 means sharing information with strong responders hurts weak responders.

The pooled-fit probabilities of the "other" subjects are also stored (`subjects$Group == "other"`), in case you want to show the effect on strong responders.

---

## Study 4. Cell-count imbalance (`sims/23_sim_imbalance.R`)

**Aims.** In real data the four assays of a subject rarely have the same number of cells. When one assay has fewer cells, the DiD estimate becomes noisier, and the beta-binomial likelihood of that assay becomes flatter in MIMOSA2. This study asks which method is hurt more. It also separates the effect of **imbalance** from the effect of simply having **fewer cells in total**.

**Data-generating mechanisms.** δ = 5×10⁻⁴, φ = 5000, P = 100 and ρ = 0.25 throughout. The scenarios are:

- **Reference:** 150,000 cells in every assay.
- **Imbalanced:** the assays in a "depletion pattern" are reduced to a fraction f of 150,000, and the others keep 150,000.
  - Patterns: s1, u1, s0, u0 (one assay); s1+s0, u1+u0, s1+u1, s0+u0, s1+u0, s0+u1 (two); s1+u1+u0, s1+u1+s0, s1+u0+s0, u1+u0+s0 (three).
  - f ∈ {0.90, 0.75, 0.50, 0.25, 0.10, 0.05, 0.025}.
- **Balanced controls:** for k depleted assays, the same total, (4 − k + kf)×150,000, is spread equally, so every assay has (4 − k + kf)/4 × 150,000 cells.

This gives 1 + 21 + 98 = 120 scenarios. Standard: reference, balanced k = 1 and 2, the four single-assay patterns, and s1+u1 and s0+u0, giving 57 scenarios × 100 repetitions. Extended: all 120 × 200.

**Targets and methods.** As in Section 0.

**Performance measures.**

- Mean per-dataset AUC against f for imbalanced, matched balanced and reference (`cell_imbalance_plot`).
- **Imbalanced − matched balanced** AUC with MCSE (`imbalance_vs_balanced`). This is the effect of imbalance at a fixed total number of cells.
- The AUC odds ratio against the reference, as in version 2, with a delta-method CI built from the MCSEs (`auc_or_plot`).

---

## Study 5. Is maxit = 30 enough? (NEW; `sims/24_sim_maxit_check.R`)

**Aims.**

1. Justify stopping MIMOSA2's EM algorithm at 30 iterations; the package default is 100. MIMOSA2 does not report whether it converged, so the only direct check is to refit with more iterations.
2. Demonstrate that the stored random-number states reproduce datasets exactly (Morris Section 4.1).

**Data-generating mechanisms.** Twelve baseline scenarios: ρ = 0.5; P ∈ {20, 50, 100}; δ ∈ {2.5×10⁻⁴, 6.25×10⁻⁵}; N ∈ {100k, 15k}. The datasets are regenerated from the **baseline** streams, so they are exactly the baseline datasets with the same Task_ID. The run uses 50 repetitions (standard) or 100 (extended).

**Methods.** MIMOSA2 with maxit = 30 and with maxit = 100 on each dataset, plus DiD.

**Performance measures** (`analysis/35_maxit_check.R`).

- Paired differences (maxit 100 − maxit 30) in AUC, TPR and FDP, with MCSE.
- Agreement of the 1% calls.
- The maximum change in any subject's P(responder).
- Fit times.
- The reproducibility check: the maxit = 30 results must equal the baseline results for the same Task_ID.

---

## Where the results go

| Location | Contents |
|---|---|
| `_simulations/<study>/tasks/task_XXXXXX.rds` | one file per dataset (resumable) |
| `_simulations/<study>/design_full.rds` | every planned dataset: scenario, repetition, Task_ID, **RNG stream seed** |
| `_simulations/<study>/<study>_results_<profile>.rds` | the combined tables `design`, `datasets`, `fits`, `estimates`, `subjects` (the "estimates" and "states" datasets of Morris Table 5) |
| `_tables/*.csv` | performance tables with MCSE, failure tables, realised DGM tables, precision report |
| `_fig/*.pdf` | all figures |
