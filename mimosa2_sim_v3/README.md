# MIMOSA2 simulation study: version 3 (re-run version)

Isabella Lethbridge & Tayyeb Abrahams. Prepared 2–3 October 2026 for the re-run on Saturday 3 October.

This folder replaces the version-2 simulation and plotting scripts. It keeps your designs (the same factors and levels). It fixes the errors found in the audit and brings the study in line with Morris, White & Crowther (2019).

*Morris, White & Crowther (2019) - Using simulation studies to evaluate statistical methods* 

There are three documents:

| File | What it is for |
|---|---|
| `README.md` (this file) | how to run everything tomorrow |
| `STUDY_DESCRIPTIONS.md` | what each of the five studies does, in ADEMP order (use for Chapter 3) |
| `CHANGES_AND_AUDIT.md` | file-by-file audit of version 2, what changed and why, a Morris checklist, limitations and write-up notes |

> **Important: this code has not been run yet.** R could not be installed where it was written. The numerical parts (calibration, allocation, AUC, ROC averaging) were checked by re-implementing them in Python. The R code was then reviewed line by line, independently, and the problems found were fixed. **Run the smoke test first (step 2).** It takes about 15–30 minutes, exercises every script end to end on a tiny grid, and stops with a clear FAIL message if anything is wrong.

---

## Folder layout

```
mimosa2_sim_v3/
  R/config.R          profile, seeds, nsim, time limit, paths      <- the only file with settings
  R/dgm.R             data-generating mechanism (fixed)
  R/methods.R         DiD, MIMOSA2 wrapper with a time limit, per-dataset measures
  R/runner.R          seeds/streams, parallel scheduler, resume, combine
  R/scenarios.R       the scenario grids of all studies
  sims/10_dgm_checks.R      checks the DGM (Morris 4.2)                     ~5 min
  sims/11_smoke_test.R      checks everything + projects run times           ~15-30 min
  sims/20_sim_baseline.R    Study 1
  sims/21_sim_prior.R       Study 2
  sims/22_sim_heterogeneity.R  Study 3
  sims/23_sim_imbalance.R   Study 4
  sims/24_sim_maxit_check.R Study 5 (new, small)
  sims/29_combine.R         merges task files (only needed after job arrays)
  analysis/30_performance_tables.R   all tables with Monte Carlo SEs
  analysis/31-35_*.R                 figures + tables per study
  slurm/*.sh                job scripts (edit the #SBATCH and module lines)
```

Always run scripts **from inside `mimosa2_sim_v3/`**, e.g. `Rscript sims/20_sim_baseline.R`. On the cluster, copy the folder into your repository or scratch space. Results go to `_simulations/`, `_tables/` and `_fig/` inside the folder.

---

## Run sheet

### Step 0. Copy and edit (5 min)

1. Copy `mimosa2_sim_v3/` to the cluster.
2. In `slurm/*.sh`, edit the `#SBATCH --account/--partition` lines and the `module load` line. Copy them from the job script you used for version 2.
3. If your R library is not `/scratch/abrmoe030/R_libs`, run `export MIMOSA2_R_LIBS=/path/to/R_libs` (or edit `R/config.R`).

### Step 1. Check the DGM and the pipeline (about 30 min)

```bash
sbatch slurm/smoke_test.sh        # runs sims/10_dgm_checks.R then sims/11_smoke_test.R
```

Read the end of `slurm_mimosa2_smoke_<jobid>.out`. It should finish with **"ALL DGM CHECKS PASSED"** and **"All checks passed"**.

The smoke test also prints a **projected run-time table** for both profiles, measured on your cluster. Use it to choose an option in step 2.

When it has passed, delete the `_smoke/` folder.

If something fails, the message names the check. The logs of each script are in `_smoke/*.log`.

### Step 2. Run the simulations

**Option 1: about 20 cores (standard profile, ~375 core-hours ≈ 19 h).**

If you can run several jobs at the same time, submit each study separately. This is faster: the baseline alone takes about 11 h.

```bash
sbatch slurm/run_study.sh prior
sbatch slurm/run_study.sh imbalance
sbatch slurm/run_study.sh heterogeneity
sbatch slurm/run_study.sh maxit_check
sbatch slurm/run_study.sh baseline
```

If you can run only one job at a time:

```bash
sbatch slurm/run_all_standard.sh      # all studies, smallest first, baseline last
```

**Option 2: more cores (extended profile, ~1,650 core-hours ≈ 26 h on 64 cores).**

Use job arrays, where each array element is one 20-core job:

```bash
SIM_PROFILE=extended sbatch --array=0-3 slurm/run_study_array.sh baseline
SIM_PROFILE=extended sbatch --array=0-5 slurm/run_study_array.sh prior
SIM_PROFILE=extended sbatch --array=0-7 slurm/run_study_array.sh heterogeneity
SIM_PROFILE=extended sbatch --array=0-3 slurm/run_study_array.sh imbalance
SIM_PROFILE=extended sbatch --array=0-1 slurm/run_study_array.sh maxit_check
# when the arrays have finished:
SIM_PROFILE=extended sbatch slurm/combine.sh baseline prior heterogeneity imbalance maxit_check
```

You can also run standard first and extended later. The standard datasets are an exact subset of the extended ones (same task IDs and random-number streams), so the extended run only computes the tasks that are missing. The two profiles' results are combined into separate files (`*_results_standard.rds`, `*_results_extended.rds`).

**If a job hits its time limit:** submit exactly the same command again. Finished tasks are skipped.

**If you must stop early** (for example the baseline on Sunday morning): `scancel` the job, then run

```bash
SIM_PROFILE=standard Rscript sims/29_combine.R baseline
```

Tasks run repetition by repetition, so the combine step keeps the repetitions completed in every scenario (for example 1–85) and tells you the n_sim it used. Report that n_sim and its MCSE.

**Progress:** `_simulations/<study>/progress_chunkXX.log`. The `.out` file also prints the tasks per hour and the estimated hours left every 500 tasks.

### Step 3. Tables and figures

On the cluster, or locally after copying `_simulations/<study>/<study>_results_<profile>.rds`:

```bash
SIM_PROFILE=standard Rscript analysis/30_performance_tables.R   # all tables + precision_report.csv
SIM_PROFILE=standard Rscript analysis/31_plots_baseline.R
SIM_PROFILE=standard Rscript analysis/32_plots_prior.R
SIM_PROFILE=standard Rscript analysis/33_plots_heterogeneity.R
SIM_PROFILE=standard Rscript analysis/34_plots_imbalance.R
SIM_PROFILE=standard Rscript analysis/35_maxit_check.R
```

In RStudio on Windows: `setwd("C:/.../mimosa2_sim_v3")`, then `Sys.setenv(SIM_PROFILE = "standard")`, then `source("analysis/31_plots_baseline.R")`. The analysis scripts do not fork, so they work on Windows. The simulation scripts also run on Windows, sequentially and without the time limit, which is useful only for testing.

---

## What changed, in one list

See `CHANGES_AND_AUDIT.md` for the details and the evidence.

1. **Responder allocation fixed.** Realised ρ was up to 0.60 for a nominal 0.25.
2. **Prior families matched in mean and variance at every condition mean.** In version 2, EG and LN were 2–3× more variable than Beta at the responder mean, and the LN mean was 1.28μ.
3. **Effect-size heterogeneity: effects were swapped.** The evaluated group had the larger effect.
4. **Imbalance plots: every scenario was mislabelled or dropped** by the name matching.
5. **FDR calibration used the wrong rule** (local posterior instead of q-values) and pooled datasets.
6. **AUC computed per dataset**, not by pooling subjects across datasets and P. The invalid Hanley–McNeil ΔAUC intervals are replaced by paired per-dataset differences with MCSE.
7. **Monte Carlo SEs for every measure.** n_sim = 100/200 is justified by MCSE, and a precision report is produced.
8. **Seeds:** one master seed per study, one stored random-number stream per dataset, reproducible and resumable.
9. **One time limit (600 s per fit) everywhere.** Failures are stored with reasons and reported.
10. **New:** null scenarios (ρ = 0), the DiD-BH comparator (FDR-controlling), specificity and FDR figures, realised-Δ tables, and Study 5 (is maxit = 30 enough? plus a reproducibility check).
11. **DiD computed in closed form** (the same statistic as the GLM). This removes the ties at 0.05 for negative or non-converged subjects.
12. **Smaller fixes:** "Extremely Low" no longer dropped from the prior plots, axis-label typos, swapped table columns, undefined objects.

---

## Re-creating a single dataset (e.g. one where MIMOSA2 failed)

```r
setwd("mimosa2_sim_v3"); for (f in c("R/config.R","R/dgm.R","R/methods.R","R/runner.R","R/scenarios.R")) source(f)
res <- readRDS("_simulations/baseline/baseline_results_standard.rds")
bad <- subset(res$fits, Status != "ok")[1, ]
row <- subset(res$design, Task_ID == bad$Task_ID)
sim <- regenerate_dataset(row, simulate_baseline)   # exactly the dataset of that task
library(MIMOSA2); fit <- MIMOSA2(Ntot = sim$Ntot, ns1 = sim$ns1, nu1 = sim$nu1, ns0 = sim$ns0, nu0 = sim$nu0, maxit = 30)
```

Use `simulate_prior`, `simulate_heterogeneity` or `simulate_imbalance` for the other studies.
