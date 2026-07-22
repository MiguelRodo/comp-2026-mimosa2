# Assessing MIMOSA2

This repository provides infastructure for the project to assess MIMOSA2.

## Contact

For more information, please contact:  
- Miguel J. Rodo:
  - miguel.rodo@uct.ac.za
  - miguel.rodo@outlook.com
- 

## Links

- [URLs to data sources (e.g. OneDrive), GitHub repositories, publications, etc.]

## Details

### Installing `R` packages

```r
if (!requireNamespace("pak", quietly = TRUE)) {
  utils::install.packages("pak")
}
pak::local_install_dev_deps()
```

### HPC

#### Running `R` interactively
`
Once you have cloned the repo down to `/scratch/$USER/projects/mimosa2`:

- Either, start a VS Code session in that directory (via `ondemand.uct.ac.za`), or `cd` to it.
- Open a terminal: `Ctrl + Shift + backtick`
- Run `apptainer-run -f mimosa` in that terminal
- Run `R`

#### Running `R` in batch mode

- Log into hpc using powershell: `ssh <username>@hex.uct.ac.za`
- Switch to project directory: `cd /scratch/$USER/projects/mimosa2`
- Run `slurm-sbatch <path/to/script>`, e.g. `slurm-sbatch scripts/slurm/test.sh`
- To check job is running, run `slurm-squeue`
- To check job history, look in `_tmp/log/sbatch/<name_of_script>`
  - `run_<date>_<time>_`:
    - `out.txt`
      - Job output
    - `script.txt`
      - Actual script executed


