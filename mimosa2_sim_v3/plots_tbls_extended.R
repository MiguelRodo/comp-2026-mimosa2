setwd("C:/GitHub/comp-2026-mimosa2/mimosa2_sim_v3")    # folder containing R/, analysis/, _simulations/
install.packages(c("dplyr", "tidyr", "purrr", "ggplot2", "knitr"))   # first time only

for (f in c("30_performance_tables", "31_plots_baseline", "32_plots_prior",
            "33_plots_heterogeneity", "34_plots_imbalance", "35_maxit_check")) {
  message("\n######## ", f)
  tryCatch(source(file.path("analysis", paste0(f, ".R")), local = new.env()),
           error = function(e) message("FAILED ", f, ": ", conditionMessage(e)))
  gc()
}


fig <- file.info(list.files("_fig", "\\.pdf$", full.names = TRUE))
fig$name <- basename(rownames(fig))
fig$today <- as.Date(fig$mtime) == Sys.Date()

expected <- c("base_plot", "base_plot_clean", "base_specificity_010", "base_specificity_050",
              "base_fdr_010", "base_fdr_050", "base_auc", "base_dAUC", "ROC_plot_clean", "fdr_plot_clean",
              "prior_density", "prior_v2_vs_v3", "prior_ROC_comp_plot",
              "het_dAUC_pooled_vs_separate",
              "cell_imbalance_plot", "imbalance_vs_balanced", "auc_or_plot")
cat("Missing:", setdiff(paste0(expected, ".pdf"), fig$name), "\n")

# figures made once per P or per scenario
for (pat in c("^prior_plot_P", "^fdr_prior_plot_P", "^prior_auc_P", "^prior_dAUC_vs_beta_P",
              "^het_auc_", "^het_tpr_", "^het_fdr_"))
  cat(sprintf("%-22s %d files\n", pat, sum(grepl(pat, fig$name))))

cat("\nNot made today (old leftovers):\n"); print(fig$name[!fig$today])
cat("\nSuspiciously small (< 5 KB, probably empty):\n"); print(fig$name[fig$size < 5000])

