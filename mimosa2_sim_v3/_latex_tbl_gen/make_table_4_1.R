library(knitr)
library(kableExtra)

# ---- Settings ----
csv_file <- "baseline_summary_P50.csv"
out_dir  <- "."          # folder for the .tex files ("." = working directory)
rho      <- "Prop_0.50"
digits   <- 3

# Label stem from the CSV file name, e.g. "tbl_baseline_summary_P50"
label_stem <- paste0("tbl_", tools::file_path_sans_ext(basename(csv_file)))

# ---- Read the CSV (as text, so "1.2e-04" etc. stay as labels) ----
d <- read.csv(csv_file, colClasses = "character", check.names = FALSE)
d <- d[d$Res_prop == rho, ]

est  <- function(x) as.numeric(sub(" .*", "", x))         # "0.842 (0.007)" -> 0.842
mcse <- function(x) as.numeric(gsub(".*\\(|\\)", "", x))  # "0.842 (0.007)" -> 0.007
fmt  <- function(x) formatC(x, format = "f", digits = digits)

effects <- c("1.0e-03", "5.0e-04", "2.5e-04", "1.2e-04", "6.2e-05")
labels  <- c("$10^{-3}$", "$5\\times10^{-4}$", "$2.5\\times10^{-4}$",
             "$1.25\\times10^{-4}$", "$6.25\\times10^{-5}$")
cells   <- c("High", "Medium", "Low")
rho_txt <- sub("Prop_", "", rho)

# ---- Function: write the LaTeX table for one measure to <label>.tex ----
# 'what' is the last part of the label, e.g. "sensitivity" -> tbl_<csv>_sensitivity
write_table <- function(m2_col, did_col, measure_txt, what) {
  tab <- data.frame(delta = labels, stringsAsFactors = FALSE)
  for (cc in cells) {
    m2 <- did <- character(length(effects))
    for (k in seq_along(effects)) {
      r <- d[d$Effect == effects[k] & d$Cell_range == cc, ]
      if (nrow(r) != 1) stop("Expected one row for ", effects[k], ", ", cc)
      m2[k]  <- fmt(est(r[[m2_col]]))
      did[k] <- fmt(est(r[[did_col]]))
    }
    tab[[paste(cc, "M2")]]  <- m2
    tab[[paste(cc, "DiD")]] <- did
  }
  names(tab) <- c("$\\delta$", rep(c("MIMOSA2", "DiD-BH"), length(cells)))
  
  lab     <- paste0(label_stem, "_", what)
  max_se  <- max(mcse(d[[m2_col]]), mcse(d[[did_col]]))
  caption <- sprintf(paste0(
    "Study~1 ($\\Nsub = 50$, $\\rho = %s$): %s of MIMOSA2 and DiD-BH by nominal ",
    "effect size and cell count, averaged over 100 datasets. All MCSEs are at ",
    "most %.3f.\\label{%s}"), rho_txt, measure_txt, max_se, lab)
  
  latex <- kable(tab, format = "latex", booktabs = TRUE, escape = FALSE,
                 align = "lcccccc", caption = caption,
                 row.names = FALSE) |>
    add_header_above(c(" " = 1, "High" = 2, "Medium" = 2, "Low" = 2)) |>
    kable_styling(latex_options = "hold_position")
  
  out_file <- file.path(out_dir, paste0(lab, ".tex"))
  writeLines(as.character(latex), out_file)
  cat("Written:", out_file, "\n")
  cat("  Insert with:    \\input{", lab, ".tex}\n", sep = "")
  cat("  Reference with: Table~\\ref{", lab, "}\n\n", sep = "")
}

# ---- Sensitivity table (main text) ----
write_table("TPR_MIMOSA2", "TPR_DiD (BH)", "sensitivity at $\\alpha = 0.01$",
            "sensitivity")

# ---- AUC table (appendix, optional) ----
write_table("AUC_MIMOSA2", "AUC_DiD (BH)", "AUC", "auc")