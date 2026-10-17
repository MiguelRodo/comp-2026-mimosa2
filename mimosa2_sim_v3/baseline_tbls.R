# =============================================================================
# tbl_baseline_performance_fdr.R
# Study 1: observed FDR of MIMOSA2 at alpha = 0.05 by rho (rows) x P (columns).
# rho > 0: mean over the 15 Effect x Cell_range scenarios [max].
# rho = 0: mean over the 3 Cell_range null scenarios [max].
# Writes tbl_baseline_performance_fdr.txt and tbl_baseline_performance_fdr.tex
# =============================================================================

library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

perf_file <- "_tables/baseline_performance.csv"   # edit path if needed
out_txt   <- "tbl_baseline_performance_fdr.txt"
out_tex   <- "tbl_baseline_performance_fdr.tex"

perf <- read.csv(perf_file)   # numeric columns (FDR, FDR_mcse are separate)

# ---- Summarise ---------------------------------------------------------------
s <- perf |>
  filter(Method == "MIMOSA2", Alpha == 0.05) |>
  mutate(rho = as.numeric(sub("Prop_", "", Res_prop))) |>
  group_by(rho, P) |>
  summarise(mean     = mean(FDR),
            max      = max(FDR),
            mcse_max = max(FDR_mcse),
            n        = n(),
            .groups  = "drop")

# Check: 3 scenarios per cell for rho = 0, 15 otherwise
stopifnot(all(s$n[s$rho == 0] == 3), all(s$n[s$rho > 0] == 15))

max_mcse <- max(s$mcse_max)

# ---- Build the table -----------------------------------------------------------
P_levels <- c(10, 20, 30, 50, 75, 100)

tab <- s |>
  mutate(cell = sprintf("%.3f [%.3f]", mean, max)) |>
  select(rho, P, cell) |>
  pivot_wider(names_from = P, values_from = cell) |>
  arrange(rho)
tab <- tab[, c("rho", as.character(P_levels))]

row_lab <- ifelse(tab$rho == 0, "$\\rho = 0$ (null)", sprintf("$\\rho = %.2f$", tab$rho))
out <- data.frame(Row = row_lab, tab[, -1], check.names = FALSE)

write.table(out, out_txt, sep = "\t", row.names = FALSE, quote = FALSE)

# ---- LaTeX -----------------------------------------------------------------------
cap_long <- paste0(
  "Observed FDR of MIMOSA2 at $\\alpha = 0.05$ in Study~1 for each proportion of ",
  "responders and number of subjects. For $\\rho > 0$, each value is the mean ",
  "over the 15 combinations of effect size and cell count, with the largest ",
  "value in brackets. For the null scenarios ($\\rho = 0$), the mean is over the ",
  "three cell counts. Since every responder call in a null scenario is a false ",
  "positive, the observed FDR is equal to the proportion of datasets in which ",
  "at least one subject is classified as a responder. The largest MCSE is ",
  sprintf("%.3f", max_mcse), ".",
  "\\label{tbl_baseline_performance_fdr}")

kt <- kable(out, format = "latex", booktabs = TRUE, escape = FALSE,
            col.names = c("", as.character(P_levels)),
            align = c("l", rep("r", length(P_levels))),
            caption = cap_long,
            caption.short = "Study~1: observed FDR of MIMOSA2 at $\\alpha = 0.05$",
            linesep = c("\\midrule", rep("", nrow(out) - 1))) |>  # rule under null row
  add_header_above(c(" " = 1, "Number of subjects ($\\\\Nsub$)" = length(P_levels)),
                   escape = FALSE) |>
  kable_styling(latex_options = "hold_position")

tex <- as.character(kt)
tex <- sub("\\\\centering", "\\\\centering\n\\\\footnotesize", tex)  # smaller font

writeLines(tex, out_txt |> sub(pattern = "\\.txt$", replacement = ".tex"))
cat(tex)
