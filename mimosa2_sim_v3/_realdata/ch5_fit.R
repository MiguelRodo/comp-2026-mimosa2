## ch5_fit.R — Chapter 5: fit MIMOSA2 and DiD-BH once, save everything to ch5_results.rds
library(MIMOSA2)
library(DataTidy26VaccQFTDose)
library(tidyverse)
data("data_tidy_counts_ag", package = "DataTidy26VaccQFTDose")
dat <- data_tidy_counts_ag

MAXIT  <- 100      # MIMOSA2 EM iterations (package default)
DID_CC <- 0.5      # continuity correction in the DiD variance, as in the simulations
SEED   <- 20261026
OUT    <- "ch5_results.rds"

## One row per subject: four parent counts and four positive counts.
## Drops subjects with a zero parent count or more positive than parent cells.
prep_wide <- function(dataset_name, parent, pop_name,
                      stim_sel = "mtb", pre = "d0", post = "d224") {
  d <- dat |>
    filter(dataset == dataset_name, n_dose == "1", pop == pop_name,
           stim %in% c("uns", stim_sel), visit %in% c(pre, post))
  if (!is.na(parent)) d <- filter(d, parent_pop == parent)
  d |>
    mutate(assay = case_when(stim == stim_sel & visit == post ~ "s1",
                             stim == "uns"    & visit == post ~ "u1",
                             stim == stim_sel & visit == pre  ~ "s0",
                             stim == "uns"    & visit == pre  ~ "u0")) |>
    select(pid, qft, assay, parent_count, count) |>
    pivot_wider(names_from = assay, values_from = c(parent_count, count)) |>
    drop_na() |>
    filter(if_all(starts_with("parent_count"), ~ .x > 0),
           count_s1 <= parent_count_s1, count_u1 <= parent_count_u1,
           count_s0 <= parent_count_s0, count_u0 <= parent_count_u0)
}

## DiD Wald test with continuity-corrected variance (as DiD_wald() in the simulations)
did_test <- function(w, cc = DID_CC) {
  v   <- function(x, N) { pt <- (x + cc) / (N + 2 * cc); pt * (1 - pt) / N }
  est <- (w$count_s1 / w$parent_count_s1 - w$count_u1 / w$parent_count_u1) -
    (w$count_s0 / w$parent_count_s0 - w$count_u0 / w$parent_count_u0)
  se  <- sqrt(v(w$count_s1, w$parent_count_s1) + v(w$count_u1, w$parent_count_u1) +
                v(w$count_s0, w$parent_count_s0) + v(w$count_u0, w$parent_count_u0))
  z   <- ifelse(se == 0, 0, est / se)
  list(est = est, p = pnorm(z, lower.tail = FALSE))
}

## Fit one population; return one row per subject with both methods' results
fit_one <- function(dataset, parent, pop) {
  w <- prep_wide(dataset, parent, pop); P <- nrow(w)
  parent_lab <- ifelse(is.na(parent), "cd4", parent)
  cat(sprintf("%-10s %-26s P = %2d\n", parent_lab, pop, P))
  
  Ntot <- as.matrix(w[, c("parent_count_s1", "parent_count_s0",
                          "parent_count_u1", "parent_count_u0")])
  colnames(Ntot) <- c("ns1", "ns0", "nu1", "nu0")
  fit <- if (P >= 5) tryCatch(
    MIMOSA2(Ntot, w$count_s1, w$count_u1, w$count_s0, w$count_u0, maxit = MAXIT),
    error = function(e) NULL) else NULL
  ok <- !is.null(fit)
  
  d    <- did_test(w)
  p_bh <- p.adjust(d$p, method = "BH")
  
  tibble(parent_pop = rep(parent_lab, P), pop = rep(pop, P), pid = w$pid,
         s1 = w$count_s1, u1 = w$count_u1, s0 = w$count_s0, u0 = w$count_u0,
         Ns1 = w$parent_count_s1, Nu1 = w$parent_count_u1,
         Ns0 = w$parent_count_s0, Nu0 = w$parent_count_u0,
         mim_prob    = if (ok) rowSums(fit$z[, 1:4, drop = FALSE]) else rep(NA_real_, P),
         mim_resp_01 = if (ok) as.logical(getResponse(fit, 0.01)) else rep(NA, P),
         mim_resp_05 = if (ok) as.logical(getResponse(fit, 0.05)) else rep(NA, P),
         did_est     = d$est,
         did_p_bh    = p_bh,
         did_resp_01 = p_bh <= 0.01,
         did_resp_05 = p_bh <= 0.05)
}

## All populations: CD4+ and every phenotypic subset
pop_list <- bind_rows(
  tibble(dataset = "any_cyt_parent_cd4", parent = NA_character_, pop = "any_cyt"),
  dat |>
    filter(dataset == "any_cyt_pheno_parent_single") |>
    distinct(parent = parent_pop, pop) |>
    mutate(dataset = "any_cyt_pheno_parent_single"))

set.seed(SEED)
subjects <- pmap_dfr(pop_list, fit_one)

## Summary table (same columns as chapter5_full_table)
tab <- subjects |>
  group_by(parent_pop, pop) |>
  summarise(P = n(),
            median_count_s1  = median(s1),  median_parent_s1 = median(Ns1),
            median_count_u1  = median(u1),  median_parent_u1 = median(Nu1),
            median_count_s0  = median(s0),  median_parent_s0 = median(Ns0),
            median_count_u0  = median(u0),  median_parent_u0 = median(Nu0),
            rho_hat       = mean(mim_prob),
            n_called_05   = sum(mim_resp_05), n_called_01 = sum(mim_resp_01),
            n_did_05      = sum(did_resp_05), n_did_01    = sum(did_resp_01),
            rho_did_bh_01 = mean(did_resp_01),
            median_raw_did = median(did_est),
            .groups = "drop")
tab <- pop_list |>
  transmute(parent_pop = coalesce(parent, "cd4"), pop) |>
  left_join(tab, by = c("parent_pop", "pop"))

## What was run, so the results can be reproduced
run_info <- list(
  date                  = Sys.time(),
  R                     = R.version.string,
  MIMOSA2               = as.character(packageVersion("MIMOSA2")),
  MIMOSA2_arguments     = names(formals(MIMOSA2)),
  DataTidy26VaccQFTDose = as.character(packageVersion("DataTidy26VaccQFTDose")),
  data_commit           = packageDescription("DataTidy26VaccQFTDose")$RemoteSha,
  MAXIT = MAXIT, DID_CC = DID_CC, SEED = SEED)

saveRDS(list(run_info = run_info, table = tab, subjects = subjects), OUT)
cat("\nSaved", OUT, "\n")
print(tab, n = Inf, width = Inf)