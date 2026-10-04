# =============================================================================
# dgm.R : data-generating mechanism (DGM)
# Replaces simulate_MIMOSA2_alt_prior() in Non_beta_simulations.R
# =============================================================================
#
# The structure of the DGM is unchanged: eight response profiles (R1-R4 are
# vaccine responders; NR1, NR2, NR3, NSR are non-responders), condition means
# built from mu_0, delta_0, gamma and delta, rejection sampling for the
# inequality constraints, and binomial cell counts. The changes are:
#
# [CHANGE 1] RESPONDER ALLOCATION (bug; audit A1).
#   Version 2 used n = floor(P * pis) and then added the missing subjects one
#   at a time to the smallest component, which is always R1/R2 first. The
#   realised proportion of responders was therefore not the nominal one: e.g.
#   at P = 10 the "25%" and "75%" designs were both 60% responders, and at
#   P = 20 "10%" was 20%. Now the number of responders is rounded ONCE,
#   n_R = round(rho * P) (halves rounded up), and the responders and
#   non-responders are then split across their four profiles, with any
#   remainder given to randomly chosen profiles. The realised n_R is stored
#   for every dataset so it can be reported.
#
# [CHANGE 2] VARIANCE MATCHING AT EVERY CONDITION MEAN (validity; audit A7).
#   Version 2 calibrated the exponential-gamma (k = 113) and logit-normal
#   (phi = 2.05) families once, at mu_0 = 1e-4, and then used the same
#   dispersion for every condition mean. Checked numerically (see
#   sims/10_dgm_checks.R): at mu = 1.1e-3 (the stimulated responder mean for
#   delta = 1e-3) the Beta(phi = 10000) SD is 3.3e-4, but the old EG SD was
#   7.6e-4 (2.3x) and the old LN SD was 1.1e-3 (3.3x). The non-Beta families
#   were therefore much MORE heterogeneous than the Beta reference, so the
#   prior-family comparison mixed up "different shape" with "more variance".
#   Now each family is calibrated separately at each condition mean so that
#   its mean AND variance equal those of Beta(mu * phi, (1 - mu) * phi).
#
# [CHANGE 3] LOGIT-NORMAL MEAN (bug; audit A7).
#   Version 2 put the logit-scale location at logit(mu), which makes the
#   MEDIAN equal to mu; the mean was 1.28 * mu. The location is now solved so
#   that the mean equals mu.
#
# [CHANGE 4] BIMODAL BETA APPLIED CONSISTENTLY.
#   Version 2 drew the stimulated post-vaccination proportion from the bimodal
#   mixture, but the stimulated baseline proportion p(s,0) was bimodal in
#   profiles R3 and NSR and Beta in R1, R4 and NR2. Now every stimulated draw
#   (p(s,1) and p(s,0)) is bimodal and every unstimulated draw is Beta, which
#   matches the description in Chapter 3 ("applied to the stimulated
#   proportions"). Where a single draw is shared by a stimulated and an
#   unstimulated condition (NR1, NR2 share p(s,1) = p(u,1) at the
#   unstimulated mean; NR3 shares all four) it is treated as unstimulated.
#
# [CHANGE 5] The 'simplex' sampler (inverse-Gaussian on the odds scale) is
#   unchanged: with lambda = phi * mu^2 * (1 - mu)^2 its variance is
#   mu / (phi (1 - mu)^5), which equals the Beta variance to within 0.02% at
#   every mean used here, and its mean is 0.9999 mu (checked numerically).
#   It is NOT the Barndorff-Nielsen & Jorgensen simplex distribution; call it
#   "inverse-Gaussian (odds)" or "simplex-type" in the write-up.
#
# [CHANGE 6] Unused families (uniform, odds-exponential, odds-gamma,
#   unit-lognormal) were removed, and the large if/else blocks were replaced
#   by one sampler draw_prop(), so each profile is written out once. The
#   constraints and the order of random draws inside each profile are as in
#   version 2 (and as in MIMOSA2::simulate_MIMOSA2).
#
# [CHANGE 7] Safety: every rejection loop stops with an informative error
#   after 1e6 rounds instead of hanging, and the components check uses a
#   tolerance (version 2 used sum(components) != 1, which can fail for
#   floating-point sums).
#
# [CHANGE 9] STUCK SUBJECTS IN R1 AND R4 (found by the smoke test, 4 Oct).
#   In R1 and R4 the unstimulated proportions are drawn once and only the
#   stimulated ones are redrawn (as in MIMOSA2::simulate_MIMOSA2). If a
#   subject happens to get a large p(u,1), the constraints can be (almost)
#   impossible to satisfy: with the bimodal family the stimulated draws have
#   light tails, and at delta = 6.25e-5 about 0.7% of R1 subjects had
#   acceptance probability 0 (checked by simulation), so a dataset of 50
#   subjects hung in about 1 in 20 cases. Version 2 would have looped
#   forever. Now, if a subject is still rejected after 1e4 rounds, its
#   unstimulated proportion(s) are redrawn as well. This only touches
#   subjects that are essentially infeasible (it practically never happens
#   for Beta, EG, LN, SX), and the number of such redraws is stored per
#   dataset (n_stuck_redraws) so it can be reported.
#
# [CHANGE 8] The function also returns the true proportions and the true
#   difference-in-differences of every subject (true_delta), so the realised
#   effect size can be reported (rejection sampling makes it larger than the
#   nominal delta; Morris 3.2, pitfall 4).
# =============================================================================

PROFILE_LABELS   <- c("R1", "R2", "R3", "R4", "NR1", "NR2", "NR3", "NSR")
RESPONDER_LABELS <- c("R1", "R2", "R3", "R4")

# ---- Family names --------------------------------------------------------------
normalise_family <- function(prior) {
  p <- tolower(gsub("[ _]", "", prior))
  if (p %in% c("beta", "b"))                                          return("Beta")
  if (p %in% c("eg", "exponentialgamma"))                             return("EG")
  if (p %in% c("ln", "logitnormal", "logitnorm"))                     return("LN")
  if (p %in% c("sx", "simplex", "s", "invgauss", "ig"))               return("SX")
  if (p %in% c("bb", "bernoullibeta", "bimodal", "bimodalbeta"))      return("BB")
  stop("Unknown prior family: ", prior)
}

beta_var <- function(mu, phi) mu * (1 - mu) / (1 + phi)

# =============================================================================
# Calibration of the non-Beta families (mean and variance matched to Beta)
# =============================================================================

# ---- Exponential-gamma: p = exp(-G), G ~ Gamma(shape k, rate k / m) ----------
# E[p] = (1 + m/k)^(-k) = mu  <=>  m/k = mu^(-1/k) - 1   (exact)
# E[p^2] = (1 + 2m/k)^(-k)
# Var[p] = mu^2 * expm1( k * (2 log1p(r) - log1p(2r)) ),  r = m/k
# (written this way to avoid cancellation when Var is tiny)
eg_var <- function(mu, k) {
  r <- expm1(-log(mu) / k)
  mu^2 * expm1(k * (2 * log1p(r) - log1p(2 * r)))
}
calib_eg <- function(mu, phi_beta) {
  target <- beta_var(mu, phi_beta)
  f <- function(lk) log(eg_var(mu, exp(lk))) - log(target)
  lk <- uniroot(f, lower = log(0.5), upper = log(1e9), tol = 1e-12)$root
  k <- exp(lk)
  list(k = k, m = k * expm1(-log(mu) / k))
}
rgamma_exp <- function(n, par) {
  if (n <= 0) return(numeric(0))
  exp(-rgamma(n, shape = par$k, rate = par$k / par$m))
}

# ---- Logit-normal: logit(p) ~ N(loc, s^2) -----------------------------------
# Moments by Gauss-Hermite quadrature (probabilists' version, 200 nodes),
# nodes and weights from the Golub-Welsch algorithm (no extra packages).
.gh <- local({
  n <- 200
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1))
  J[cbind(1:(n - 1), 2:n)] <- off
  J[cbind(2:n, 1:(n - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  w <- e$vectors[1, ]^2
  list(x = e$values, w = w / sum(w))
})
ln_moments <- function(loc, s) {
  p <- plogis(loc + s * .gh$x)
  c(m1 = sum(.gh$w * p), m2 = sum(.gh$w * p^2))
}
ln_location <- function(mu, s) {
  f <- function(loc) ln_moments(loc, s)[["m1"]] - mu
  uniroot(f, lower = qlogis(mu) - 5 * s - 5, upper = qlogis(mu) + 1, tol = 1e-13)$root
}
calib_ln <- function(mu, phi_beta) {
  target <- beta_var(mu, phi_beta)
  f <- function(s) {
    loc <- ln_location(mu, s)
    m <- ln_moments(loc, s)
    log(m[["m2"]] - m[["m1"]]^2) - log(target)
  }
  s <- uniroot(f, lower = 1e-3, upper = 5, tol = 1e-12)$root
  list(sd = s, loc = ln_location(mu, s))
}
rlogitnorm <- function(n, par) {
  if (n <= 0) return(numeric(0))
  plogis(rnorm(n, mean = par$loc, sd = par$sd))
}

# ---- Simplex-type: odds X ~ inverse-Gaussian(mu/(1-mu), lambda) ----------------
# (unchanged from version 2; already matched, see header CHANGE 5)
rsimplex <- function(n, mu, phi_beta) {
  if (n <= 0) return(numeric(0))
  mu_ig     <- mu / (1 - mu)
  lambda_ig <- phi_beta * mu^2 * (1 - mu)^2
  v <- rnorm(n)^2
  x <- mu_ig + (mu_ig^2 * v) / (2 * lambda_ig) -
    (mu_ig / (2 * lambda_ig)) * sqrt(4 * mu_ig * lambda_ig * v + mu_ig^2 * v^2)
  z <- runif(n)
  idx <- z > (mu_ig / (mu_ig + x))
  x[idx] <- mu_ig^2 / x[idx]
  x / (1 + x)
}

# ---- Bimodal Beta: two-component Beta mixture with Beta mean and variance ------
# (unchanged from version 2 except that it is used consistently, CHANGE 4)
# B ~ Bernoulli(0.25); within-component variance = 1/8 of the target variance;
# separation chosen so the marginal mean and variance equal the Beta ones.
rbimodal_beta <- function(n, mu, phi_beta, p_bernoulli = 0.25) {
  if (n <= 0) return(numeric(0))
  var_marginal <- beta_var(mu, phi_beta)
  var_cond     <- 0.125 * var_marginal
  delta        <- sqrt((var_marginal - var_cond) / (p_bernoulli * (1 - p_bernoulli)))
  mu1 <- mu + (1 - p_bernoulli) * delta
  mu0 <- mu - p_bernoulli * delta
  if (mu0 <= 0 || mu1 >= 1) stop(sprintf("Bimodal Beta infeasible at mu = %g, phi = %g", mu, phi_beta))
  phi1 <- mu1 * (1 - mu1) / var_cond - 1
  phi0 <- mu0 * (1 - mu0) / var_cond - 1
  z <- rbinom(n, 1, p_bernoulli)
  out <- numeric(n)
  n1 <- sum(z == 1)
  if (n1 > 0)     out[z == 1] <- rbeta(n1, mu1 * phi1, (1 - mu1) * phi1)
  if (n - n1 > 0) out[z == 0] <- rbeta(n - n1, mu0 * phi0, (1 - mu0) * phi0)
  out
}

# ---- Calibration cache (computed once per family x mean x phi) ---------------
# Call prewarm_calibration() in the parent R process BEFORE forking workers so
# that every worker inherits the solved parameters.
.calib_cache <- new.env(parent = emptyenv())
family_param <- function(family, mu, phi_beta) {
  key <- sprintf("%s|%.12g|%.12g", family, mu, phi_beta)
  val <- .calib_cache[[key]]
  if (is.null(val)) {
    val <- switch(family,
                  EG = calib_eg(mu, phi_beta),
                  LN = calib_ln(mu, phi_beta),
                  NULL)
    assign(key, val, envir = .calib_cache)
  }
  val
}
condition_means <- function(effect,
                            bg_effect = BG_EFFECT,
                            baseline_stim_effect = BASELINE_STIM_EFFECT,
                            baseline_background = BASELINE_BACKGROUND) {
  MU0 <- baseline_background
  MS0 <- baseline_stim_effect + MU0
  MU1 <- MU0 + bg_effect
  MS1 <- MU1 + effect
  c(u1 = MU1, u0 = MU0, s0 = MS0, s1 = MS1)
}
prewarm_calibration <- function(families, effects, phis) {
  for (f in unique(vapply(families, normalise_family, ""))) {
    if (!f %in% c("EG", "LN")) next
    for (e in unique(effects)) for (ph in unique(phis)) {
      mus <- condition_means(e)
      for (mu in unique(mus)) family_param(f, mu, ph)
    }
  }
  invisible(length(ls(.calib_cache)))
}

# ---- One sampler for all families -------------------------------------------
# cond is "u1", "u0", "s0" or "s1" (which condition the draw is for).
draw_prop <- function(n, mu, cond, family, phi_beta) {
  if (n <= 0) return(numeric(0))
  switch(family,
         Beta = rbeta(n, mu * phi_beta, (1 - mu) * phi_beta),
         EG   = rgamma_exp(n, family_param("EG", mu, phi_beta)),
         LN   = rlogitnorm(n, family_param("LN", mu, phi_beta)),
         SX   = rsimplex(n, mu, phi_beta),
         BB   = if (cond %in% c("s1", "s0")) rbimodal_beta(n, mu, phi_beta)
                else rbeta(n, mu * phi_beta, (1 - mu) * phi_beta),
         stop("unknown family ", family))
}

# =============================================================================
# Allocation of subjects to profiles  [CHANGE 1]
# =============================================================================
split_group <- function(n, w) {
  k <- length(w)
  if (n <= 0 || sum(w) <= 0) return(rep(0L, k))
  share <- n * w / sum(w)
  base  <- floor(share + 1e-9)
  rem   <- n - sum(base)
  if (rem > 0) {
    frac <- share - base
    ord  <- order(-round(frac, 9), runif(k))   # largest remainder, random tie-break
    base[ord[seq_len(rem)]] <- base[ord[seq_len(rem)]] + 1
  }
  as.integer(base)
}
allocate_profiles <- function(P, components) {
  rho <- round(sum(components[1:4]), 10)
  n_r <- floor(rho * P + 0.5)                   # round half up, once
  n_n <- P - n_r
  n <- c(split_group(n_r, components[1:4]), split_group(n_n, components[5:8]))
  if (n_n > 0 && sum(components[5:8]) <= 0) stop("Non-responders requested but non-responder weights are zero")
  stopifnot(sum(n) == P)
  n
}

# =============================================================================
# Main simulator
# =============================================================================
simulate_MIMOSA2_alt_prior <- function(effect = 5e-4,
                                       bg_effect = BG_EFFECT,
                                       baseline_stim_effect = BASELINE_STIM_EFFECT,
                                       baseline_background = BASELINE_BACKGROUND,
                                       phi = 10000,
                                       P = 100,
                                       rng = c(100000, 150000),
                                       prior = "Beta",
                                       components = rep(1/8, 8),
                                       max_rounds = 1e6,
                                       stuck_rounds = 1e4) {
  if (effect < 0) stop("'effect' must be nonnegative.")
  if (baseline_stim_effect < 0) stop("'baseline_stim_effect' must be nonnegative.")
  if (bg_effect < 0) stop("'bg_effect' must be nonnegative.")
  if (length(components) != 8 || any(components < 0) || abs(sum(components) - 1) > 1e-8)
    stop("'components' must be 8 nonnegative proportions summing to one.")
  if (length(rng) %% 2 != 0) stop("'rng' must have even length (min, max pairs).")
  family <- normalise_family(prior)

  # ---- total cell counts (unchanged): columns nu1, ns1, nu0, ns0 ----
  rng_8 <- rep_len(rng, 8)
  mins <- rng_8[c(1, 3, 5, 7)]
  maxs <- rng_8[c(2, 4, 6, 8)]
  Ntot <- matrix(0, nrow = P, ncol = 4)
  for (d in 1:4) Ntot[, d] <- round(runif(P, min = mins[d], max = maxs[d]))
  colnames(Ntot) <- c("nu1", "ns1", "nu0", "ns0")

  mus <- condition_means(effect, bg_effect, baseline_stim_effect, baseline_background)
  MU1 <- mus[["u1"]]; MU0 <- mus[["u0"]]; MS0 <- mus[["s0"]]; MS1 <- mus[["s1"]]
  PHI <- if (length(phi) == 4) phi else rep(phi, 4)
  names(PHI) <- c("u1", "u0", "s0", "s1")     # same order as version 2: PHI[1..4]

  d <- function(n, mu, cond, phi_cond = cond) draw_prop(n, mu, cond, family, PHI[[phi_cond]])
  guard <- function(r, k) if (r > max_rounds)
    stop(sprintf("Rejection sampling did not finish for profile %s (family %s, effect %g)", PROFILE_LABELS[k], family, effect))

  n <- allocate_profiles(P, components)
  stuck <- c(0, 0, 0, 0)   # redraws of unstimulated proportions (CHANGE 9), by profile R1-R4
  PU1 <- PU0 <- PS1 <- PS0 <- numeric(0)
  add <- function(pu1, pu0, ps1, ps0) {
    PU1 <<- c(PU1, pu1); PU0 <<- c(PU0, pu0); PS1 <<- c(PS1, ps1); PS0 <<- c(PS0, ps0)
  }

  # ---- R1: all four different; DiD > 0 and s1 > u1 ----
  k <- 1
  if (n[k] > 0) {
    ps0 <- ps1 <- rep(0, n[k])
    pu1 <- d(n[k], MU1, "u1")
    pu0 <- d(n[k], MU0, "u0")
    r <- 0; stuck_sub <- rep(FALSE, n[k])
    while (any(bad <- (ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1))) {
      r <- r + 1; guard(r, k); m <- sum(bad)
      ps1[bad] <- d(m, MS1, "s1")
      ps0[bad] <- d(m, MS0, "s0")
      if (r > stuck_rounds) {               # [CHANGE 9] see header
        pu1[bad] <- d(m, MU1, "u1"); pu0[bad] <- d(m, MU0, "u0")
        stuck_sub[bad] <- TRUE
      }
    }
    stuck[k] <- sum(stuck_sub)
    if (effect == 0) ps1 <- ps0
    add(pu1, pu0, ps1, ps0)
  }

  # ---- R2: s0 = u0; s1 > u1 ----
  k <- 2
  if (n[k] > 0) {
    pu1 <- d(n[k], MU1, "u1")
    pu0 <- d(n[k], MU0, "u0")
    ps0 <- pu0
    ps1 <- d(n[k], MS1, "s1")
    r <- 0
    while (any(bad <- (ps1 - pu1 <= 0))) {
      r <- r + 1; guard(r, k); m <- sum(bad)
      ps1[bad] <- d(m, MS1, "s1")
      pu1[bad] <- d(m, MU1, "u1")
    }
    add(pu1, pu0, ps1, ps0)
  }

  # ---- R3: s1 = s0 (mean MS1); DiD > 0, s1 > u1, u0 > u1 ----
  k <- 3
  if (n[k] > 0) {
    ps0 <- ps1 <- d(n[k], MS1, "s1")
    pu0 <- d(n[k], MU0, "u0")
    pu1 <- d(n[k], MU1, "u1")
    r <- 0
    while (any(bad <- (ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | pu0 <= pu1))) {
      r <- r + 1; guard(r, k); m <- sum(bad)
      pu0[bad] <- d(m, MU0, "u0")
      pu1[bad] <- d(m, MU1, "u1")
    }
    add(pu1, pu0, ps1, ps0)
  }

  # ---- R4: u1 = u0; DiD > 0, s1 > u1, s1 > s0 ----
  k <- 4
  if (n[k] > 0) {
    pu0 <- pu1 <- d(n[k], MU0, "u0")
    ps1 <- d(n[k], MS1, "s1")
    ps0 <- d(n[k], MS0, "s0")
    r <- 0; stuck_sub <- rep(FALSE, n[k])
    while (any(bad <- (ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | ps1 <= ps0))) {
      r <- r + 1; guard(r, k); m <- sum(bad)
      ps0[bad] <- d(m, MS0, "s0")
      ps1[bad] <- d(m, MS1, "s1")
      if (r > stuck_rounds) {               # [CHANGE 9] see header
        pu0[bad] <- pu1[bad] <- d(m, MU0, "u0")
        stuck_sub[bad] <- TRUE
      }
    }
    stuck[k] <- sum(stuck_sub)
    add(pu1, pu0, ps1, ps0)
  }

  # ---- NR1: s0 = u0, s1 = u1 ----
  k <- 5
  if (n[k] > 0) {
    ps1 <- pu1 <- d(n[k], MU1, "u1")
    ps0 <- pu0 <- d(n[k], MU0, "u0")
    add(pu1, pu0, ps1, ps0)
  }

  # ---- NR2: s1 = u1; s0 >= u0 ----
  k <- 6
  if (n[k] > 0) {
    ps1 <- pu1 <- d(n[k], MU1, "u1")
    ps0 <- d(n[k], MS0, "s0")
    pu0 <- d(n[k], MU0, "u0")
    r <- 0
    while (any(bad <- (ps0 < pu0))) {
      r <- r + 1; guard(r, k); m <- sum(bad)
      ps0[bad] <- d(m, MS0, "s0")
      pu0[bad] <- d(m, MU0, "u0")
    }
    add(pu1, pu0, ps1, ps0)
  }

  # ---- NR3: all four equal ----
  k <- 7
  if (n[k] > 0) {
    ps1 <- ps0 <- pu1 <- pu0 <- d(n[k], MU0, "u0")
    add(pu1, pu0, ps1, ps0)
  }

  # ---- NSR: s1 = s0 (mean MS0), u1 = u0 ----
  k <- 8
  if (n[k] > 0) {
    ps0 <- ps1 <- d(n[k], MS0, "s0")
    pu0 <- pu1 <- d(n[k], MU0, "u0")
    add(pu1, pu0, ps1, ps0)
  }

  # ---- binomial counts (same order of draws as version 2) ----
  nu1 <- rbinom(P, Ntot[, "nu1"], PU1)
  ns1 <- rbinom(P, Ntot[, "ns1"], PS1)
  nu0 <- rbinom(P, Ntot[, "nu0"], PU0)
  ns0 <- rbinom(P, Ntot[, "ns0"], PS0)

  truth <- rep(PROFILE_LABELS, n)
  list(Ntot = Ntot, ns0 = ns0, ns1 = ns1, nu0 = nu0, nu1 = nu1, truth = truth,
       p = cbind(pu1 = PU1, ps1 = PS1, pu0 = PU0, ps0 = PS0),
       true_delta = (PS1 - PU1) - (PS0 - PU0),
       n_profile = setNames(n, PROFILE_LABELS),
       n_stuck_redraws = sum(stuck),
       family = family)
}

is_responder <- function(truth) as.integer(truth %in% RESPONDER_LABELS)

# Combine two simulated datasets (effect-heterogeneity study)
bind_sims <- function(a, b) {
  list(Ntot = rbind(a$Ntot, b$Ntot),
       ns1 = c(a$ns1, b$ns1), nu1 = c(a$nu1, b$nu1),
       ns0 = c(a$ns0, b$ns0), nu0 = c(a$nu0, b$nu0),
       truth = c(a$truth, b$truth),
       p = rbind(a$p, b$p),
       true_delta = c(a$true_delta, b$true_delta),
       n_profile = a$n_profile + b$n_profile,
       n_stuck_redraws = a$n_stuck_redraws + b$n_stuck_redraws,
       family = a$family)
}
