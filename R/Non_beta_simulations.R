# Non-beta generative mechanisms
# Isabella Lethbridge and Tayyeb Abrahams
# 30 June 2026

#-------------------------------------------------------------------------------
# Setup:

#-------------------------------------------------------------------------------
# Function to simulate from MIMOSA2 using a uniform prior
simulate_MIMOSA2_uniform <- function(
  effect = 5e-4,
  bg_effect = 0,
  baseline_stim_effect = 2.5e-4,
  baseline_background = 1e-4,
  phi = 5000,
  P = 100,
  rng = c(100000, 150000)
) {
  # Run 100 simulations, dropping any errors (i.e. simulations that error out)
  K <- 8
  n <- rep(0, K)
  while (sum(n) != P | any(P == 0)) {
    #proportion of vaccine specific responses
    pis <- (prop.table(runif(K)))
    R <- NULL
    D <- 4

    # Total observations simulated
    # from uniform with between 100,000 and 150,000 cells
    Ntot <- matrix(round(runif(P * D, rng[1], rng[2])), ncol = D, nrow = P)

    # Hyperprior configurations
    # No longer used for beta!!
    MU0 <- baseline_background
    MS0 <- baseline_stim_effect + MU0
    # Non Stimulated post vaccine
    MU1 <- MU0 + bg_effect
    # Mean of hyperprior (stimulated time 1) $\alpha/(\alpha+\beta)$
    MS1 <- MU1 + effect

    #Precision of hyperprior $\alpha+\beta = \phi$
    if (length(phi) != 4) {
      PHI <- rep(phi, 4)
    } else {
      PHI <- phi
    }

    # There are 8 model components we simulate from:
    # 1. all different
    # 2. s0=u0
    # 3. s1 = s0
    # 4. u1 = u0
    # 5. s0 = u0, s1 = u1
    # 6. s1 = u1
    # 7. s1 = u1 = s0 = u0
    # 8. s1 = s0, u1 = u0

    # number of responders / non-responders / non-specific responders
    n <- round(P * pis)
  }

  PS0 <- PS1 <- PU0 <- PU1 <- NULL

  # Simulate from component 1.
  k <- 1
  ps0 <- rep(0, n[k])
  ps1 <- rep(0, n[k])
  pu1 <- runif(n[k], 0, 1)
  pu0 <- runif(n[k], 0, 1)
  # pu1 = rbeta(n[k], MU1 * PHI[1], (1 - MU1) * PHI[1])
  # pu0 = rbeta(n[k], MU0 * PHI[2], (1-MU0) * PHI[2])

  while (any(ps1 <= pu1 | ps0 <= pu0 | ps1 - pu1 <= ps0 - pu0)) {
    bar <- ps1 <= pu1 | ps0 <= pu0 | ps1 - pu1 <= ps0 - pu0
    foo <- sum(bar)
    ps0[bar] <- runif(foo, 0, 1)
    ps1[bar] <- runif(foo, 0, 1)
    # ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
    # ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
  }

  if (effect == 0) {
    ps1 <- ps0
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 2.
  k <- 2
  pu1 <- runif(n[k], 0, 1)
  pu0 <- runif(n[k], 0, 1)
  ps0 <- pu0
  ps1 <- runif(n[k], 0, 1)
  # pu1 = rbeta(n[k], MU1 * PHI[1], (1 - MU1) * PHI[1])
  # pu0 = rbeta(n[k], MU0 * PHI[2], (1 - MU0) * PHI[2])
  # ps0 = pu0
  # ps1 = rbeta(n[k],MS1*PHI[4], (1-MS1)*PHI[4])

  while (any(ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1)) {
    bar <- ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1
    foo <- sum(bar)
    ps1[bar] <- runif(foo, 0, 1)
    pu1[bar] <- runif(foo, 0, 1)
    # ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
    # pu1[bar] = rbeta(foo, MU1 * PHI[1], (1 - MU1) * PHI[1])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 3.
  k <- 3
  ps0 <- ps1 <- runif(n[k], 0, 1)
  pu0 <- runif(n[k], 0, 1)
  pu1 <- runif(n[k], 0, 1)
  # ps0 = ps1 = rbeta(n[k],MS1*PHI[4],(1-MS1)*PHI[4])
  # pu0 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  # pu1 = rbeta(n[k],MU1*PHI[1],(1-MU1)*PHI[1])

  while (any(ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | pu0 <= pu1)) {
    bar <- ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | pu0 <= pu1
    foo <- sum(bar)
    pu0[bar] <- runif(foo, 0, 1)
    pu1[bar] <- runif(foo, 0, 1)
    # pu0[bar] = rbeta(foo, MU0 * PHI[2], (1 - MU0) * PHI[2])
    # pu1[bar] = rbeta(foo, MU1 * PHI[1], (1 - MU1) * PHI[1])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 4.
  k <- 4
  pu0 <- pu1 <- runif(n[k], 0, 1)
  ps1 <- runif(n[k], 0, 1)
  ps0 <- runif(n[k], 0, 1)
  # pu0 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  # ps1 = rbeta(n[k],MS1*PHI[4],(1-MS1)*PHI[4])
  # ps0 = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])

  # iter = 0
  while (any(ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | ps1 <= ps0)) {
    bar <- ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | ps1 <= ps0
    foo <- sum(bar)
    ps0[bar] <- runif(foo, 0, 1)
    ps1[bar] <- runif(foo, 0, 1)
    # ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
    # ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 5.
  k <- 5
  ps1 <- pu1 <- runif(n[k], 0, 1)
  ps0 <- pu0 <- runif(n[k], 0, 1)
  # ps1 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  # ps0 = pu0 = rbeta(n[k],MU1*PHI[1],(1-MU1)*PHI[1])

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 6.
  k <- 6
  ps1 <- pu1 <- runif(n[k], 0, 1)
  ps0 <- runif(n[k], 0, 1)
  pu0 <- runif(n[k], 0, 1)
  # ps1 = pu1 = rbeta(n[k],MU1*PHI[2],(1-MU1)*PHI[2])
  # ps0  = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])
  # pu0  = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  #s0 should be greater than u0

  while (any(ps0 < pu0)) {
    bar <- ps0 < pu0
    foo <- sum(bar)
    ps0[bar] <- runif(foo, 0, 1)
    pu0[bar] <- runif(foo, 0, 1)
    # ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
    # pu0[bar] = rbeta(foo, MU0 * PHI[2], (1 - MU0) * PHI[2])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 7.
  k <- 7
  ps1 <- ps0 <- pu1 <- pu0 <- runif(n[k], 0, 1)
  # ps1=ps0=pu1=pu0 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 8.
  k <- 8
  if (n[k] > 0) {
    ps0 <- ps1 <- runif(n[k], 0, 1)
    pu0 <- pu1 <- runif(n[k], 0, 1)
    # ps0 = ps1 = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])
    # pu0 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])

    PU1 <- c(PU1, pu1)
    PU0 <- c(PU0, pu0)
    PS1 <- c(PS1, ps1)
    PS0 <- c(PS0, ps0)
  }

  # Simulate count observations from binomial.
  colnames(Ntot) <- c("nu1", "ns1", "nu0", "ns0")
  nu1 <- rbinom(P, Ntot[, "nu1"], PU1)
  ns1 <- rbinom(P, Ntot[, "ns1"], PS1)
  nu0 <- rbinom(P, Ntot[, "nu0"], PU0)
  ns0 <- rbinom(P, Ntot[, "ns0"], PS0)

  # Empirical estimates of proportions
  pu1_hat <- prop.table(cbind(Ntot[, "nu1"], nu1), 1)[, 2]
  ps1_hat <- prop.table(cbind(Ntot[, "ns1"], ns1), 1)[, 2]
  pu0_hat <- prop.table(cbind(Ntot[, "nu0"], nu0), 1)[, 2]
  ps0_hat <- prop.table(cbind(Ntot[, "ns0"], ns0), 1)[, 2]

  # True response categories
  truth <- rep(
    c("R1", "R2", "R3", "R4", "NR1", "NR2", "NR3", "NSR"),
    n
  )
  l <- list(
    "Ntot" = Ntot,
    "ns0" = ns0,
    "ns1" = ns1,
    "nu0" = nu0,
    "nu1" = nu1,
    "truth" = truth
  )
  names(l) <- c("Ntot", "ns0", "ns1", "nu0", "nu1", "truth")
  return(l)
}

#-------------------------------------------------------------------------------
# Function to simulate from MIMOSA2 using a logit Normal prior
simulate_MIMOSA2_logitnorm <- function(
  effect = 5e-4,
  bg_effect = 0,
  baseline_stim_effect = 2.5e-4,
  baseline_background = 1e-4,
  phi = 5000,
  P = 100,
  rng = c(100000, 150000)
) {
  # Function to generate logit-normal proportions:
  rlogitnorm <- function(n, mean, phi) {
    if (n <= 0) {
      return(numeric(0))
    }

    sdlogit <- 1 / sqrt(phi)
    meanlogit <- qlogis(mean)

    draws <- plogis(rnorm(n, mean = meanlogit, sd = sdlogit))
    return(draws)
  }

  # Run 100 simulations, dropping any errors (i.e. simulations that error out)
  K <- 8
  n <- rep(0, K)
  while (sum(n) != P | any(P == 0)) {
    #proportion of vaccine specific responses
    pis <- (prop.table(runif(K)))
    R <- NULL
    D <- 4

    # Total observations simulated
    # from uniform with between 100,000 and 150,000 cells
    Ntot <- matrix(round(runif(P * D, rng[1], rng[2])), ncol = D, nrow = P)

    # Hyperprior configurations
    # No longer used for beta!!
    MU0 <- baseline_background
    MS0 <- baseline_stim_effect + MU0
    # Non Stimulated post vaccine
    MU1 <- MU0 + bg_effect
    # Mean of hyperprior (stimulated time 1) $\alpha/(\alpha+\beta)$
    MS1 <- MU1 + effect

    #Precision of hyperprior $\alpha+\beta = \phi$
    if (length(phi) != 4) {
      PHI <- rep(phi, 4)
    } else {
      PHI <- phi
    }

    # There are 8 model components we simulate from:
    # 1. all different
    # 2. s0=u0
    # 3. s1 = s0
    # 4. u1 = u0
    # 5. s0 = u0, s1 = u1
    # 6. s1 = u1
    # 7. s1 = u1 = s0 = u0
    # 8. s1 = s0, u1 = u0

    # number of responders / non-responders / non-specific responders
    n <- round(P * pis)
  }

  PS0 <- PS1 <- PU0 <- PU1 <- NULL

  # Simulate from component 1.
  k <- 1
  ps0 <- rep(0, n[k])
  ps1 <- rep(0, n[k])
  pu1 <- rlogitnorm(n[k], MU1, PHI[1])
  pu0 <- rlogitnorm(n[k], MU0, PHI[2])
  # pu1 = rbeta(n[k], MU1 * PHI[1], (1 - MU1) * PHI[1])
  # pu0 = rbeta(n[k], MU0 * PHI[2], (1-MU0) * PHI[2])

  while (any(ps1 <= pu1 | ps0 <= pu0 | ps1 - pu1 <= ps0 - pu0)) {
    bar <- ps1 <= pu1 | ps0 <= pu0 | ps1 - pu1 <= ps0 - pu0
    foo <- sum(bar)
    ps0[bar] <- rlogitnorm(foo, MS0, PHI[3])
    ps1[bar] <- rlogitnorm(foo, MS1, PHI[4])
    # ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
    # ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
  }

  if (effect == 0) {
    ps1 <- ps0
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 2.
  k <- 2
  pu1 <- rlogitnorm(n[k], MU1, PHI[1])
  pu0 <- rlogitnorm(n[k], MU0, PHI[2])
  ps0 <- pu0
  ps1 <- rlogitnorm(n[k], MS1, PHI[4])
  # pu1 = rbeta(n[k], MU1 * PHI[1], (1 - MU1) * PHI[1])
  # pu0 = rbeta(n[k], MU0 * PHI[2], (1 - MU0) * PHI[2])
  # ps0 = pu0
  # ps1 = rbeta(n[k],MS1*PHI[4], (1-MS1)*PHI[4])

  while (any(ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1)) {
    bar <- ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1
    foo <- sum(bar)
    ps1[bar] <- rlogitnorm(foo, MS1, PHI[4])
    pu1[bar] <- rlogitnorm(foo, MU1, PHI[1])
    # ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
    # pu1[bar] = rbeta(foo, MU1 * PHI[1], (1 - MU1) * PHI[1])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 3.
  k <- 3
  ps0 <- ps1 <- rlogitnorm(n[k], MS1, PHI[4])
  pu0 <- rlogitnorm(n[k], MU0, PHI[2])
  pu1 <- rlogitnorm(n[k], MU1, PHI[1])
  # ps0 = ps1 = rbeta(n[k],MS1*PHI[4],(1-MS1)*PHI[4])
  # pu0 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  # pu1 = rbeta(n[k],MU1*PHI[1],(1-MU1)*PHI[1])

  while (any(ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | pu0 <= pu1)) {
    bar <- ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | pu0 <= pu1
    foo <- sum(bar)
    pu0[bar] <- rlogitnorm(foo, MU0, PHI[2])
    pu1[bar] <- rlogitnorm(foo, MU1, PHI[1])
    # pu0[bar] = rbeta(foo, MU0 * PHI[2], (1 - MU0) * PHI[2])
    # pu1[bar] = rbeta(foo, MU1 * PHI[1], (1 - MU1) * PHI[1])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 4.
  k <- 4
  pu0 <- pu1 <- rlogitnorm(n[k], MU0, PHI[2])
  ps1 <- rlogitnorm(n[k], MS1, PHI[4])
  ps0 <- rlogitnorm(n[k], MS0, PHI[3])
  # pu0 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  # ps1 = rbeta(n[k],MS1*PHI[4],(1-MS1)*PHI[4])
  # ps0 = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])

  # iter = 0
  while (any(ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | ps1 <= ps0)) {
    bar <- ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | ps1 <= ps0
    foo <- sum(bar)
    ps0[bar] <- rlogitnorm(foo, MS0, PHI[3])
    ps1[bar] <- rlogitnorm(foo, MS1, PHI[4])
    # ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
    # ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 5.
  k <- 5
  ps1 <- pu1 <- rlogitnorm(n[k], MU0, PHI[2])
  ps0 <- pu0 <- rlogitnorm(n[k], MU1, PHI[1])
  # ps1 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  # ps0 = pu0 = rbeta(n[k],MU1*PHI[1],(1-MU1)*PHI[1])

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 6.
  k <- 6
  ps1 <- pu1 <- rlogitnorm(n[k], MU1, PHI[2])
  ps0 <- rlogitnorm(n[k], MS0, PHI[3])
  pu0 <- rlogitnorm(n[k], MU0, PHI[2])
  # ps1 = pu1 = rbeta(n[k],MU1*PHI[2],(1-MU1)*PHI[2])
  # ps0  = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])
  # pu0  = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  #s0 should be greater than u0

  while (any(ps0 < pu0)) {
    bar <- ps0 < pu0
    foo <- sum(bar)
    ps0[bar] <- rlogitnorm(foo, MS0, PHI[3])
    pu0[bar] <- rlogitnorm(foo, MU0, PHI[2])
    # ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
    # pu0[bar] = rbeta(foo, MU0 * PHI[2], (1 - MU0) * PHI[2])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 7.
  k <- 7
  ps1 <- ps0 <- pu1 <- pu0 <- rlogitnorm(n[k], MU0, PHI[2])
  # ps1=ps0=pu1=pu0 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 8.
  k <- 8
  if (n[k] > 0) {
    ps0 <- ps1 <- rlogitnorm(n[k], MS0, PHI[3])
    pu0 <- pu1 <- rlogitnorm(n[k], MU0, PHI[2])
    # ps0 = ps1 = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])
    # pu0 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])

    PU1 <- c(PU1, pu1)
    PU0 <- c(PU0, pu0)
    PS1 <- c(PS1, ps1)
    PS0 <- c(PS0, ps0)
  }

  # Simulate count observations from binomial.
  colnames(Ntot) <- c("nu1", "ns1", "nu0", "ns0")
  nu1 <- rbinom(P, Ntot[, "nu1"], PU1)
  ns1 <- rbinom(P, Ntot[, "ns1"], PS1)
  nu0 <- rbinom(P, Ntot[, "nu0"], PU0)
  ns0 <- rbinom(P, Ntot[, "ns0"], PS0)

  # Empirical estimates of proportions
  pu1_hat <- prop.table(cbind(Ntot[, "nu1"], nu1), 1)[, 2]
  ps1_hat <- prop.table(cbind(Ntot[, "ns1"], ns1), 1)[, 2]
  pu0_hat <- prop.table(cbind(Ntot[, "nu0"], nu0), 1)[, 2]
  ps0_hat <- prop.table(cbind(Ntot[, "ns0"], ns0), 1)[, 2]

  # True response categories
  truth <- rep(
    c("R1", "R2", "R3", "R4", "NR1", "NR2", "NR3", "NSR"),
    n
  )
  l <- list(
    "Ntot" = Ntot,
    "ns0" = ns0,
    "ns1" = ns1,
    "nu0" = nu0,
    "nu1" = nu1,
    "truth" = truth
  )
  names(l) <- c("Ntot", "ns0", "ns1", "nu0", "nu1", "truth")
  return(l)
}

#-------------------------------------------------------------------------------
# Function to simulate from MIMOSA2 using X/(1+X), X~exp()
simulate_MIMOSA2_rexp <- function(
  effect = 5e-4,
  bg_effect = 0,
  baseline_stim_effect = 2.5e-4,
  baseline_background = 1e-4,
  phi = 5000,
  P = 100,
  rng = c(100000, 150000)
) {
  # Function to generate proportions from X/(1+X), X~exponential
  rexp_prop <- function(n, mean) {
    if (n <= 0) {
      return(numeric(0))
    }

    lambda <- (1 - mean) / mean

    x <- rexp(n, rate = lambda)
    return(x / (1 + x))
  }

  # Run 100 simulations, dropping any errors (i.e. simulations that error out)
  K <- 8
  n <- rep(0, K)
  while (sum(n) != P | any(P == 0)) {
    #proportion of vaccine specific responses
    pis <- (prop.table(runif(K)))
    R <- NULL
    D <- 4

    # Total observations simulated
    # from uniform with between 100,000 and 150,000 cells
    Ntot <- matrix(round(runif(P * D, rng[1], rng[2])), ncol = D, nrow = P)

    # Hyperprior configurations
    MU0 <- baseline_background
    MS0 <- baseline_stim_effect + MU0
    # Non Stimulated post vaccine
    MU1 <- MU0 + bg_effect
    # Mean of hyperprior (stimulated time 1) $\alpha/(\alpha+\beta)$
    MS1 <- MU1 + effect

    #Precision of hyperprior $\alpha+\beta = \phi$
    if (length(phi) != 4) {
      PHI <- rep(phi, 4)
    } else {
      PHI <- phi
    }

    # There are 8 model components we simulate from:
    # 1. all different
    # 2. s0=u0
    # 3. s1 = s0
    # 4. u1 = u0
    # 5. s0 = u0, s1 = u1
    # 6. s1 = u1
    # 7. s1 = u1 = s0 = u0
    # 8. s1 = s0, u1 = u0

    # number of responders / non-responders / non-specific responders
    n <- round(P * pis)
  }

  PS0 <- PS1 <- PU0 <- PU1 <- NULL

  # Simulate from component 1.
  k <- 1
  ps0 <- rep(0, n[k])
  ps1 <- rep(0, n[k])
  pu1 <- rexp_prop(n[k], MU1)
  pu0 <- rexp_prop(n[k], MU0)
  # pu1 = rbeta(n[k], MU1 * PHI[1], (1 - MU1) * PHI[1])
  # pu0 = rbeta(n[k], MU0 * PHI[2], (1-MU0) * PHI[2])

  while (any(ps1 <= pu1 | ps0 <= pu0 | ps1 - pu1 <= ps0 - pu0)) {
    bar <- ps1 <= pu1 | ps0 <= pu0 | ps1 - pu1 <= ps0 - pu0
    foo <- sum(bar)
    ps0[bar] <- rexp_prop(foo, MS0)
    ps1[bar] <- rexp_prop(foo, MS1)
    # ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
    # ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
  }

  if (effect == 0) {
    ps1 <- ps0
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 2.
  k <- 2
  pu1 <- rexp_prop(n[k], MU1)
  pu0 <- rexp_prop(n[k], MU0)
  ps0 <- pu0
  ps1 <- rexp_prop(n[k], MS1)
  # pu1 = rbeta(n[k], MU1 * PHI[1], (1 - MU1) * PHI[1])
  # pu0 = rbeta(n[k], MU0 * PHI[2], (1 - MU0) * PHI[2])
  # ps0 = pu0
  # ps1 = rbeta(n[k],MS1*PHI[4], (1-MS1)*PHI[4])

  while (any(ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1)) {
    bar <- ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1
    foo <- sum(bar)
    ps1[bar] <- rexp_prop(foo, MS1)
    pu1[bar] <- rexp_prop(foo, MU1)
    # ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
    # pu1[bar] = rbeta(foo, MU1 * PHI[1], (1 - MU1) * PHI[1])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 3.
  k <- 3
  ps0 <- ps1 <- rexp_prop(n[k], MS1)
  pu0 <- rexp_prop(n[k], MU0)
  pu1 <- rexp_prop(n[k], MU1)
  # ps0 = ps1 = rbeta(n[k],MS1*PHI[4],(1-MS1)*PHI[4])
  # pu0 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  # pu1 = rbeta(n[k],MU1*PHI[1],(1-MU1)*PHI[1])

  while (any(ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | pu0 <= pu1)) {
    bar <- ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | pu0 <= pu1
    foo <- sum(bar)
    pu0[bar] <- rexp_prop(foo, MU0)
    pu1[bar] <- rexp_prop(foo, MU1)
    # pu0[bar] = rbeta(foo, MU0 * PHI[2], (1 - MU0) * PHI[2])
    # pu1[bar] = rbeta(foo, MU1 * PHI[1], (1 - MU1) * PHI[1])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 4.
  k <- 4
  pu0 <- pu1 <- rexp_prop(n[k], MU0)
  ps1 <- rexp_prop(n[k], MS1)
  ps0 <- rexp_prop(n[k], MS0)
  # pu0 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  # ps1 = rbeta(n[k],MS1*PHI[4],(1-MS1)*PHI[4])
  # ps0 = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])

  # iter = 0
  while (any(ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | ps1 <= ps0)) {
    bar <- ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1 | ps1 <= ps0
    foo <- sum(bar)
    ps0[bar] <- rexp_prop(foo, MS0)
    ps1[bar] <- rexp_prop(foo, MS1)
    # ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
    # ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 5.
  k <- 5
  ps1 <- pu1 <- rexp_prop(n[k], MU0)
  ps0 <- pu0 <- rexp_prop(n[k], MU1)
  # ps1 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  # ps0 = pu0 = rbeta(n[k],MU1*PHI[1],(1-MU1)*PHI[1])

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 6.
  k <- 6
  ps1 <- pu1 <- rexp_prop(n[k], MU1)
  ps0 <- rexp_prop(n[k], MS0)
  pu0 <- rexp_prop(n[k], MU0)
  # ps1 = pu1 = rbeta(n[k],MU1*PHI[2],(1-MU1)*PHI[2])
  # ps0  = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])
  # pu0  = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  #s0 should be greater than u0

  while (any(ps0 < pu0)) {
    bar <- ps0 < pu0
    foo <- sum(bar)
    ps0[bar] <- rexp_prop(foo, MS0)
    pu0[bar] <- rexp_prop(foo, MU0)
    # ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
    # pu0[bar] = rbeta(foo, MU0 * PHI[2], (1 - MU0) * PHI[2])
  }

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 7.
  k <- 7
  ps1 <- ps0 <- pu1 <- pu0 <- rexp_prop(n[k], MU0)
  # ps1=ps0=pu1=pu0 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])

  PU1 <- c(PU1, pu1)
  PU0 <- c(PU0, pu0)
  PS1 <- c(PS1, ps1)
  PS0 <- c(PS0, ps0)

  # Simulate from component 8.
  k <- 8
  if (n[k] > 0) {
    ps0 <- ps1 <- rexp_prop(n[k], MS0)
    pu0 <- pu1 <- rexp_prop(n[k], MU0)
    # ps0 = ps1 = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])
    # pu0 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])

    PU1 <- c(PU1, pu1)
    PU0 <- c(PU0, pu0)
    PS1 <- c(PS1, ps1)
    PS0 <- c(PS0, ps0)
  }

  # Simulate count observations from binomial.
  colnames(Ntot) <- c("nu1", "ns1", "nu0", "ns0")
  nu1 <- rbinom(P, Ntot[, "nu1"], PU1)
  ns1 <- rbinom(P, Ntot[, "ns1"], PS1)
  nu0 <- rbinom(P, Ntot[, "nu0"], PU0)
  ns0 <- rbinom(P, Ntot[, "ns0"], PS0)

  # Empirical estimates of proportions
  pu1_hat <- prop.table(cbind(Ntot[, "nu1"], nu1), 1)[, 2]
  ps1_hat <- prop.table(cbind(Ntot[, "ns1"], ns1), 1)[, 2]
  pu0_hat <- prop.table(cbind(Ntot[, "nu0"], nu0), 1)[, 2]
  ps0_hat <- prop.table(cbind(Ntot[, "ns0"], ns0), 1)[, 2]

  # True response categories
  truth <- rep(
    c("R1", "R2", "R3", "R4", "NR1", "NR2", "NR3", "NSR"),
    n
  )
  l <- list(
    "Ntot" = Ntot,
    "ns0" = ns0,
    "ns1" = ns1,
    "nu0" = nu0,
    "nu1" = nu1,
    "truth" = truth
  )
  names(l) <- c("Ntot", "ns0", "ns1", "nu0", "nu1", "truth")
  return(l)
}
