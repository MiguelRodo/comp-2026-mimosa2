my_libs <- "/scratch/abrmoe030/R_libs"
.libPaths(c(my_libs, .libPaths()))

library(MIMOSA2)
library(ggplot2)

simulate_MIMOSA2_alt_prior = function(effect = 5e-4, 
                                      bg_effect = 0,
                                      baseline_stim_effect=2.5e-4,
                                      baseline_background=1e-4,
                                      phi = 5000,
                                      P = 100,
                                      rng = c(100000,150000), prior="beta",
                                      components=rep(1/8,8)) {
  if (effect < 0) stop("'effect' must be nonnegative.")
  if (baseline_stim_effect < 0) stop("'baseline_stim_effect' must be nonnegative.")
  if (bg_effect < 0) stop("'bg_effect' must be nonnegative.")
  if (sum(components)!=1) stop("Component proportions must sum to one")
  
  K = 8
  n = rep(0, K)
  
  pis = components
  R = NULL
  D = 4
  
  Ntot = matrix(round(runif(P * D, rng[1], rng[2])), ncol = D, nrow = P)
  MU0 = baseline_background
  MS0 = baseline_stim_effect + MU0
  MU1 = MU0 + bg_effect
  MS1 = MU1 + effect
  
  if(length(phi) != 4) {
    PHI = rep(phi, 4)
  } else {
    PHI = phi
  }
  
n = round(P * pis)
  n[8] = max(P - sum(n[1:7]), 0)
  
  # Shave off counts one by one from the largest components until the sum equals P
  while (sum(n) > P) {
    idx = which.max(n)
    n[idx] = n[idx] - 1
  }
  while (sum(n) < P) {
    idx = which.min(n)
    n[idx] = n[idx] + 1
  }
  
  PS0=PS1=PU0=PU1=NULL
  is_beta = prior %in% c("beta","b","B","Beta")
  is_unif = prior %in% c("Uniform","uniform","u","U")
  is_logit_normal = prior %in% c("Logit_normal","logit_normal","Logit normal",
                                 "logit normal","Logit norm","logit norm",
                                 "Logit_Normal","Logit Normal","Logit Norm",
                                 "ln","Ln","LN")
  is_odds_exponential = prior %in% c("Odds_exponential","Odds exponential",
                                     "odds_exponential","odds exponential",
                                     "Odds Exponential","Odds_Exponential",
                                     "Odds exp","Odds Exp","Odds_exp","Odds_Exp",
                                     "odds exp","odds_exp",
                                     "Oe","oe","OE")
  is_exponential_gamma = prior %in% c("Exponential gamma","Exponential_gamma",
                                      "exponential gamma","exponential_gamma",
                                      "Exponential Gamma","Exponential_Gamma",
                                      "Eg","EG","eg")
  is_odds_gamma = prior %in% c("Odds_gamma", "Odds gamma", "odds_gamma", "odds gamma", 
                               "Og", "OG", "og")
  is_unit_lognormal = prior %in% c("Unit_lognormal", "Unit lognormal", "unit_lognormal", 
                                   "unit lognormal", "Uln", "ULN", "uln")
  is_simplex = prior %in% c("Simplex", "simplex", "S", "s", "Sx", "SX", "sx")
  
  # ==========================================
  # Samplers Suite
  # ==========================================
  rlogitnorm = function(num_draws, target_mean, component_phi) {
    if (num_draws <= 0) return(numeric(0))
    sdlogit   = 1 / sqrt(component_phi)
    meanlogit = qlogis(target_mean)
    draws     = plogis(rnorm(num_draws, mean = meanlogit, sd = sdlogit))
    return(draws)
  }
  
  rodds_exp = function(num_draws, target_mean) {
    obs = rexp(num_draws, rate = 1/target_mean)
    return(obs / (1 + obs))
  }
  
  rgamma_exp = function(num_draws, target_mean, gamma_dispersion){
    gamma_mean <- gamma_dispersion*(target_mean^(-1/gamma_dispersion)-1)
    outs <- exp(-rgamma(num_draws, gamma_dispersion, gamma_dispersion/gamma_mean))
    return(outs)
  }
  
  rodds_gamma = function(num_draws, target_mean, component_phi) {
    if (num_draws <= 0) return(numeric(0))
    lambda = target_mean / (1 - target_mean)
    x = rgamma(num_draws, shape = component_phi, rate = component_phi / lambda)
    return(x / (1 + x))
  }
  
  runit_lognormal = function(num_draws, target_mean, component_phi) {
    if (num_draws <= 0) return(numeric(0))
    sigma2 = 1 / component_phi
    theta = -log(-log(target_mean) / (1 + 0.5 * sigma2))
    norm_draws = rnorm(num_draws, mean = theta, sd = sqrt(sigma2))
    return(exp(-exp(norm_draws)))
  }
  
  rsimplex = function(num_draws, target_mean, component_phi) {
    if (num_draws <= 0) return(numeric(0))
    mu_ig = target_mean / (1 - target_mean)
    lambda_ig = component_phi * (target_mean^2) * ((1 - target_mean)^2)
    
    v = rnorm(num_draws)^2
    x = mu_ig + (mu_ig^2 * v) / (2 * lambda_ig) - 
      (mu_ig / (2 * lambda_ig)) * sqrt(4 * mu_ig * lambda_ig * v + mu_ig^2 * v^2)
    z = runif(num_draws)
    indices = z > (mu_ig / (mu_ig + x))
    x[indices] = (mu_ig^2) / x[indices]
    
    return(x / (1 + x))
  }
  
  # Calculate uniform hyperprior ranges if selected
  if (is_unif) {
    mus <- c(MU1, MU0, MS0, MS1)
    names(mus) <- c("MU1", "MU0", "MS0", "MS1")
    phi_min <- 3 / pmin(mus, 1 - mus)^2
    bad <- PHI < phi_min
    
    if (any(bad)) {
      msg <- paste(sprintf("%s: supplied = %.0f, minimum = %.0f", names(mus)[bad], PHI[bad], ceiling(phi_min[bad])), collapse = "\n")
      stop(paste("Precisions too small for a Uniform prior.", "Required minimum precisions are:", msg, sep = "\n"))
    }
    radius <- sqrt(3/PHI)
  }
  
  # ==========================================
  # Component 1: All different
  # ==========================================
  k = 1
  if(n[k] > 0) {
    ps0 = rep(0, n[k])
    ps1 = rep(0, n[k])
    
    if (is_beta) {
      pu1 = rbeta(n[k], MU1 * PHI[1], (1 - MU1) * PHI[1])
      pu0 = rbeta(n[k], MU0 * PHI[2], (1 - MU0) * PHI[2])
    } else if (is_unif) {
      pu1 = runif(n[k], min = MU1-radius[1], max = MU1+radius[1])
      pu0 = runif(n[k], min = MU0-radius[2], max = MU0+radius[2])
    } else if (is_logit_normal) {
      pu1 = rlogitnorm(n[k], MU1, PHI[1])
      pu0 = rlogitnorm(n[k], MU0, PHI[2])
    } else if (is_odds_exponential) {
      pu1 = rodds_exp(n[k], MU1)
      pu0 = rodds_exp(n[k], MU0)
    } else if (is_exponential_gamma) {
      pu1 = rgamma_exp(n[k], MU1, PHI[1])
      pu0 = rgamma_exp(n[k], MU0, PHI[2])
    } else if (is_odds_gamma) {
      pu1 = rodds_gamma(n[k], MU1, PHI[1])
      pu0 = rodds_gamma(n[k], MU0, PHI[2])
    } else if (is_unit_lognormal) {
      pu1 = runit_lognormal(n[k], MU1, PHI[1])
      pu0 = runit_lognormal(n[k], MU0, PHI[2])
    } else if (is_simplex) {
      pu1 = rsimplex(n[k], MU1, PHI[1])
      pu0 = rsimplex(n[k], MU0, PHI[2])
    }
    
    while (any(ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1)) {
      bar <- ps1 - pu1 <= ps0 - pu0 | ps1 <= pu1
      foo <- sum(bar)
      if (is_beta) {
        ps1[bar] <- rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
        ps0[bar] <- rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
      } else if (is_unif) {
        ps0[bar] <- runif(foo, min = MS0 - radius[3], max = MS0 + radius[3])
        ps1[bar] <- runif(foo, min = MS1 - radius[4], max = MS1 + radius[4])
      } else if (is_logit_normal) {
        ps1[bar] <- rlogitnorm(foo, MS1, PHI[4])
        ps0[bar] <- rlogitnorm(foo, MS0, PHI[3])
      } else if (is_odds_exponential) {
        ps1[bar] <- rodds_exp(foo, MS1)
        ps0[bar] <- rodds_exp(foo, MS0)
      } else if (is_exponential_gamma) {
        ps1[bar] <- rgamma_exp(foo, MS1, PHI[4])
        ps0[bar] <- rgamma_exp(foo, MS0, PHI[3])
      } else if (is_odds_gamma) {
        ps1[bar] <- rodds_gamma(foo, MS1, PHI[4])
        ps0[bar] <- rodds_gamma(foo, MS0, PHI[3])
      } else if (is_unit_lognormal) {
        ps1[bar] <- runit_lognormal(foo, MS1, PHI[4])
        ps0[bar] <- runit_lognormal(foo, MS0, PHI[3])
      } else if (is_simplex) {
        ps1[bar] <- rsimplex(foo, MS1, PHI[4])
        ps0[bar] <- rsimplex(foo, MS0, PHI[3])
      }
    }
    if (effect == 0) { ps1 = ps0 }
    PU1=c(PU1,pu1); PU0=c(PU0,pu0); PS1=c(PS1,ps1); PS0=c(PS0,ps0)
  }
  
  # ==========================================
  # Component 2: s0 = u0
  # ==========================================
  k = 2
  if(n[k] > 0) {
    if (is_beta) {
      pu1 = rbeta(n[k], MU1 * PHI[1], (1 - MU1) * PHI[1])
      pu0 = rbeta(n[k], MU0 * PHI[2], (1 - MU0) * PHI[2])
      ps0 = pu0
      ps1 = rbeta(n[k], MS1 * PHI[4], (1 - MS1) * PHI[4])
    } else if (is_unif) {
      pu1 = runif(n[k], min = MU1-radius[1], max = MU1+radius[1])
      pu0 = runif(n[k], min = MU0-radius[2], max = MU0+radius[2])
      ps0 = pu0
      ps1 = runif(n[k], min = MS1-radius[4], max = MS1+radius[4])
    } else if (is_logit_normal) {
      pu1 = rlogitnorm(n[k], MU1, PHI[1])
      pu0 = rlogitnorm(n[k], MU0, PHI[2])
      ps0 = pu0
      ps1 = rlogitnorm(n[k], MS1, PHI[4])
    } else if (is_odds_exponential) {
      pu1 = rodds_exp(n[k], MU1)
      pu0 = rodds_exp(n[k], MU0)
      ps0 = pu0
      ps1 = rodds_exp(n[k], MS1)
    } else if (is_exponential_gamma) {
      pu1 = rgamma_exp(n[k], MU1, PHI[1])
      pu0 = rgamma_exp(n[k], MU0, PHI[2])
      ps0 = pu0
      ps1 = rgamma_exp(n[k], MS1, PHI[4])
    } else if (is_odds_gamma) {
      pu1 = rodds_gamma(n[k], MU1, PHI[1])
      pu0 = rodds_gamma(n[k], MU0, PHI[2])
      ps0 = pu0
      ps1 = rodds_gamma(n[k], MS1, PHI[4])
    } else if (is_unit_lognormal) {
      pu1 = runit_lognormal(n[k], MU1, PHI[1])
      pu0 = runit_lognormal(n[k], MU0, PHI[2])
      ps0 = pu0
      ps1 = runit_lognormal(n[k], MS1, PHI[4])
    } else if (is_simplex) {
      pu1 = rsimplex(n[k], MU1, PHI[1])
      pu0 = rsimplex(n[k], MU0, PHI[2])
      ps0 = pu0
      ps1 = rsimplex(n[k], MS1, PHI[4])
    }
    
    while(any(ps1 - pu1 <= 0)) { 
      bar = ps1 - pu1 <= 0
      foo = sum(bar)
      if (is_beta) {
        ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
        pu1[bar] = rbeta(foo, MU1 * PHI[1], (1 - MU1) * PHI[1])
      } else if (is_unif) {
        pu1[bar] = runif(foo, min = MU1-radius[1], max = MU1+radius[1])
        ps1[bar] = runif(foo, min = MS1-radius[4], max = MS1+radius[4])
      } else if (is_logit_normal) {
        pu1[bar] = rlogitnorm(foo, MU1, PHI[1])
        ps1[bar] = rlogitnorm(foo, MS1, PHI[4])
      } else if (is_odds_exponential) {
        pu1[bar] = rodds_exp(foo, MU1)
        ps1[bar] = rodds_exp(foo, MS1)
      } else if (is_exponential_gamma) {
        pu1[bar] = rgamma_exp(foo, MU1, PHI[1])
        ps1[bar] = rgamma_exp(foo, MS1, PHI[4])
      } else if (is_odds_gamma) {
        pu1[bar] = rodds_gamma(foo, MU1, PHI[1])
        ps1[bar] = rodds_gamma(foo, MS1, PHI[4])
      } else if (is_unit_lognormal) {
        pu1[bar] = runit_lognormal(foo, MU1, PHI[1])
        ps1[bar] = runit_lognormal(foo, MS1, PHI[4])
      } else if (is_simplex) {
        pu1[bar] = rsimplex(foo, MU1, PHI[1])
        ps1[bar] = rsimplex(foo, MS1, PHI[4])
      }
    }
    PU1=c(PU1,pu1); PU0=c(PU0,pu0); PS1=c(PS1,ps1); PS0=c(PS0,ps0)
  }
  
  # ==========================================
  # Component 3: s1 = s0
  # ==========================================
  k = 3
  if(n[k] > 0) {
    if (is_beta) {
      ps0 = ps1 = rbeta(n[k], MS1*PHI[4], (1-MS1)*PHI[4])
      pu0 = rbeta(n[k], MU0*PHI[2], (1-MU0)*PHI[2])
      pu1 = rbeta(n[k], MU1*PHI[1], (1-MU1)*PHI[1])
    } else if (is_unif) {
      ps0 = ps1 = runif(n[k], min = MS1-radius[4], max = MS1+radius[4])
      pu0 = runif(n[k], min = MU0-radius[2], max = MU0+radius[2])
      pu1 = runif(n[k], min = MU1-radius[1], max = MU1+radius[1])
    } else if (is_logit_normal) {
      ps0 = ps1 = rlogitnorm(n[k], MS1, PHI[4])
      pu0 = rlogitnorm(n[k], MU0, PHI[2])
      pu1 = rlogitnorm(n[k], MU1, PHI[1])
    } else if (is_odds_exponential) {
      ps0 = ps1 = rodds_exp(n[k], MS1)
      pu0 = rodds_exp(n[k], MU0)
      pu1 = rodds_exp(n[k], MU1)
    } else if (is_exponential_gamma) {
      ps0 = ps1 = rgamma_exp(n[k], MS1, PHI[4])
      pu0 = rgamma_exp(n[k], MU0, PHI[2])
      pu1 = rgamma_exp(n[k], MU1, PHI[1])
    } else if (is_odds_gamma) {
      ps0 = ps1 = rodds_gamma(n[k], MS1, PHI[4])
      pu0 = rodds_gamma(n[k], MU0, PHI[2])
      pu1 = rodds_gamma(n[k], MU1, PHI[1])
    } else if (is_unit_lognormal) {
      ps0 = ps1 = runit_lognormal(n[k], MS1, PHI[4])
      pu0 = runit_lognormal(n[k], MU0, PHI[2])
      pu1 = runit_lognormal(n[k], MU1, PHI[1])
    } else if (is_simplex) {
      ps0 = ps1 = rsimplex(n[k], MS1, PHI[4])
      pu0 = rsimplex(n[k], MU0, PHI[2])
      pu1 = rsimplex(n[k], MU1, PHI[1])
    }
    
    while(any(ps1-pu1 <= ps0 - pu0 | ps1<=pu1 | pu0<=pu1)){
      bar = ps1-pu1 <= ps0 - pu0 | ps1<=pu1 | pu0<=pu1
      foo = sum(bar)
      if (is_beta) {
        pu0[bar] = rbeta(foo, MU0 * PHI[2], (1 - MU0) * PHI[2])
        pu1[bar] = rbeta(foo, MU1 * PHI[1], (1 - MU1) * PHI[1])
      } else if (is_unif) {
        pu0[bar] = runif(foo, min = MU0-radius[2], max = MU0+radius[2])
        pu1[bar] = runif(foo, min = MU1-radius[1], max = MU1+radius[1]) 
      } else if (is_logit_normal) {
        pu0[bar] = rlogitnorm(foo, MU0, PHI[2])
        pu1[bar] = rlogitnorm(foo, MU1, PHI[1])
      } else if (is_odds_exponential) {
        pu0[bar] = rodds_exp(foo, MU0)
        pu1[bar] = rodds_exp(foo, MU1)
      } else if (is_exponential_gamma) {
        pu0[bar] = rgamma_exp(foo, MU0, PHI[2])
        pu1[bar] = rgamma_exp(foo, MU1, PHI[1])
      } else if (is_odds_gamma) {
        pu0[bar] = rodds_gamma(foo, MU0, PHI[2])
        pu1[bar] = rodds_gamma(foo, MU1, PHI[1])
      } else if (is_unit_lognormal) {
        pu0[bar] = runit_lognormal(foo, MU0, PHI[2])
        pu1[bar] = runit_lognormal(foo, MU1, PHI[1])
      } else if (is_simplex) {
        pu0[bar] = rsimplex(foo, MU0, PHI[2])
        pu1[bar] = rsimplex(foo, MU1, PHI[1])
      }
    }
    PU1=c(PU1,pu1); PU0=c(PU0,pu0); PS1=c(PS1,ps1); PS0=c(PS0,ps0)
  }
  
  # ==========================================
  # Component 4: u1 = u0
  # ==========================================
  k = 4
  if(n[k] > 0) {
    if (is_beta) {
      pu0 = pu1 = rbeta(n[k], MU0*PHI[2], (1-MU0)*PHI[2])
      ps1 = rbeta(n[k], MS1*PHI[4], (1-MS1)*PHI[4])
      ps0 = rbeta(n[k], MS0*PHI[3], (1-MS0)*PHI[3])
    } else if (is_unif) {
      pu0 = pu1 = runif(n[k], min = MU0-radius[2], max = MU0+radius[2])
      ps1 = runif(n[k], min = MS1-radius[4], max = MS1+radius[4])
      ps0 = runif(n[k], min = MS0-radius[3], max = MS0+radius[3])
    } else if (is_logit_normal) {
      pu0 = pu1 = rlogitnorm(n[k], MU0, PHI[2])
      ps1 = rlogitnorm(n[k], MS1, PHI[4])
      ps0 = rlogitnorm(n[k], MS0, PHI[3])
    } else if (is_odds_exponential) {
      pu0 = pu1 = rodds_exp(n[k], MU0)
      ps1 = rodds_exp(n[k], MS1)
      ps0 = rodds_exp(n[k], MS0)
    } else if (is_exponential_gamma) {
      pu0 = pu1 = rgamma_exp(n[k], MU0, PHI[2])
      ps1 = rgamma_exp(n[k], MS1, PHI[4])
      ps0 = rgamma_exp(n[k], MS0, PHI[3])
    } else if (is_odds_gamma) {
      pu0 = pu1 = rodds_gamma(n[k], MU0, PHI[2])
      ps1 = rodds_gamma(n[k], MS1, PHI[4])
      ps0 = rodds_gamma(n[k], MS0, PHI[3])
    } else if (is_unit_lognormal) {
      pu0 = pu1 = runit_lognormal(n[k], MU0, PHI[2])
      ps1 = runit_lognormal(n[k], MS1, PHI[4])
      ps0 = runit_lognormal(n[k], MS0, PHI[3])
    } else if (is_simplex) {
      pu0 = pu1 = rsimplex(n[k], MU0, PHI[2])
      ps1 = rsimplex(n[k], MS1, PHI[4])
      ps0 = rsimplex(n[k], MS0, PHI[3])
    }
    
    while(any(ps1-pu1 <= ps0 - pu0 | ps1<=pu1 | ps1<=ps0)){
      bar = ps1-pu1 <= ps0 - pu0 | ps1<=pu1 | ps1<=ps0
      foo = sum(bar)
      if (is_beta) {
        ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
        ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
      } else if (is_unif) {
        ps0[bar] = runif(foo, min = MS0-radius[3], max = MS0+radius[3])
        ps1[bar] = runif(foo, min = MS1-radius[4], max = MS1+radius[4])
      } else if (is_logit_normal) {
        ps0[bar] = rlogitnorm(foo, MS0, PHI[3])
        ps1[bar] = rlogitnorm(foo, MS1, PHI[4])
      } else if (is_odds_exponential) {
        ps0[bar] = rodds_exp(foo, MS0)
        ps1[bar] = rodds_exp(foo, MS1)
      } else if (is_exponential_gamma) {
        ps0[bar] = rgamma_exp(foo, MS0, PHI[3])
        ps1[bar] = rgamma_exp(foo, MS1, PHI[4])
      } else if (is_odds_gamma) {
        ps0[bar] = rodds_gamma(foo, MS0, PHI[3])
        ps1[bar] = rodds_gamma(foo, MS1, PHI[4])
      } else if (is_unit_lognormal) {
        ps0[bar] = runit_lognormal(foo, MS0, PHI[3])
        ps1[bar] = runit_lognormal(foo, MS1, PHI[4])
      } else if (is_simplex) {
        ps0[bar] = rsimplex(foo, MS0, PHI[3])
        ps1[bar] = rsimplex(foo, MS1, PHI[4])
      }
    }
    PU1=c(PU1,pu1); PU0=c(PU0,pu0); PS1=c(PS1,ps1); PS0=c(PS0,ps0)
  }
  
  # ==========================================
  # Component 5: s0 = u0, s1 = u1
  # ==========================================
  k = 5
  if(n[k] > 0) {
    if (is_beta) {
      ps1 = pu1 = rbeta(n[k], MU1*PHI[1], (1-MU1)*PHI[1])
      ps0 = pu0 = rbeta(n[k], MU0*PHI[2], (1-MU0)*PHI[2])
    } else if (is_unif) {
      ps1 = pu1 = runif(n[k], min = MU1-radius[1], max = MU1+radius[1])
      ps0 = pu0 = runif(n[k], min = MU0-radius[2], max = MU0+radius[2])
    } else if (is_logit_normal) {
      ps1 = pu1 = rlogitnorm(n[k], MU1, PHI[1])
      ps0 = pu0 = rlogitnorm(n[k], MU0, PHI[2])
    } else if (is_odds_exponential) {
      ps1 = pu1 = rodds_exp(n[k], MU1)
      ps0 = pu0 = rodds_exp(n[k], MU0)
    } else if (is_exponential_gamma) {
      ps1 = pu1 = rgamma_exp(n[k], MU1, PHI[1])
      ps0 = pu0 = rgamma_exp(n[k], MU0, PHI[2])
    } else if (is_odds_gamma) {
      ps1 = pu1 = rodds_gamma(n[k], MU1, PHI[1])
      ps0 = pu0 = rodds_gamma(n[k], MU0, PHI[2])
    } else if (is_unit_lognormal) {
      ps1 = pu1 = runit_lognormal(n[k], MU1, PHI[1])
      ps0 = pu0 = runit_lognormal(n[k], MU0, PHI[2])
    } else if (is_simplex) {
      ps1 = pu1 = rsimplex(n[k], MU1, PHI[1])
      ps0 = pu0 = rsimplex(n[k], MU0, PHI[2])
    }
    PU1=c(PU1,pu1); PU0=c(PU0,pu0); PS1=c(PS1,ps1); PS0=c(PS0,ps0)
  }
  
  # ==========================================
  # Component 6: s1 = u1
  # ==========================================
  k = 6
  if(n[k] > 0) {
    if (is_beta) {
      ps1 = pu1 = rbeta(n[k], MU1*PHI[1], (1-MU1)*PHI[1])
      ps0 = rbeta(n[k], MS0*PHI[3], (1-MS0)*PHI[3])
      pu0 = rbeta(n[k], MU0*PHI[2], (1-MU0)*PHI[2])
    } else if (is_unif) {
      ps1 = pu1 = runif(n[k], min = MU1-radius[1], max = MU1+radius[1])
      ps0 = runif(n[k], min = MS0-radius[3], max = MS0+radius[3])
      pu0 = runif(n[k], min = MU0-radius[2], max = MU0+radius[2])
    } else if (is_logit_normal) {
      ps1 = pu1 = rlogitnorm(n[k], MU1, PHI[1])
      ps0 = rlogitnorm(n[k], MS0, PHI[3])
      pu0 = rlogitnorm(n[k], MU0, PHI[2])
    } else if (is_odds_exponential) {
      ps1 = pu1 = rodds_exp(n[k], MU1)
      ps0 = rodds_exp(n[k], MS0)
      pu0 = rodds_exp(n[k], MU0)
    } else if (is_exponential_gamma) {
      ps1 = pu1 = rgamma_exp(n[k], MU1, PHI[1])
      ps0 = rgamma_exp(n[k], MS0, PHI[3])
      pu0 = rgamma_exp(n[k], MU0, PHI[2])
    } else if (is_odds_gamma) {
      ps1 = pu1 = rodds_gamma(n[k], MU1, PHI[1])
      ps0 = rodds_gamma(n[k], MS0, PHI[3])
      pu0 = rodds_gamma(n[k], MU0, PHI[2])
    } else if (is_unit_lognormal) {
      ps1 = pu1 = runit_lognormal(n[k], MU1, PHI[1])
      ps0 = runit_lognormal(n[k], MS0, PHI[3])
      pu0 = runit_lognormal(n[k], MU0, PHI[2])
    } else if (is_simplex) {
      ps1 = pu1 = rsimplex(n[k], MU1, PHI[1])
      ps0 = rsimplex(n[k], MS0, PHI[3])
      pu0 = rsimplex(n[k], MU0, PHI[2])
    }
    
    while(any(ps0 < pu0)){
      bar = ps0 < pu0
      foo = sum(bar)
      if (is_beta) {
        ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
        pu0[bar] = rbeta(foo, MU0 * PHI[2], (1 - MU0) * PHI[2])
      } else if (is_unif) {
        ps0[bar] = runif(foo, min = MS0-radius[3], max = MS0+radius[3])
        pu0[bar] = runif(foo, min = MU0-radius[2], max = MU0+radius[2])
      } else if (is_logit_normal) {
        ps0[bar] = rlogitnorm(foo, MS0, PHI[3])
        pu0[bar] = rlogitnorm(foo, MU0, PHI[2])
      } else if (is_odds_exponential) {
        ps0[bar] = rodds_exp(foo, MS0)
        pu0[bar] = rodds_exp(foo, MU0)
      } else if (is_exponential_gamma) {
        ps0[bar] = rgamma_exp(foo, MS0, PHI[3])
        pu0[bar] = rgamma_exp(foo, MU0, PHI[2])
      } else if (is_odds_gamma) {
        ps0[bar] = rodds_gamma(foo, MS0, PHI[3])
        pu0[bar] = rodds_gamma(foo, MU0, PHI[2])
      } else if (is_unit_lognormal) {
        ps0[bar] = runit_lognormal(foo, MS0, PHI[3])
        pu0[bar] = runit_lognormal(foo, MU0, PHI[2])
      } else if (is_simplex) {
        ps0[bar] = rsimplex(foo, MS0, PHI[3])
        pu0[bar] = rsimplex(foo, MU0, PHI[2])
      }
    }
    PU1=c(PU1,pu1); PU0=c(PU0,pu0); PS1=c(PS1,ps1); PS0=c(PS0,ps0)
  }
  
  # ==========================================
  # Component 7: s1 = u1 = s0 = u0
  # ==========================================
  k = 7
  if(n[k] > 0) {
    if (is_beta) {
      ps1=ps0=pu1=pu0 = rbeta(n[k], MU0*PHI[2], (1-MU0)*PHI[2])
    } else if (is_unif) {
      ps1=ps0=pu1=pu0 = runif(n[k], min = MU0-radius[2], max = MU0+radius[2])
    } else if (is_logit_normal) {
      ps1=ps0=pu1=pu0 = rlogitnorm(n[k], MU0, PHI[2])
    } else if (is_odds_exponential) {
      ps1=ps0=pu1=pu0 = rodds_exp(n[k], MU0)
    } else if (is_exponential_gamma) {
      ps1=ps0=pu1=pu0 = rgamma_exp(n[k], MU0, PHI[2])
    } else if (is_odds_gamma) {
      ps1=ps0=pu1=pu0 = rodds_gamma(n[k], MU0, PHI[2])
    } else if (is_unit_lognormal) {
      ps1=ps0=pu1=pu0 = runit_lognormal(n[k], MU0, PHI[2])
    } else if (is_simplex) {
      ps1=ps0=pu1=pu0 = rsimplex(n[k], MU0, PHI[2])
    }
    PU1=c(PU1,pu1); PU0=c(PU0,pu0); PS1=c(PS1,ps1); PS0=c(PS0,ps0)
  }
  
  # ==========================================
  # Component 8: s1 = s0, u1 = u0
  # ==========================================
  k = 8
  if(n[k] > 0){
    if (is_beta) {
      ps0 = ps1 = rbeta(n[k], MS0*PHI[3], (1-MS0)*PHI[3])
      pu0 = pu1 = rbeta(n[k], MU0*PHI[2], (1-MU0)*PHI[2])
    } else if (is_unif) {
      ps0 = ps1 = runif(n[k], min = MS0-radius[3], max = MS0+radius[3])
      pu0 = pu1 = runif(n[k], min = MU0-radius[2], max = MU0+radius[2])
    } else if (is_logit_normal) {
      ps0 = ps1 = rlogitnorm(n[k], MS0, PHI[3])
      pu0 = pu1 = rlogitnorm(n[k], MU0, PHI[2])
    } else if (is_odds_exponential) {
      ps0 = ps1 = rodds_exp(n[k], MS0)
      pu0 = pu1 = rodds_exp(n[k], MU0)
    } else if (is_exponential_gamma) {
      ps0 = ps1 = rgamma_exp(n[k], MS0, PHI[3])
      pu0 = pu1 = rgamma_exp(n[k], MU0, PHI[2])
    } else if (is_odds_gamma) {
      ps0 = ps1 = rodds_gamma(n[k], MS0, PHI[3])
      pu0 = pu1 = rodds_gamma(n[k], MU0, PHI[2])
    } else if (is_unit_lognormal) {
      ps0 = ps1 = runit_lognormal(n[k], MS0, PHI[3])
      pu0 = pu1 = runit_lognormal(n[k], MU0, PHI[2])
    } else if (is_simplex) {
      ps0 = ps1 = rsimplex(n[k], MS0, PHI[3])
      pu0 = pu1 = rsimplex(n[k], MU0, PHI[2])
    }
    PU1=c(PU1,pu1); PU0=c(PU0,pu0); PS1=c(PS1,ps1); PS0=c(PS0,ps0)
  }
  
  # ==========================================
  # Count Simulation & Mapping
  # ==========================================
  colnames(Ntot) = c("nu1", "ns1", "nu0", "ns0")
  nu1 = rbinom(P, Ntot[, "nu1"], PU1)
  ns1 = rbinom(P, Ntot[, "ns1"], PS1)
  nu0 = rbinom(P, Ntot[, "nu0"], PU0)
  ns0 = rbinom(P, Ntot[, "ns0"], PS0)
  
  truth = rep(c("R1","R2","R3","R4","NR1","NR2","NR3","NSR"), n)
  
  return(list(Ntot=Ntot, ns0=ns0, ns1=ns1, nu0=nu0, nu1=nu1, truth=truth))
}

DiD_GLM <- function(Ntot, ns1, nu1, ns0, nu0) {
  P <- nrow(Ntot)
  ind_ids <- 1:P
  df <- data.frame(
    Individual = rep(ind_ids, times = 4),
    Time_Point = rep(c("Active", "Active", "Baseline", "Baseline"), each = P),
    Stimulation = rep(c("Stimulated", "Unstimulated", "Stimulated", "Unstimulated"), each = P),
    Positives = c(ns1, nu1, ns0, nu0),
    Total = c(Ntot[, 2], Ntot[, 1], Ntot[, 4], Ntot[, 3])
  )
  df$Negatives <- df$Total - df$Positives
  df$Time_Point  <- factor(df$Time_Point, levels = c("Baseline", "Active"))
  df$Stimulation <- factor(df$Stimulation, levels = c("Unstimulated", "Stimulated"))
  
  responder_probs <- rep(NA,P)
  for (i in 1:P) {
    test_mod <- glm(cbind(Positives,Negatives)~Time_Point*Stimulation,family = "binomial",data = df[df$Individual==i,]) |> summary()
    
    if (test_mod$coefficients["Time_PointActive:StimulationStimulated","Estimate"]>=0){
      responder_probs[i] <- (test_mod$coefficients["Time_PointActive:StimulationStimulated",c("Pr(>|z|)")]/2)
    } else{
      responder_probs[i] <- 0.95
    }
  }
  
  responder_probs <- 1-responder_probs
  return(responder_probs)
}
