my_libs <- "/scratch/abrmoe030/R_libs"
.libPaths(c(my_libs, .libPaths()))

library(MIMOSA2)
library(ggplot2)
library(ashr)

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
  
n = floor(P * pis)
  
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
    Individual  = rep(ind_ids, times = 4),
    Time_Point  = rep(c("Active", "Active", "Baseline", "Baseline"), each = P),
    Stimulation = rep(c("Stimulated", "Unstimulated", "Stimulated", "Unstimulated"), each = P),
    Positives   = c(ns1, nu1, ns0, nu0),
    Total       = c(Ntot[, "ns1"], Ntot[, "nu1"], Ntot[, "ns0"], Ntot[, "nu0"])
  )
  df$Negatives   <- df$Total - df$Positives
  df$Time_Point  <- factor(df$Time_Point, levels = c("Baseline", "Active"))
  df$Stimulation <- factor(df$Stimulation, levels = c("Unstimulated", "Stimulated"))
  
  responder_probs <- numeric(P)
  
  for (i in 1:P) {
    df_ind <- df[df$Individual == i, ]
    
    # Compute starting values from observed cell proportions
    p_AS <- df_ind$Positives[1] / df_ind$Total[1]
    p_AU <- df_ind$Positives[2] / df_ind$Total[2]
    p_BS <- df_ind$Positives[3] / df_ind$Total[3]
    p_BU <- df_ind$Positives[4] / df_ind$Total[4]
    
    # Nudge exact 0s and 1s inside (0, 1)
    eps <- 1e-5
    p_AS <- pmax(pmin(p_AS, 1 - eps), eps)
    p_AU <- pmax(pmin(p_AU, 1 - eps), eps)
    p_BS <- pmax(pmin(p_BS, 1 - eps), eps)
    p_BU <- pmax(pmin(p_BU, 1 - eps), eps)
    
    start_vals <- c(
      "(Intercept)"                             = p_BU,
      "Time_PointActive"                        = p_AU - p_BU,
      "StimulationStimulated"                   = p_BS - p_BU,
      "Time_PointActive:StimulationStimulated" = (p_AS - p_AU) - (p_BS - p_BU)
    )
    
    # Suppress temporary step-halving boundary warnings inside the loop
    fit <- suppressWarnings(
      tryCatch({
        glm(
          cbind(Positives, Negatives) ~ Time_Point * Stimulation,
          family  = binomial(link = "identity"),
          data    = df_ind,
          start   = start_vals,
          control = glm.control(maxit = 200, epsilon = 1e-8)
        )
      }, error = function(e) NULL)
    )
    
    # Check convergence
    if (is.null(fit) || !fit$converged) {
      responder_probs[i] <- 0.95
      next
    }
    
    coef_matrix <- summary(fit)$coefficients
    param_name  <- "Time_PointActive:StimulationStimulated"
    
    if (param_name %in% rownames(coef_matrix)) {
      est  <- coef_matrix[param_name, "Estimate"]
      pval <- coef_matrix[param_name, "Pr(>|z|)"]
      
      if (!is.na(est) && est >= 0) {
        responder_probs[i] <- pval / 2
      } else {
        responder_probs[i] <- 0.95
      }
    } else {
      responder_probs[i] <- 0.95
    }
  }
  
  return(1 - responder_probs)
}

# ==============================================================================
# Helper: Numerically Stable Log-Sum-Exp for Vectorized Inputs
# ==============================================================================
log_sum_exp_2 <- function(a, b) {
  max_val <- pmax(a, b)
  max_val + log(exp(a - max_val) + exp(b - max_val))
}

# ==============================================================================
# Helper: Empirical Bayes Variance Shrinkage
# Shrinks noisy individual sample variances toward a global inverse-Gamma prior.
# Prevents low-count/zero-event gates from creating zero or infinite variances.
# ==============================================================================
shrink_variances_eb <- function(var_D, df_raw = 4, floor_val = 1e-8) {
  # Floor raw variances to prevent zeroes
  var_safe <- pmax(var_D, floor_val)
  
  # Log-transform variances for robust moment estimation
  log_v <- log(var_safe)
  mean_log_v <- mean(log_v, na.rm = TRUE)
  var_log_v  <- var(log_v, na.rm = TRUE)
  
  # Method of moments estimation for inverse-Gamma prior parameters (d0, s02)
  # Trigonometric psi-function variance approximation: trigamma(d0 / 2) ~ var_log_v
  d0 <- max(1, 2 / max(0.01, var_log_v - trigamma(df_raw / 2)))
  s02 <- exp(mean_log_v)
  
  # Moderated variance: Weighted average of raw variance and prior target
  var_shrunk <- (d0 * s02 + df_raw * var_safe) / (d0 + df_raw)
  return(pmax(var_shrunk, floor_val))
}

# ==============================================================================
# Main Model: Upgraded EM Mixture Model with Variance Shrinkage & Log-Likelihood
# ==============================================================================
fit_shrunk_mixture_EM <- function(D, var_D, max_iter = 200, tol = 1e-6, prior = 0.1, use_eb = TRUE) {
  N <- length(D)
  
  # 1. Apply Empirical Bayes shrinkage if requested
  var_clean <- if (use_eb) shrink_variances_eb(var_D) else pmax(var_D, 1e-8)
  sd_D <- sqrt(var_clean)
  
  # 2. Smart Initialization
  p <- pmin(pmax(prior, 0.01), 0.99)
  pos_diffs <- D[D > 0]
  mu <- if (length(pos_diffs) > 0) max(1e-4, quantile(pos_diffs, 0.75, na.rm = TRUE)) else 1e-3
  
  log_lik_history <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    p_old  <- p
    mu_old <- mu
    
    # --- E-Step (In Log-Space) ---
    log_f0 <- dnorm(D, mean = 0,  sd = sd_D, log = TRUE)
    log_f1 <- dnorm(D, mean = mu, sd = sd_D, log = TRUE)
    
    log_num0 <- log(1 - p) + log_f0
    log_num1 <- log(p)     + log_f1
    
    # Marginal log-likelihood per observation
    log_marginal <- log_sum_exp_2(log_num0, log_num1)
    
    # Track total log-likelihood across cohort
    log_lik_history[iter] <- sum(log_marginal)
    
    # Posterior responsibilities
    log_gamma <- log_num1 - log_marginal
    gamma     <- exp(log_gamma)
    gamma[is.na(gamma)] <- 0
    
    # --- M-Step ---
    # Update responder proportion
    p <- mean(gamma)
    p <- pmin(pmax(p, 1e-12), 1 - 1e-12) # Clamp away from 0/1 boundary
    
    # Update responder mean (Weighted Least Squares with shrunk precision weights)
    weights <- gamma / var_clean
    sum_w   <- sum(weights)
    mu      <- if (sum_w > 0) max(1e-6, sum(weights * D) / sum_w) else 1e-6
    
    # --- Convergence Check ---
    param_diff <- max(abs(p - p_old), abs(mu - mu_old))
    if (param_diff < tol) {
      log_lik_history <- log_lik_history[1:iter]
      break
    }
  }
  
  return(list(
    p               = p,
    mu              = mu,
    prob_responder  = gamma,
    var_shrunk      = var_clean,
    iterations      = iter,
    log_lik_history = log_lik_history
  ))
}

# ==============================================================================
# DiD Wrapper Function for Simulation Pipelines
# ==============================================================================
DiD_mixture_shrunk <- function(Ntot, ns1, nu1, ns0, nu0, use_eb = TRUE) {
  
  Ntot_mat <- as.matrix(Ntot)
  N_AS <- as.numeric(Ntot_mat[, "ns1"]); N_AU <- as.numeric(Ntot_mat[, "nu1"])
  N_BS <- as.numeric(Ntot_mat[, "ns0"]); N_BU <- as.numeric(Ntot_mat[, "nu0"])
  
  ns1 <- as.numeric(ns1); nu1 <- as.numeric(nu1)
  ns0 <- as.numeric(ns0); nu0 <- as.numeric(nu0)
  
  # Raw proportions
  p_AS <- ns1 / N_AS
  p_AU <- nu1 / N_AU
  p_BS <- ns0 / N_BS
  p_BU <- nu0 / N_BU
  
  DiD <- (p_AS - p_AU) - (p_BS - p_BU)
  
  # Jeffreys variance estimation for counts
  v_AS <- (ns1 + 0.5) / (N_AS + 1)
  v_AU <- (nu1 + 0.5) / (N_AU + 1)
  v_BS <- (ns0 + 0.5) / (N_BS + 1)
  v_BU <- (nu0 + 0.5) / (N_BU + 1)
  
  var_DiD <- v_AS * (1 - v_AS) / N_AS +
    v_AU * (1 - v_AU) / N_AU +
    v_BS * (1 - v_BS) / N_BS +
    v_BU * (1 - v_BU) / N_BU
  
  fit_shrunk_mixture_EM(DiD, var_DiD, use_eb = use_eb)
}

