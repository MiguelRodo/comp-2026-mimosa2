match_gamma_dispersion <- function(mu, phi_beta = 5000) {
  target_var <- (mu * (1 - mu)) / (1 + phi_beta)
  
  variance_diff <- function(k) {
    if (k <= 0) return(Inf)
    mu_G <- k * (mu^(-1/k) - 1)
    # E[X^2] for rgamma_exp formulation
    e_x2 <- (1 + 2 * mu_G / k)^(-k)
    current_var <- e_x2 - mu^2
    return(current_var - target_var)
  }
  
  res <- uniroot(variance_diff, lower = 0.1, upper = 100000)
  return(res$root)
}

# Example for baseline background mean:
mu_0 <- 1e-4
equivalent_k <- match_gamma_dispersion(mu_0, phi_beta = 5000)
print(equivalent_k)


rgamma_exp = function(num_draws, target_mean, gamma_dispersion){
  gamma_mean <- gamma_dispersion*(target_mean^(-1/gamma_dispersion)-1)
  outs <- exp(-rgamma(num_draws, gamma_dispersion, gamma_dispersion/gamma_mean))
  return(outs)
}

rgamma_exp(100000,mu_0,68) %>% var()

match_logitnorm_phi <- function(mu, phi_beta = 5000) {
  # 1. Target Beta variance given precision phi_beta
  target_var <- (mu * (1 - mu)) / (1 + phi_beta)
  
  # Logit location parameter
  meanlogit <- qlogis(mu)
  
  # Objective function: difference between Logit-Normal variance and target variance
  variance_diff <- function(sdlogit) {
    if (sdlogit <= 0) return(Inf)
    
    # Integrand for E[Y^2] under Logit-Normal
    integrand_e2 <- function(x) {
      y <- plogis(x)
      y^2 * dnorm(x, mean = meanlogit, sd = sdlogit)
    }
    
    # Numerical integration with adaptive limits around meanlogit
    e_y2 <- integrate(
      integrand_e2, 
      lower = meanlogit - 10 * sdlogit, 
      upper = meanlogit + 10 * sdlogit
    )$value
    
    current_var <- e_y2 - mu^2
    return(current_var - target_var)
  }
  
  # Search for the matching standard deviation on logit scale
  # Search interval corresponds to SD range [1e-4, 5]
  res <- uniroot(variance_diff, lower = 1e-4, upper = 5.0)
  
  # Convert SD to your function's precision parameter: component_phi = 1 / (sdlogit^2)
  equivalent_phi <- 1 / (res$root^2)
  
  return(list(
    sdlogit = res$root,
    equivalent_phi = equivalent_phi
  ))
}

match_logitnorm_phi(mu_0,5000)

rlogitnorm = function(num_draws, target_mean, component_phi) {
  if (num_draws <= 0) return(numeric(0))
  sdlogit   = 1 / sqrt(component_phi)
  meanlogit = qlogis(target_mean)
  draws     = plogis(rnorm(num_draws, mean = meanlogit, sd = sdlogit))
  return(draws)
}

rlogitnorm(100000,mu_0,1.45) %>% var()

match_simplex_phi <- function(mu, phi_beta = 5000) {
  # Target Beta variance
  target_var <- (mu * (1 - mu)) / (1 + phi_beta)
  
  # Base odds parameter
  mu_ig <- mu / (1 - mu)
  
  # Objective function comparing Simplex/IG variance to target Beta variance
  variance_diff <- function(phi) {
    if (phi <= 0) return(Inf)
    
    # Lambda scale parameter defined in your rsimplex formulation
    lambda_ig <- phi * (mu^2) * ((1 - mu)^2)
    
    # Inverse Gaussian Density function
    dig <- function(x) {
      suppressWarnings({
        dens <- sqrt(lambda_ig / (2 * pi * x^3)) * 
          exp(- (lambda_ig * (x - mu_ig)^2) / (2 * (mu_ig^2) * x))
      })
      # Replace any numerical NaN/Inf at boundary limits with 0
      dens[!is.finite(dens)] <- 0
      return(dens)
    }
    
    # Integrand for E[Y^2] where Y = X / (1 + X)
    integrand_e2 <- function(x) {
      y <- x / (1 + x)
      y^2 * dig(x)
    }
    
    # Integrate around the mode/mean on the odds scale
    # Integrates over a wide spread relative to mu_ig
    e_y2 <- integrate(
      integrand_e2, 
      lower = 0, 
      upper = mu_ig * 1000, 
      rel.tol = 1e-8,
      stop.on.error = FALSE
    )$value
    
    current_var <- e_y2 - mu^2
    return(current_var - target_var)
  }
  
  # Search for matching component_phi over a wide precision range
  res <- uniroot(variance_diff, lower = 0.1, upper = 1e7)
  
  return(res$root)
}

match_simplex_phi(mu_0,5000)

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

rsimplex(100000,mu_0,5000) %>% var()
