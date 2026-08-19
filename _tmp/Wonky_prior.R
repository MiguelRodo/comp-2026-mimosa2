rbernoulli_beta_mixture = function(num_draws, target_mean, component_phi, p_bernoulli = 0.25) {
  if (num_draws <= 0) return(numeric(0))
  
  # 1. Target Marginal Variance of Beta(target_mean, component_phi)
  var_marginal = (target_mean * (1 - target_mean)) / (component_phi + 1)
  
  # 2. Target Conditional Variance (halved)
  var_cond = 0.125 * var_marginal
  
  # 3. Between-component variance needed: Var_between = Var_marginal - Var_cond
  var_between = var_marginal - var_cond
  
  # 4. Difference between conditional means: delta = mu1 - mu0
  delta = sqrt(var_between / (p_bernoulli * (1 - p_bernoulli)))
  
  # 5. Conditional Means
  mu1 = target_mean + (1 - p_bernoulli) * delta
  mu0 = target_mean - p_bernoulli * delta
  
  # Boundary check: ensure conditional means lie strictly within (0, 1)
  if (mu0 <= 0 || mu1 >= 1) {
    stop(sprintf(
      "Variance too large relative to mean (mu = %e). Conditional means [%e, %e] exceed (0, 1).",
      target_mean, mu0, mu1
    ))
  }
  
  # 6. Solve precision parameters phi_1 and phi_0 for exact conditional variance
  # Using formula: Var_cond = mu_k * (1 - mu_k) / (phi_k + 1)
  phi1 = (mu1 * (1 - mu1) / var_cond) - 1
  phi0 = (mu0 * (1 - mu0) / var_cond) - 1
  
  # 7. Sample Z ~ Bernoulli(p)
  z = rbinom(num_draws, size = 1, prob = p_bernoulli)
  
  # 8. Sample X | Z ~ Beta(alpha_k, beta_k)
  draws = numeric(num_draws)
  n1 = sum(z == 1)
  n0 = num_draws - n1
  
  if (n1 > 0) {
    draws[z == 1] = rbeta(n1, shape1 = mu1 * phi1, shape2 = (1 - mu1) * phi1)
  }
  if (n0 > 0) {
    draws[z == 0] = rbeta(n0, shape1 = mu0 * phi0, shape2 = (1 - mu0) * phi0)
  }
  
  return(draws)
}

obj = rbernoulli_beta_mixture(num_draws=100000,target_mean=6e-4,component_phi=10000)

hist(obj,breaks = 100)
