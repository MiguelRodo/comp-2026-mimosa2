rbernoulli_beta_bimodal_wide = function(num_draws, target_mean, component_phi, 
                                        p_bernoulli = 0.20,             # 20% high responders, 80% dominant low responders
                                        variance_split_ratio = 0.65,    # Allocation to mean distance
                                        peak_ratio = 0.15) {            # Allocation of internal variance to dominant peak
  if (num_draws <= 0) return(numeric(0))
  
  # 1. Total Marginal Variance anchored to target_mean and component_phi
  var_marginal = (target_mean * (1 - target_mean)) / (component_phi + 1)
  
  # 2. Split into between-component and within-component variance
  var_between = variance_split_ratio * var_marginal
  var_within  = (1 - variance_split_ratio) * var_marginal
  
  # 3. Component Means (Distance driven by between-component variance)
  delta = sqrt(var_between / (p_bernoulli * (1 - p_bernoulli)))
  mu1   = target_mean + (1 - p_bernoulli) * delta  # High responder mean
  mu0   = target_mean - p_bernoulli * delta        # Dominant low responder mean
  
  # Boundary protection
  if (mu0 <= 1e-7 || mu1 >= (1 - 1e-7)) {
    max_delta   = min(target_mean / p_bernoulli, (1 - target_mean) / (1 - p_bernoulli)) * 0.95
    delta       = max_delta
    mu1         = target_mean + (1 - p_bernoulli) * delta
    mu0         = target_mean - p_bernoulli * delta
    var_between = (delta^2) * (p_bernoulli * (1 - p_bernoulli))
    var_within  = var_marginal - var_between
  }
  
  # 4. Asymmetric Variance Allocation for Asymmetric Peak Heights
  # dominant group (0) gets less variance (taller peak), minor group (1) gets more (flatter peak)
  var_0 = (peak_ratio * var_within) / (1 - p_bernoulli)
  var_1 = ((1 - peak_ratio) * var_within) / p_bernoulli
  
  # 5. Convert component variances into Beta shape parameters (phi)
  phi0 = (mu0 * (1 - mu0) / var_0) - 1
  phi1 = (mu1 * (1 - mu1) / var_1) - 1
  
  # Ensure numerical stability
  phi0 = max(phi0, 2)
  phi1 = max(phi1, 2)
  
  # 6. Sampling
  z     = rbinom(num_draws, size = 1, prob = p_bernoulli)
  draws = numeric(num_draws)
  n1    = sum(z == 1)
  n0    = num_draws - n1
  
  if (n1 > 0) draws[z == 1] = rbeta(n1, shape1 = mu1 * phi1, shape2 = (1 - mu1) * phi1)
  if (n0 > 0) draws[z == 0] = rbeta(n0, shape1 = mu0 * phi0, shape2 = (1 - mu0) * phi0)
  
  return(draws)
}














library(ggplot2)

n_sims      <- 200000
target_mean <- 0.00075
phi_param   <- 10000

# Compute expected target statistics
target_var <- (target_mean * (1 - target_mean)) / (phi_param + 1)

# Generate sample
wide_asymmetric_draws <- rbernoulli_beta_bimodal_wide(
  num_draws            = n_sims, 
  target_mean          = target_mean, 
  component_phi        = phi_param,
  p_bernoulli          = 0.15,  # 85% dominant peak vs 15% secondary peak
  variance_split_ratio = 0.70,  # Keeps peaks well-separated
  peak_ratio           = 0.10   # Makes dominant peak very sharp & high
)

# Verify empirical adherence to mean and variance constraints
cat("Target Mean:   ", target_mean,  "| Empirical Mean:   ", mean(wide_asymmetric_draws), "\n")
cat("Target Variance:", target_var,   "| Empirical Variance:", var(wide_asymmetric_draws), "\n")

# Plot Density
df <- data.frame(Value = wide_asymmetric_draws)

ggplot(df, aes(x = Value)) +
  geom_density(fill = "#2B5C8F", alpha = 0.6, color = "#1A365D", linewidth = 0.8) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Asymmetric Bimodal Prior",
    subtitle = "Distinct wide-separated peaks with conserved marginal mean & variance",
    x = "Proportion / Effect Size",
    y = "Density"
  )
