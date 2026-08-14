# Bimodal prior visualisation 
# 14 August 2026 

library(ggplot2)

# ==============================================================================
# 1. Bimodal Prior Sampler (Calibrated to Target Variance = 2e-8)
# ==============================================================================
rbimodal_beta_matched <- function(n, target_mean = 6e-4, target_var = 2e-8, separation_ratio = 0.8, mix_weight = 0.5) {
  if (n <= 0) return(numeric(0))
  
  sd_target <- sqrt(target_var)
  delta     <- separation_ratio * sd_target
  mu1       <- target_mean - delta
  mu2       <- target_mean + delta
  
  # Allocate variance across components (Law of Total Variance)
  var_between <- mix_weight * (1 - mix_weight) * (2 * delta)^2
  var_within  <- target_var - var_between
  
  if (var_within <= 0) stop("Separation ratio too large! Results in negative within-variance.")
  
  # Solve precisions for sub-components
  phi1 <- (mu1 * (1 - mu1)) / var_within - 1
  phi2 <- (mu2 * (1 - mu2)) / var_within - 1
  
  a1 <- mu1 * phi1; b1 <- (1 - mu1) * phi1
  a2 <- mu2 * phi2; b2 <- (1 - mu2) * phi2
  
  # Sample indicator and draws
  z <- rbinom(n, size = 1, prob = 1 - mix_weight)
  draws <- numeric(n)
  
  n1 <- sum(z == 0)
  n2 <- sum(z == 1)
  
  if (n1 > 0) draws[z == 0] <- rbeta(n1, a1, b1)
  if (n2 > 0) draws[z == 1] <- rbeta(n2, a2, b2)
  
  return(draws)
}

# ==============================================================================
# 2. Density Function
# ==============================================================================
dbimodal_beta_matched <- function(x, target_mean = 6e-4, target_var = 2e-8, separation_ratio = 0.8, mix_weight = 0.5) {
  delta       <- separation_ratio * sqrt(target_var)
  mu1         <- target_mean - delta
  mu2         <- target_mean + delta
  var_between <- mix_weight * (1 - mix_weight) * (2 * delta)^2
  var_within  <- target_var - var_between
  
  phi1 <- (mu1 * (1 - mu1)) / var_within - 1
  phi2 <- (mu2 * (1 - mu2)) / var_within - 1
  
  a1 <- mu1 * phi1; b1 <- (1 - mu1) * phi1
  a2 <- mu2 * phi2; b2 <- (1 - mu2) * phi2
  
  return(mix_weight * dbeta(x, a1, b1) + (1 - mix_weight) * dbeta(x, a2, b2))
}

# ==============================================================================
# 3. Parameters & Empirical Verification
# ==============================================================================
MS1        <- 6e-4
target_var <- 2e-8

# Beta precision derived directly from target variance
phi_beta   <- (MS1 * (1 - MS1)) / target_var - 1

set.seed(42)
samples_beta    <- rbeta(100000, MS1 * phi_beta, (1 - MS1) * phi_beta)
samples_bimodal <- rbimodal_beta_matched(100000, target_mean = MS1, target_var = target_var, separation_ratio = 0.8)

cat(sprintf("Baseline Beta  -> Mean: %.7f | Var: %.10f\n", mean(samples_beta), var(samples_beta)))
cat(sprintf("Bimodal Prior  -> Mean: %.7f | Var: %.10f\n", mean(samples_bimodal), var(samples_bimodal)))

# ==============================================================================
# 4. Visualization
# ==============================================================================
x_grid <- seq(1e-5, 0.0012, length.out = 1000)

df_beta <- data.frame(
  x = x_grid, 
  y = dbeta(x_grid, MS1 * phi_beta, (1 - MS1) * phi_beta), 
  Type = "Baseline Beta (Unimodal)"
)

df_bimodal <- data.frame(
  x = x_grid, 
  y = dbimodal_beta_matched(x_grid, target_mean = MS1, target_var = target_var, separation_ratio = 0.8), 
  Type = "Bimodal Mixture Prior"
)

df_plot <- rbind(df_beta, df_bimodal)

ggplot(df_plot, aes(x = x, y = y, color = Type, fill = Type)) +
  geom_line(linewidth = 1.1) +
  geom_area(alpha = 0.15, position = "identity") +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = c("Baseline Beta (Unimodal)" = "#2b5c8f", "Bimodal Mixture Prior" = "#d95f02")) +
  scale_fill_manual(values = c("Baseline Beta (Unimodal)" = "#2b5c8f", "Bimodal Mixture Prior" = "#d95f02")) +
  labs(
    title = "Comparison of Unimodal vs Variance-Matched Bimodal Prior",
    subtitle = expression("Both distributions share identical " ~ mu == 0.0006 ~ " and " ~ sigma^2 == 2 %*% 10^-8),
    x = expression(p[i]^"(s,1)" ~ "(Stimulated Cell Proportion)"),
    y = "Prior Density",
    color = "Prior Specification",
    fill = "Prior Specification"
  ) +
  theme(legend.position = "top")
