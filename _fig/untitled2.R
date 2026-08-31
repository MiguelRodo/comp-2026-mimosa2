library(ggplot2)
library(dplyr)
library(purrr)

set.seed(42)

responders50  = c(rep(0.50 / 4, 4), rep(0.50 / 4, 4))

# Combine into a named list for tracking:
component_list = list(
  "Prop_0.50" = responders50
)

# Cell count scenarios: 
# Combine your rng cell counts into a named list for tracking:
rng_list = list(
  "High"  = c(250000, 250000),
  "Medium" = c(100000, 100000),
  "Low"     = c(15000, 15000),
  "Very Low"= c(7000,7000)
)

# Build simulation grid:
stresstest_mat = expand.grid(
  Distribution_Phi   = list(c("Beta",10000),
                            c("EG",113),
                            c("LN",2.05),
                            c("SX",10000),
                            c("BB",10000)),
  Comp_Name      = names(component_list),
  P              = c(50),
  Effect         = c(1e-3, 2.5e-4, 1.25e-4, 6.25e-5),
  Rng_Name       = names(rng_list),
  Replication    = 1:30, 
  stringsAsFactors = FALSE
)

dist_phi_pairs <- unique(stresstest_mat$Distribution_Phi)
MU_TARGET <- 5e-4 
N_DRAWS   <- 50000

# 1. Define custom attributes for each (Distribution, Phi) pair
# Adjust titles, fill colors, and line colors here
custom_specs <- list(
  "Beta" = list(Title = "Beta",              Fill = "grey30", Color = "grey30"),
  "EG"   = list(Title = "Exponential Gamma", Fill = "steelblue3", Color = "steelblue3"),
  "LN"   = list(Title = "Logit Normal",      Fill = "purple", Color = "purple"),
  "SX"   = list(Title = "Simplex",           Fill = "deeppink", Color = "deeppink"),
  "BB"   = list(Title = "Bimodal Beta",      Fill = "orange", Color = "orange")
)

# 2. Sampler function
sample_prior <- function(prior_alias, phi_val, mu, n_draws) {
  if (prior_alias %in% c("Beta", "beta", "b", "B")) {
    rbeta(n_draws, shape1 = mu * phi_val, shape2 = (1 - mu) * phi_val)
    
  } else if (prior_alias %in% c("EG", "Exponential gamma", "eg")) {
    gamma_mean <- phi_val * (mu^(-1 / phi_val) - 1)
    exp(-rgamma(n_draws, shape = phi_val, rate = phi_val / gamma_mean))
    
  } else if (prior_alias %in% c("LN", "logit normal", "ln")) {
    sdlogit   <- 1 / sqrt(phi_val)
    meanlogit <- qlogis(mu)
    plogis(rnorm(n_draws, mean = meanlogit, sd = sdlogit))
    
  } else if (prior_alias %in% c("SX", "Simplex", "sx")) {
    mu_ig     <- mu / (1 - mu)
    lambda_ig <- phi_val * (mu^2) * ((1 - mu)^2)
    v         <- rnorm(n_draws)^2
    x         <- mu_ig + (mu_ig^2 * v) / (2 * lambda_ig) - 
      (mu_ig / (2 * lambda_ig)) * sqrt(4 * mu_ig * lambda_ig * v + mu_ig^2 * v^2)
    z         <- runif(n_draws)
    indices   <- z > (mu_ig / (mu_ig + x))
    x[indices] <- (mu_ig^2) / x[indices]
    x / (1 + x)
    
  } else if (prior_alias %in% c("BB", "Bernoulli_beta", "bb")) {
    p_bernoulli  <- 0.25
    var_marginal <- (mu * (1 - mu)) / (phi_val + 1)
    var_cond     <- 0.125 * var_marginal
    var_between  <- var_marginal - var_cond
    delta        <- sqrt(var_between / (p_bernoulli * (1 - p_bernoulli)))
    
    mu1  <- mu + (1 - p_bernoulli) * delta
    mu0  <- mu - p_bernoulli * delta
    phi1 <- (mu1 * (1 - mu1) / var_cond) - 1
    phi0 <- (mu0 * (1 - mu0) / var_cond) - 1
    
    z  <- rbinom(n_draws, size = 1, prob = p_bernoulli)
    n1 <- sum(z == 1)
    n0 <- n_draws - n1
    
    draws <- numeric(n_draws)
    if (n1 > 0) draws[z == 1] <- rbeta(n1, shape1 = mu1 * phi1, shape2 = (1 - mu1) * phi1)
    if (n0 > 0) draws[z == 0] <- rbeta(n0, shape1 = mu0 * phi0, shape2 = (1 - mu0) * phi0)
    draws
  }
}

# 3. Build data frame with custom labels
plot_data <- map_df(dist_phi_pairs, function(pair) {
  dist_name <- pair[1]
  phi_val   <- as.numeric(pair[2])
  
  draws <- sample_prior(dist_name, phi_val, MU_TARGET, N_DRAWS)
  
  # Fetch custom title; fallback if key is missing
  panel_title <- if (!is.null(custom_specs[[dist_name]])) {
    custom_specs[[dist_name]]$Title
  } else {
    paste0(dist_name, " (\u03c6 = ", phi_val, ")")
  }
  
  data.frame(
    Distribution = dist_name,
    Phi          = phi_val,
    Panel_Title  = panel_title,
    Value        = draws
  )
})

# Preserve original matrix order in facets
plot_data$Panel_Title <- factor(plot_data$Panel_Title, levels = unique(plot_data$Panel_Title))

# 4. Extract Named Vectors for Manual Color Scales
fill_palette <- setNames(
  sapply(names(custom_specs), function(x) custom_specs[[x]]$Fill),
  sapply(names(custom_specs), function(x) custom_specs[[x]]$Title)
)

color_palette <- setNames(
  sapply(names(custom_specs), function(x) custom_specs[[x]]$Color),
  sapply(names(custom_specs), function(x) custom_specs[[x]]$Title)
)

# 5. Plot
prior_density = ggplot(plot_data, aes(x = Value, fill = Panel_Title, color = Panel_Title)) +
  geom_density(alpha = 0.7, linewidth = 0.5) +
  facet_wrap(~ Panel_Title, scales = "fixed", ncol = 3) +
  scale_x_continuous(
    labels = scales::label_scientific(digits = 1),
    breaks = function(limits) {
      # Use extended breaks for all panels
      b <- scales::breaks_extended(n = 3)(limits)
      
      # If the panel range spans past 0.01 (Logit Normal panel), drop the highest tick mark
      if (max(limits) > 0.01 && length(b) > 2) {
        b <- b[-length(b)]
      }
      return(b)
    }
  ) +
  coord_cartesian(xlim = c(0, 0.002)) +
  scale_fill_manual(values = fill_palette) +
  scale_color_manual(values = color_palette) +
  labs(
    title = "Prior Distributions.",
    subtitle = paste0("Target mean (\u03bc) = ", MU_TARGET,"."),
    x = "Probability",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey30"),
    strip.text       = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )
prior_density
ggsave('_fig/prior_density.pdf', plot = prior_density, width = 8, height = 6)
#-------------------------------------------------------------------------------
# # Plot: Overlaid Density Curves
# prior_density_overlaid = ggplot(plot_data, aes(x = Value, fill = Panel_Title, color = Panel_Title)) +
#   geom_density(alpha = 0.3, linewidth = 0.8) +
#   scale_x_continuous(labels = scales::label_scientific(digits = 1)) +
#   scale_fill_manual(values = fill_palette) +
#   scale_color_manual(values = color_palette) +
#   labs(
#     title = "Prior Distributions Comparison",
#     subtitle = paste0("Target mean (\u03bc) = ", MU_TARGET, "."),
#     x = "Probability",
#     y = "Density",
#     color = "Distribution",
#     fill = "Distribution"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title        = element_text(face = "bold", hjust = 0.5),
#     plot.subtitle     = element_text(hjust = 0.5, color = "grey30"),
#     legend.position   = "bottom",
#     panel.grid.minor  = element_blank()
#   )
# 
# prior_density_overlaid
