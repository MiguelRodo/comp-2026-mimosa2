# MIMOSA2 Simulation 
# Isabella and Tayyeb
# 26 June 2026

#-------------------------------------------------------------------------------
# Proportion of true responders: 
# True responders distributed evenly across components 1-4
# Non-responders distributed evenly across components 5-8
responders25 = c(rep(0.25 / 4, 4), rep(0.75 / 4, 4))
responders50 = c(rep(0.50 / 4, 4), rep(0.50 / 4, 4))
responders75 = c(rep(0.75 / 4, 4), rep(0.25 / 4, 4))

# Combine into a named list for tracking:
component_list = list(
  "Prop_0.25" = responders25,
  "Prop_0.50" = responders50,
  "Prop_0.75" = responders75
)

# Cell count scenarios: 
# Combine your rng cell counts into a named list for tracking:
rng_list = list(
  "Wide_High"  = c(10000, 150000),
  "Medium_Low" = c(5000, 10000),
  "Sparse"     = c(2000, 5000)
)

# Build simulation grid:
stresstest_mat = expand.grid(
  Distribution   = c("Beta"),
  Comp_Name      = names(component_list),
  P              = c(20, 50, 100),
  Effect         = c(1e-10, 5e-4, 1e-2),
  Rng_Name       = names(rng_list),
  Phi            = c(2000),
  Replication    = 1:5, 
  stringsAsFactors = FALSE
)

# Initialize objects to store results:
results_summary = data.frame()
master_obs_list = list()

#-------------------------------------------------------------------------------
# Baseline simulation
#-------------------------------------------------------------------------------
# For loop to simulation from stress test matrix: 
for (i in 1:nrow(stresstest_mat)){
  dist     = stresstest_mat$Distribution[i]
  comp_nm  = stresstest_mat$Comp_Name[i]
  p        = stresstest_mat$P[i]
  eff      = stresstest_mat$Effect[i]
  rng_nm   = stresstest_mat$Rng_Name[i]
  phi      = stresstest_mat$Phi[i]
  rep_id   = stresstest_mat$Replication[i]
  
  active_components = component_list[[comp_nm]]
  active_rng        = rng_list[[rng_nm]]
  
  # Print out scenario we are running:
  cat(paste0("Running simulation ", i, "/", nrow(stresstest_mat),
             ": Distribution = ", dist,
             ", Responder Prop = ", comp_nm,
             ", P = ", p,
             ", Effect = ", format(eff, scientific=FALSE, drop0trailing=TRUE),
             ", Rep = ", rep_id,
             ", Cell Range = ", rng_nm,
             ", phi = ", phi, "\n"))
  
  # Simulate:
  sim = simulate_MIMOSA2_alt_prior(
    effect     = eff,
    phi        = phi,
    P          = p,
    prior      = dist,
    components = active_components,
    rng        = active_rng
  )
  
  # Define truth responder logic:
  true_responder = as.numeric(sim$truth %in% c("R1","R2","R3","R4"))
  
  # Fisher's exact test (one-sided):
  fisher_p = sapply(1:p, function(j) {
    mat = matrix(c(sim$ns1[j], sim$Ntot[j]-sim$ns1[j], 
                   sim$nu1[j], sim$Ntot[j]-sim$nu1[j]), 
                 nrow = 2)
    fisher.test(mat, alternative="greater")$p.value
  })
  
  # Log-fold change (LFC):
  prop_s = sim$ns1/sim$Ntot
  prop_u = sim$nu1/sim$Ntot
  log_fold_change = log2((prop_s+1e-5)/(prop_u+1e-5))
  
  # Likelihood ratio test (LRT):
  lrt_p = sapply(1:p, function(j) { 
    fit_null = glm(cbind(c(sim$ns1[j], sim$nu1[j]),
                         c(sim$Ntot[j]-sim$ns1[j], sim$Ntot[j]-sim$nu1[j])) ~ 1, 
                   family = binomial)
    group = factor(c("S","U"))
    fit_alt  = glm(cbind(c(sim$ns1[j], sim$nu1[j]),
                         c(sim$Ntot[j]-sim$ns1[j], sim$Ntot[j]-sim$nu1[j])) ~ group, 
                   family = binomial)
    lrt_stat = 2*(logLik(fit_alt)-logLik(fit_null))
    p_val = pchisq(as.numeric(lrt_stat), df=1, lower.tail=FALSE)
    
    if (prop_s[j] < prop_u[j]) {
      p_val = 1 - (p_val/2)
    } else {
      p_val = p_val / 2
    }
    return(p_val)
  })
  
  # Fit model:
  fit_error = FALSE        # Assumes model runs successfully
  fit       = tryCatch({   # tryCatch() prevents whole simulation from crashing
    MIMOSA2(
      Ntot    = sim$Ntot,
      ns1     = sim$ns1,
      nu1     = sim$nu1,
      ns0     = sim$ns0,
      nu0     = sim$nu0,
      maxit   = 100,
      verbose = FALSE
    )
  },
  error = function(e) {
    fit_error <<- TRUE     # Mark run as failure
    return(NULL)
  })
  
  # Initialize point metrics:
  status     = "Success"
  iterations = NA_real_
  TPR_001    = NA_real_
  tFDR_001   = NA_real_
  tFDR_005   = NA_real_
  
  # Evaluate metrics:
  if (fit_error == FALSE && is.null(fit) == FALSE) {
    iterations = length(fit$inds)
    mimosa_prob = rowSums(fit$z[, 1:4, drop = FALSE])
    
    # Extract response calls:
    rescall_001 = getResponse(fit, threshold=0.01) # 1% FDR
    rescall_005 = getResponse(fit, threshold=0.05) # 5% FDR
    
    # Calculate sensitivity:
    if (any(true_responder==1)) {
      TPR_001 = sum(rescall_001 & true_responder==1) / sum(true_responder==1)
    } else {
      TPR_001 = NA
    } 
    
    # Calculate empirical FDR:
    tFDR_001 = if (sum(rescall_001) > 0) {
      sum(rescall_001 & true_responder!=1) / sum(rescall_001)
    } else {
      0
    }
    
    tFDR_005 = if (sum(rescall_005) > 0) {
      sum(rescall_005 & true_responder!=1) / sum(rescall_005)
    } else {
      0
    }
    
  } else {
    status = "Optimisation crash"
    mimosa_prob = rep(NA_real_, p)
  }
  
  # Store results:
  row_res = data.frame(
    Distribution = dist,
    Res_prop     = comp_nm,
    P            = p,
    Effect       = eff,
    Cell_range   = rng_nm,
    Phi          = phi,
    Replication  = rep_id,
    Status       = status,
    Iterations   = iterations,
    TPR_001      = TPR_001,
    tFDR_001     = tFDR_001,
    tFDR_005     = tFDR_005,
    stringsAsFactors = FALSE
  )
  
  # Combine results:
  results_summary = rbind(results_summary, row_res)
  
  # Store results:
  scenario_obs = data.frame(
    Distribution = dist,
    Res_prop     = comp_nm,
    P            = p,
    Effect       = eff,
    Cell_range   = rng_nm,
    Phi          = phi,
    Replication  = rep_id,
    Subject_id   = 1:p,
    Truth        = true_responder,
    MIMOSA2_prob = mimosa_prob,
    Fisher_p     = fisher_p,
    LRT_P        = lrt_p,
    Log2_FC      = log_fold_change
  )
  
  # Store results:
  master_obs_list[[i]] = scenario_obs
}

# Combine:
results_continuous = do.call(rbind, master_obs_list)

# Save results: 
save(results_summary, 
     results_continuous,
     file = '_simulations/Simulation_1.1.Rdata')