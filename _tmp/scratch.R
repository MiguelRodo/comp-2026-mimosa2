#-------------------------------------------------------------------------------
# Stress-test matrix: 
# This creates all possible combinations to test
stresstest_mat = expand.grid(
  Distribution = c("Beta","Uniform","Logit-Normal","X/(1+X)"),
  P            = c(20,50,100,500),
  Effect       = c(1e-10,5e-4,1e-2),
  Phi          = c(100,2000),
  stringsAsFactors = FALSE
)[1:2,]
  
# Store results:
results_summary = data.frame()

for (i in 1:nrow(stresstest_mat))
{
  dist = stresstest_mat$Distribution[i]
  p    = stresstest_mat$P[i]
  eff  = stresstest_mat$Effect[i]
  phi  = stresstest_mat$Phi[i]
  
  # Print out scenario we are running:
  cat(paste0("Running simulation ", i, "/", nrow(stresstest_mat),
             ": Distribution = ", dist, 
             ", P = ", p, 
             ", effect = ", format(eff,scientific=FALSE,drop0trailing=TRUE), 
             ", phi = ", phi, "\n"))
  
  # Simulate:
  sim = switch(dist,
               "Beta"         = simulate_MIMOSA2(effect=eff,
                                                 phi=phi,
                                                 P=p),
               "Uniform"      = simulate_MIMOSA2_uniform(effect=eff,
                                                         phi=phi,
                                                         P=p),
               "Logit-Normal" = simulate_MIMOSA2_logitnorm(effect=eff,
                                                           phi=phi,
                                                           P=p),
               "X/(1+X)"      = simulate_MIMOSA2_rexp(effect=eff,
                                                      phi=phi,
                                                      P=p)
               )
  
  # Define truth responder logic:
  # ^R: starts with R (Responders = TRUE)
  true_responder = grepl("^R",sim$truth)
  
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
  error = function(e)
    {
    fit_error <<- TRUE     # Mark run as failure 
    return(NULL)
  })
  
  # Initialise metrics:
  status     = "Success"
  iterations = NA_real_
  TPR_001    = NA_real_
  tFDR_001   = NA_real_
  tFDR_005   = NA_real_
  
  # Evaluate metrics:
  if (fit_error == FALSE && is.null(fit) == FALSE)
  {
    iterations = length(fit$inds)
    
    # Extract response calls:
    rescall_001 = getResponse(fit,threshold=0.01) # 1% FDR
    rescall_005 = getResponse(fit,threshold=0.05) # 5% FDR
    
    # Calculate sensitivity:
    if (any(true_responder))
    {
      TPR_001 = sum(rescall_001 & true_responder)/sum(true_responder)
    } else 
    {
      TPR_001 = NA
    } 
    
    # Calculate empirical FDR:
    tFDR_001 = if (sum(rescall_001)>0)
    {
      sum(rescall_001&!true_responder)/sum(rescall_001)
    } else
    {
      0
    }
    
    tFDR_005 = if (sum(rescall_005)>0)
    {
      sum(rescall_005&!true_responder)/sum(rescall_005)
    } else
    {
      0
    }
    
  } else 
  {
    status = "Optimisation crash"
  }
  
  # Store results:
  row_res = data.frame(
    Distribution = dist,
    P            = p,
    Effect       = eff,
    Phi          = phi,
    Status       = status,
    Iterations   = iterations,
    TPR_001      = TPR_001,
    tFDR_001     = tFDR_001,
    tFDR_005     = tFDR_005,
    stringsAsFactors = FALSE
  )
  
  # Combine results: 
  results_summary = rbind(results_summary,row_res)
}