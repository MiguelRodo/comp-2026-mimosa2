setwd("C:/GitHub/comp-2026-mimosa2/mimosa2_sim_v3")

#-------------------------------------------------------------------------------
# Effect size heterogeneity 
#-------------------------------------------------------------------------------
# checks
length(list.files("_simulations/heterogeneity/tasks"))      # should be 40500
file.exists("_simulations/heterogeneity/design_full.rds")   # should be TRUE

# combine
Sys.setenv(SIM_PROFILE = "extended")
system2(file.path(R.home("bin"), "Rscript"),
        c("sims/29_combine.R", "heterogeneity"))

list.files("_simulations/heterogeneity", pattern = "results")
res <- readRDS("_simulations/heterogeneity/heterogeneity_results_extended.rds")
length(unique(res$design$Scenario_ID))   # should be 405
table(res$fits$Status)
res$run_info

#-------------------------------------------------------------------------------
# Cell count imbalance
#-------------------------------------------------------------------------------
length(list.files("_simulations/imbalance/tasks"))      # should be 24000
file.exists("_simulations/imbalance/design_full.rds")   # should be TRUE
list.files("_simulations/imbalance", pattern = "results")  # see what already exists

Sys.setenv(SIM_PROFILE = "extended")
system2(file.path(R.home("bin"), "Rscript"),
        c("sims/29_combine.R", "imbalance"))

res <- readRDS("_simulations/imbalance/imbalance_results_extended.rds")
length(unique(res$design$Scenario_ID))   # should be 120
table(res$fits$Status)
sapply(res$run_info, function(x) x$started)   # when each chunk ran
any(sapply(res$run_info, function(x) !is.null(x$did_cc)))  # TRUE = DiD fix included

#-------------------------------------------------------------------------------
# Prior
#-------------------------------------------------------------------------------
length(list.files("_simulations/prior/tasks"))         # should be 60000
file.exists("_simulations/prior/design_full.rds")      # should be TRUE
list.files("_simulations/prior", pattern = "results")  # see what already exists

Sys.setenv(SIM_PROFILE = "extended")
system2(file.path(R.home("bin"), "Rscript"),
        c("sims/29_combine.R", "prior"))

res <- readRDS("_simulations/prior/prior_results_extended.rds")
length(unique(res$design$Scenario_ID))   # should be 300
table(res$fits$Status)
any(sapply(res$run_info, function(x) !is.null(x$did_cc)))  # FALSE = needs DiD fix
#-------------------------------------------------------------------------------
