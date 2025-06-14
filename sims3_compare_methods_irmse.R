
## ======================================================================
## Script: sims3_compare_methods_irmse.R
## Purpose: Compare IRMSEs for the time varying AR curves and the
## time varying SDF for our methods vs LSTS package based on
## Palma, Olea (2010)
## ======================================================================

source("functions/LSB_AR2_sim_functions.R")
source("functions/sim_stats.R")


# Basis parameters for LSB-AR(2)
true.beta.tt1 <- c(0.61, 1.71,-1.27)
true.beta.tt2 <- c(-3.52, 5.5,-3)
true.sigma.sq <- 1
true.beta.tt3 <- log(sqrt(true.sigma.sq))

true.param <- c(true.beta.tt1, true.beta.tt2, true.beta.tt3)

# Number of sims
REPS <- 500

# List of time series lengths
ns <- c(128, 256, 512, 1024, 2048, 4096)

# Build orthogonal polynomial basis for our methods
the.X <- sapply(ns, function(n) poly.basis.for.sim(n))

# Build raw polynomial basis for Palma Olea method
the.X.lsts <- sapply(ns, function(n) poly.basis(n, deg = 2, 
                                                orthogonal = FALSE,
                                                normalized = FALSE))
# True AR curves
true.phis <- lapply(1:length(ns),
                    function(n) LSB.AR2.PARtoAR(true.param, the.X[[n]]))

# Frequencies
fs <- seq(0, 0.5, len = 100)

# Get simulation data

ar2.method.comp.sims <-
  lapply(1:length(ns), function(n)
    LSB.AR2.est.comp.combine.reps(
      n = ns[n],
      beta.tt1 = true.beta.tt1,
      beta.tt2 = true.beta.tt2,
      sigma.sq = true.sigma.sq,
      X = the.X[[n]],
      X.lsts = the.X.lsts[[n]],
      REPS = REPS
    ))


save(ar2.method.comp.sims, 
     file = "results/AR2.est.method.comp.sims.lsts.RData")

load("results/AR2.est.method.comp.sims.lsts.RData")

# Calculate IRMSE stats for each method

irmse.stats.exact <-
  lapply(1:length(ns), function(n)
    calc.sim.irmses(
      true.param = true.param,
      sims = ar2.method.comp.sims[[n]]$exact.ests,
      X = the.X[[n]],
      fs = fs,
      REPS = REPS
    ))

irmse.stats.whittle <-
  lapply(1:length(ns), function(n)
    calc.sim.irmses(
      true.param = true.param,
      sims = ar2.method.comp.sims[[n]]$whittle.ests,
      X = the.X[[n]],
      fs = fs,
      REPS = REPS
    ))

irmse.stats.lsts <-
  lapply(1:length(ns), function(n)
    calc.sim.irmses(
      true.param = true.param,
      sims = ar2.method.comp.sims[[n]]$lsts.ests,
      X = the.X[[n]],
      X.lsts = the.X.lsts[[n]],
      fs = fs,
      REPS = REPS,
      lsts = TRUE
    ))

# IRMSEs for exact likelihood
phi.21.irmse.exact <- 
  sapply(1:length(ns), 
         function(x) irmse.stats.exact[[x]]$IRMSEs[, "IRMSE.phi.21"])

phi.22.irmse.exact <- 
  sapply(1:length(ns), 
         function(x) irmse.stats.exact[[x]]$IRMSEs[, "IRMSE.phi.22"])

log.sdf.irmse.exact <- 
  sapply(1:length(ns), 
         function(x) irmse.stats.exact[[x]]$IRMSEs[, "IRMSE.log.sdf"])

# IRMSEs for block Whittle likelihood

phi.21.irmse.whittle <- 
  sapply(1:length(ns), 
         function(x) irmse.stats.whittle[[x]]$IRMSEs[, "IRMSE.phi.21"])

phi.22.irmse.whittle <- 
  sapply(1:length(ns), 
         function(x) irmse.stats.whittle[[x]]$IRMSEs[, "IRMSE.phi.22"])

log.sdf.irmse.whittle <- 
  sapply(1:length(ns), 
         function(x) irmse.stats.whittle[[x]]$IRMSEs[, "IRMSE.log.sdf"])

# IRMSEs for Palma Olea method using LSTS package 

phi.21.irmse.lsts <- 
  sapply(1:length(ns), 
         function(x) irmse.stats.lsts[[x]]$IRMSEs[, "IRMSE.phi.21"])

phi.22.irmse.lsts <- 
  sapply(1:length(ns), 
         function(x) irmse.stats.lsts[[x]]$IRMSEs[, "IRMSE.phi.22"])

log.sdf.irmse.lsts <- 
  sapply(1:length(ns), 
         function(x) irmse.stats.lsts[[x]]$IRMSEs[, "IRMSE.log.sdf"])

save(phi.21.irmse.exact, phi.21.irmse.whittle, phi.21.irmse.lsts,
     phi.22.irmse.exact, phi.22.irmse.whittle, phi.22.irmse.lsts,
     log.sdf.irmse.exact, log.sdf.irmse.whittle, log.sdf.irmse.lsts,
     file = "results/IRMSE_data.RData")


