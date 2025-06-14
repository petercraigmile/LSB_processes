
## ======================================================================
## Script: sims1_LSB_AR2_sims.R
## Purpose: Carry out simulations for the verification of asymptotic
## properties of LSB-AR(2) models
## ======================================================================


source("functions/LSB_AR2_sim_functions.R")
source("functions/sim_stats.R")


# Basis parameters for LSB-AR(2)
true.beta.tt1 <- c(0.61, 1.71,-1.27)
true.beta.tt2 <- c(-3.52, 5.5,-3)
true.sigma.sq <- 1
true.beta.tt3 <- log(sqrt(true.sigma.sq))

true.param <- c(true.beta.tt1, true.beta.tt2, true.beta.tt3)


REPS <- 500


ns <- c(128, 256, 512, 1024, 2048, 4096)

the.X <- sapply(ns, function(n) poly.basis.for.sim(n))

asymp.sim.data.exact <-
  lapply(1:length(ns), function(n)
    LSB.AR2.est.combine.reps(
      n = ns[n],
      beta.tt1 = true.beta.tt1,
      beta.tt2 = true.beta.tt2,
      sigma.sq = true.sigma.sq,
      X = the.X[[n]],
      lik.method = "exact",
      REPS = REPS
    ))


asymp.sim.data.whittle <-
  lapply(1:length(ns), function(n)
    LSB.AR2.est.combine.reps(
      n = ns[n],
      beta.tt1 = true.beta.tt1,
      beta.tt2 = true.beta.tt2,
      sigma.sq = true.sigma.sq,
      X = the.X[[n]],
      lik.method = "whittle",
      REPS = REPS
    ))

save(asymp.sim.data.exact, 
     file = "results/sims_exact_hermite_poly_sigsq.RData")
save(asymp.sim.data.whittle, 
     file = "results/sims_whittle_hermite_poly_sigsq.RData")
