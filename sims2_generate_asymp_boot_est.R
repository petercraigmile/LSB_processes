
## ======================================================================
## Script: sims2_generate_asymp_boot_est.R
## Purpose: Calculate estimates and bootstrapped CIs for the LSB-AR(2)
## simulations
## ======================================================================

source("functions/LSB_AR2_sim_functions.R")
source("functions/sim_stats.R")

load("results/sims_exact_hermite_poly_sigsq.RData")
load("results/sims_whittle_hermite_poly_sigsq.RData")

# Basis parameters for LSB-AR(2)
true.beta.tt1 <- c(0.61, 1.71,-1.27)
true.beta.tt2 <- c(-3.52, 5.5,-3)
true.sigma.sq <- 1
true.beta.tt3 <- log(sqrt(true.sigma.sq))

true.param <- c(true.beta.tt1, true.beta.tt2, true.beta.tt3)

REPS <- 500

ns <- c(128, 256, 512, 1024, 2048, 4096)

the.X <- sapply(ns, function(n) poly.basis.for.sim(n))

# Calculate stats from the sims 
exact.stats <- lapply(1:length(ns), function(n) 
  calc.stats(ns[n], asymp.sim.data.exact[[n]], true.param))

whittle.stats <- lapply(1:length(ns), function(n) 
  calc.stats(ns[n], asymp.sim.data.whittle[[n]], true.param))

# Number of bootstrapped samples
boot.reps <- 200

# Get bootstrapped samples 
bootsamp.exact = bootsamp.whittle <- list()

bootsamp.exact <- 
  lapply(1:length(ns),
         function(x) {
           replicate(boot.reps,
                     asymp.sim.data.exact[[x]][, sample(1:REPS, 
                                                             replace = TRUE)])
         })


bootsamp.whittle <- 
  lapply(1:length(ns),
         function(x) {
           replicate(boot.reps,
                     asymp.sim.data.whittle[[x]][, sample(1:REPS, 
                                                               replace = TRUE)])
         })

# Calculate bootstrap statistics

boot.stat.exact = boot.stat.whittle <- list()

for(i in 1:length(ns)) {
  boot.stat.exact[[i]] <-
    lapply(1:boot.reps, function(x)
      calc.stats(n = ns[i],
                 samp = bootsamp.exact[[i]][, , x], 
                 theta = true.param))
  
  boot.stat.whittle[[i]] <-
    lapply(1:boot.reps, function(x)
      calc.stats(n = ns[i],
                 samp = bootsamp.whittle[[i]][, , x],
                 theta = true.param))
  
}

# Calculate bootstrapped CIs

boot.ci.exact <- boot.ci(asymp.stats = exact.stats,
                         boot.stats = boot.stat.exact)


boot.ci.whittle <- boot.ci(asymp.stats = whittle.stats,
                           boot.stats = boot.stat.whittle)

save(true.param, ns, boot.ci.exact, boot.ci.whittle,
     file = "results/boot_cis.RData")

