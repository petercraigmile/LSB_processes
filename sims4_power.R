## ======================================================================
## Script: sim4_power.R
## Purpose: Calculate power for test of stationarity described in 
## Section 4 for and AR(1) process
## ======================================================================

source("functions/LSB_links.R")
source("functions/basis_functions.R")
source("functions/LSB_AR1_testfunc.R")


REPS <- 10000
beta1 <- seq(0, 1.5, 0.3)

test.n.beta <- function(n, beta, REPS){
  
  us <- (0:(n-1))/n
  
  X.0 <- poly.basis(n, deg = 0)
  
  X.1 <-  poly.basis(n, deg = 1, 
                     orthogonal = FALSE,
                     normalized = FALSE)
  
  the.beta <- c(0.5, beta)
  
  sims <- pbmcapply::pbmclapply(1:REPS, function (the.rep) {
    one.sim.test(the.beta, X.0, X.1)
  }, mc.cores = parallel::detectCores() - 1)
  
  T.stat <- unlist(sims)
  
  sum(T.stat > qchisq(0.95, 1))/REPS
  
}

sim.300 <- sapply(beta1, function(b) test.n.beta(300, b, REPS))
sim.500 <- sapply(beta1, function(b) test.n.beta(500, b, REPS))
sim.1000 <- sapply(beta1, function(b) test.n.beta(1000, b, REPS))
sim.2000 <- sapply(beta1, function(b) test.n.beta(2000, b, REPS))

