

## ======================================================================
## Script: sims5_Paparoditis.R
##
## Purpose: A simulation to estimate the power of the Paparoditis
## (2010) test for stationarity in the case of a LSB-AR(1) process, as
## the sample length T (N in this code) and slope beta_{11} is varied.
## ======================================================================


library(spectral)
library(parallel)

source("functions/LSB_links.R")
source("functions/basis_functions.R")
source("functions/LSB_AR1.R")

source("functions/local_pgram.R")
source("functions/lag_window_spectra.R")
source("functions/Paparoditis_statistics.R")
source("functions/load_save.R")




Papa.sim <- function (N, m, beta1, REPS=10000, CORES=6, use.Z=FALSE) {

    us <- u.range(N, m, N-m)
    
    the.beta <- c(0.5, beta1)
    
    X.phi <- poly.basis(N, 1)

    phi <- calc.phi(X.phi, the.beta)

    h <- 1.5 * N^(-1/5)
    b <- (m/N)^(1/5) * h

    if (!use.Z) {
        
        sims <- sim.Q.stats(us, m, h=h, b=b,
                            function (x) LSB.AR1.sim(N, phi, 1),
                            REPS=REPS, CORES=CORES)
    } else {

        sims <- sim.Q.stats(us, m, h=h, b=b,
                            function (x) rnorm(N),
                            REPS=REPS, CORES=CORES)
    }
    
    list(N=N, m=m, beta1=beta1, the.beta=the.beta, phi=phi,
         sims=sims, h=h, b=b, us=us)
}



beta1s <- seq(0, 1.5, 0.3)

Q300 <- lapply(beta1s, function (b1) Papa.sim(300, 30, b1, REPS=2000))
save.it("Q300")

Q300.0 <- Papa.sim(300, 30, 0, REPS=10000)
save.it("Q300.0")

Q300.Z <- Papa.sim(300, 30, 0, REPS=2000, use.Z=TRUE)
save.it("Q300.Z")



beta1s <- seq(0, 1.5, 0.3)

Q500 <- lapply(beta1s, function (b1) Papa.sim(500, 40, b1, REPS=2000))
save.it("Q500")

Q500.0 <- Papa.sim(500, 40, 0, REPS=10000)
save.it("Q500.0")

Q500.Z <- Papa.sim(500, 40, 0, REPS=2000, use.Z=TRUE)
save.it("Q500.Z")



beta1s <- seq(0, 1.5, 0.3)

Q1000 <- lapply(beta1s, function (b1) Papa.sim(1000, 64, b1, REPS=10000))
save.it("Q1000")

Q1000.0 <- Papa.sim(1000, 64, 0, REPS=10000)
save.it("Q1000.0")

Q1000.Z <- Papa.sim(1000, 64, 0, REPS=2000, use.Z=TRUE)
save.it("Q1000.Z")



beta1s <- seq(0, 1.5, 0.3)

Q2000 <- lapply(beta1s, function (b1) Papa.sim(2000, 96, b1, REPS=10000))
save.it("Q2000")

Q2000.0 <- Papa.sim(2000, 96, 0, REPS=10000)
save.it("Q2000.0")

Q2000.Z <- Papa.sim(2000, 96, 0, REPS=2000, use.Z=TRUE)
save.it("Q2000.Z")

