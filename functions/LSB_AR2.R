
library(spectral)

source("functions/tapered_pgram.R")
source("functions/LSB_links.R")

LSB.AR2.sim <- function(beta.tt1, beta.tt2, sigma.sq, x){
  
  ## ======================================================================
  ## Simulate LSB.AR2 processes with time-varying partial AR parameter
  ## 'phi's and variance 'sigma2'
  ## ======================================================================
  n <- nrow(x)
  
  # time varying partial autocorrelation of order 1
  phi.11 <- eta.to.phi(calc.eta(x[-1, ], beta.tt1))
  
  # time varying partial autocorrelation of order 2
  phi.22 <- eta.to.phi(calc.eta(x[-c(1, 2), ], beta.tt2))
  
  y <- array(NA, n)
  
  y[1] <-
    rnorm(1, 0, sqrt(sigma.sq / ((1 - phi.11[1] ^ 2) * (1 - phi.22[1] ^ 2))))
  y[2] <-
    rnorm(1, phi.11[1] * y[1], sqrt(sigma.sq / (1 - phi.22[1] ^ 2)))
  
  ts <- 3:length(y)
  
  phi.21 <- array(NA, length(phi.22))
  phi.21[ts - 2] <- phi.11[ts - 2] - phi.22[ts - 2] * phi.11[ts - 2]
  
  for (t in ts) {
    
    y[t] <-
      rnorm(1, phi.21[t - 2] * y[t - 1] + phi.22[t - 2] * y[t - 2], 
            sqrt(sigma.sq))
    
  }
  
  as.vector(y)
  
}

LSB.AR2.sim.exp <- function(beta.tt1, beta.tt2, sigma.sq, x){
  
  ## ======================================================================
  ## Simulate LSB.AR2 processes with time-varying partial AR parameter
  ## 'phi's and variance 'sigma2'
  ## ======================================================================
  n <- nrow(x)
  
  # time varying partial autocorrelation of order 1
  phi.11 <- eta.to.phi(calc.eta(x[-1, ], beta.tt1))
  
  # time varying partial autocorrelation of order 2
  phi.22 <- eta.to.phi(calc.eta(x[-c(1, 2), ], beta.tt2))
  
  y <- array(NA, n)
  
  ts <- 3:length(y)
  
  sd1 <- sqrt(sigma.sq / ((1 - phi.11[1] ^ 2) * (1 - phi.22[1] ^ 2)))
  sd2 <- sqrt(sigma.sq / (1 - phi.22[1] ^ 2))
  
  y[1] <- (rexp(1) - 1) * sd1
  y[2] <- phi.11[1] * y[1] + (rexp(1) - 1) * sd2
  
  phi.21 <- array(NA, length(phi.22))
  phi.21[ts - 2] <- phi.11[ts - 2] - phi.22[ts - 2] * phi.11[ts - 2]
  
  for (t in ts) {
    
    y[t] <- phi.21[t - 2] * y[t - 1] + phi.22[t - 2] * y[t - 2] + 
      (rexp(1) - 1) * sqrt(sigma.sq)
    
  }
  
  as.vector(y)
  
}



LSB.AR2.exact.cond.lik <- function(y, X, theta){
  ## ======================================================================
  ## Calculate the exact conditional likelihood for LSB-AR(2) process with 
  ## time invariant sigma.sq
  ## ======================================================================
  
  beta.tt1 <- theta[1:ncol(X)]
  beta.tt2 <- theta[(ncol(X) + 1):(2 * (ncol(X)))]
  beta.tt3 <- theta[length(theta)]
  
  
  phi.11 <- eta.to.phi(calc.eta(X[-1, ], beta.tt1))
  phi.22 <- eta.to.phi(calc.eta(X[-c(1, 2), ], beta.tt2))
  sigma.sq <- eta.to.sigma2(beta.tt3)
  
  ts <- 3:length(y)
  
  phi.21 <- array(NA, length(phi.22))
  
  phi.21[ts - 2] <- phi.11[ts - 2] - phi.22[ts - 2] * phi.11[ts - 2]
  
  
  
  # minus log-likelihood
  
  - dnorm(y[1], 0, 
          sqrt(sigma.sq / ((1 - phi.11[1] ^ 2) * (1 - phi.22[1] ^ 2))),
          log = T)  -
    dnorm(y[2], phi.11[1] * y[1],
          sqrt(sigma.sq / (1 - phi.22[1] ^ 2)), log = T) -
    sum(dnorm(y[ts], phi.21[ts - 2] * y[ts - 1] + phi.22[ts - 2] * y[ts - 2], 
              sqrt(sigma.sq), log = T))
}



LSB.AR2.whit.lik.setup <- function (y, N=NULL, S=NULL) {
  ## ======================================================================
  ## Setup code to calculate Whittle likelihood using tapered pgram with
  ## block length N and step size S. 
  ## ======================================================================
  
  n <- length(y)
  
  if (is.null(N)) {
    
    N <- trunc(n ^ (0.6))
    
  }
  
  if (is.null(S)) {
    
    S <- trunc(N * 0.35)
    
  }
  
  within <- 0:(N - 1)
  
  ks <- seq(1, n - N, S)
  
  inds <- lapply(ks, function (j) j + within)
  
  mids <- ks + trunc(S / 2)
  
  pgs <- lapply(inds, function (is) tap.prdgram(y[is]))
  
  exp.terms.1 <- lapply(pgs, function (pg) exp(-1i * 2 * pi * pg$freq))
  exp.terms.2 <- lapply(pgs, function (pg) exp(-2i * 2 * pi * pg$freq))
  
  specs <- lapply(pgs, function (pg) pg$spec)
  
  return(
    list(
      y = y,
      n = n,
      N = N,
      S = S,
      inds = inds,
      mids = mids,
      exp.terms.1 = exp.terms.1,
      exp.terms.2 = exp.terms.2,
      specs = specs
    )
  )
}


LSB.AR2.whit.lik <- function(theta, X, obj) {
  ## ======================================================================
  ## Calculate Whittle likelihood for LSB-AR(2) process with 
  ## time invariant sigma.sq using LSB.AR2.whit.lik.setup() object
  ## ======================================================================
  
  beta.tt1 <- theta[1:ncol(X)]
  beta.tt2 <- theta[(ncol(X) + 1):(2 * (ncol(X)))]
  beta.tt3 <- theta[length(theta)]
  sigma.sq <- eta.to.sigma2(beta.tt3)
  
  
  all.phi.11 <- eta.to.phi(calc.eta(X[-1, ], beta.tt1))
  all.phi.22 <- eta.to.phi(calc.eta(X[-c(1, 2), ], beta.tt2))
  
  ts <- 3:nrow(X)
  
  all.phi.21 <- array(NA, length(all.phi.22))
  
  all.phi.21[ts - 2] <-
    all.phi.11[ts - 2] - all.phi.22[ts - 2] * all.phi.11[ts - 2]
  
  
  phi.11 <- all.phi.11[obj$mids + 1]
  phi.22 <- all.phi.22[obj$mids + 2]
  phi.21 <- all.phi.21[obj$mids + 2]
  
  mean(sapply(1:length(phi.11), function (k) {
    sdf.k = sigma.sq / (2 * pi * Mod(1 - phi.21[k] * obj$exp.terms.1[[k]] -
                                       phi.22[k] * obj$exp.terms.2[[k]]) ^ 2)
    
    (sum(log(sdf.k) + obj$specs[[k]] / sdf.k)) / (4 * pi)
    
  })) 
}

LSB.AR2.PARtoAR <- function(theta, X){
  ## ======================================================================
  ## Calculate AR2 parameter curves phi21 and phi22 from LSB basis parameters
  ## ======================================================================
  
  n <- nrow(X)
  beta.tt1 <- theta[1:ncol(X)]
  beta.tt2 <- theta[(ncol(X) + 1):(2 * (ncol(X)))]
  
  # time varying partial autocorrelation of order 1
  phi.11 <- eta.to.phi(calc.eta(X[-1,], beta.tt1))
  
  # time varying partial autocorrelation of order 2
  phi.22 <- eta.to.phi(calc.eta(X[-c(1,2),], beta.tt2))
  
  ts <- 3:n
  
  phi.21 <- array(NA, length(phi.22))
  phi.21[ts-2] <- phi.11[ts-2] - phi.22[ts-2]*phi.11[ts-2]
  
  list(phi.21 = phi.21,
       phi.22 = phi.22)
}

LSTS.AR2.beta.to.phi <- function(theta, X){
  ## ======================================================================
  ## Calculate AR2 parameter curves phi21 and phi22 from LSB basis parameters
  ## ======================================================================
  
  n <- nrow(X)
  beta.tt1 <- theta[1:ncol(X)]
  beta.tt2 <- theta[(ncol(X) + 1):(2 * (ncol(X)))]
  
  # time varying partial autocorrelation of order 1
  phi.21 <- calc.eta(X[-c(1,2),], beta.tt1)
  
  # time varying partial autocorrelation of order 2
  phi.22 <- calc.eta(X[-c(1,2),], beta.tt2)
  
  list(phi.21 = phi.21,
       phi.22 = phi.22)
}

