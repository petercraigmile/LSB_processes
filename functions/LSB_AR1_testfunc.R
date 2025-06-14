source("functions/LSB_links.R")

LSB.AR1.sim <- function (beta, X, sigma.sq) {
  
  n <- nrow(X)

  phi <- calc.phi(X, beta)
  sigma <- sqrt(sigma.sq)
  
  y <- rep(NA, n)
  
  y[1] <- rnorm(1, 0, sqrt(sigma.sq/(1-phi[1]^2)))
  
  for (t in 2:n)
    y[t] <- rnorm(1, phi[t-1]*y[t-1], sigma)
  
  y                  
}


LSB.cond.mll2 <- function (beta, y, X) {
  
  p <- ncol(X)
  
  phi <- calc.phi(X, beta[1:p])
  
  sigma.sq <- eta.to.sigma2(beta[p+1])
  
  ts <- 2:length(y)
  
  -sum(dnorm(y[ts], phi[ts-1]*y[ts-1], sqrt(sigma.sq), log=TRUE))
}

one.sim.test <- function (beta, X.0, X.1, sigma.sq=1) {
  
  y <- LSB.AR1.sim(beta, X.1, sigma.sq) 
  
  lh_1 = -optim(
    rep(0.01, ncol(X.1) + 1),
    LSB.cond.mll2,
    y = y,
    X = X.1,
    method = "L-BFGS-B"
  )$value
  
  lh_0 = -optim(
    rep(0.01, ncol(X.0) + 1),
    LSB.cond.mll2,
    y = y,
    X = X.0,
    method = "L-BFGS-B"
  )$value
  
  2*(lh_1-lh_0)
  
}
