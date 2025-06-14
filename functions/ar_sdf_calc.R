LSB.ARp.sdf <- function(par.sim, X){
  require(spectral)
  the.beta = the.phi = list()
  
  the.beta[[1]] <- par.sim[1:ncol(as.matrix(X[[1]]))]
  for(i in 2:(p+1)){
    the.beta[[i]] <- 
      par.sim[(sum(sapply(1:(i-1), 
                          function(i) 
                            ncol(as.matrix(X[[i]])))) + 1) : 
                (sum(sapply(1:i, 
                            function(i) ncol(as.matrix(X[[i]])))))]
  }
  
  for(i in 1:p){
    the.phi[[i]] <- matrix(NA, i, (T.-i))
    the.phi[[i]][i,] <- 
      h(drop(as.matrix(X[[i]][-seq(1, i, 1),]) %*% the.beta[[i]]))
  }
  
  for(i in 2:p){
    for(j in 1: (i-1)){
      the.phi[[i]][j,] <- 
        the.phi[[i-1]][j,-1] - the.phi[[i]][i,]*the.phi[[i-1]][i-j,-1]
    }
  }
  
  the.sigma.sq <- exp(drop(X[[p+1]] %*% the.beta[[p+1]]))
  
  fs <- seq(0, 0.5, len = 256)
  
  ar.sdfs <- 
    lapply(1:(T.-p), 
           function(k) dB(ar.sdf(freqs = fs, 
                                 phi = as.vector(the.phi[[p]][1:p,k]), 
                                 sigma2 = the.sigma.sq[k])))
  
  t(sapply(ar.sdfs, function (x) x))
  
  
}


LSB.AR2.log.sdf <- function(phi.21, phi.22, sigma.sq, fs){
  ## ======================================================================
  ## Calculate log(SDF) from AR2 parameter curves phi21 and phi22
  ## for fixed sigma.sq and frequencies fs
  ## ======================================================================
  
  log.ar.sdfs <- 
    lapply(1:length(phi.21), 
           function(k) log(ar.sdf(freqs = fs, 
                                  phi = c(phi.21[k], phi.22[k]), 
                                  sigma2 = sigma.sq)))
  
  t(sapply(log.ar.sdfs, function (x) x))
  
}



LSB.AR2.phi.log.sdf <- function(X, theta, fs){
  ## ======================================================================
  ## Calculate log(SDF) and AR2 parameter curves phi21 and phi22
  ## from LSB basis parameters 
  ## ======================================================================
  
  n = nrow(X)
  
  beta.tt1 <- theta[1:ncol(X)]
  beta.tt2 <- theta[(ncol(X) + 1):(2 * (ncol(X)))]
  beta.tt3 <- theta[length(theta)]
  
  
  phi.11 <- eta.to.phi(calc.eta(X[-1, ], beta.tt1))
  phi.22 <- eta.to.phi(calc.eta(X[-c(1, 2), ], beta.tt2))
  sigma.sq <- exp(beta.tt3) ^ 2
  
  ts <- 3:n
  
  phi.21 <- array(NA, length(phi.22))
  
  phi.21[ts - 2] <- phi.11[ts - 2] - phi.22[ts - 2] * phi.11[ts - 2]
  
  log.ar.sdfs <- 
    lapply(1:(n-2), 
           function(k) log(ar.sdf(freqs = fs, 
                                 phi = c(phi.21[k], phi.22[k]), 
                                 sigma2 = sigma.sq)))
  
  the.log.ar.sdfs <- t(sapply(log.ar.sdfs, function (x) x))
  
  list(phi.21 = phi.21,
       phi.22 = phi.22,
       log.lsb2.sdf = the.log.ar.sdfs)
  
}
