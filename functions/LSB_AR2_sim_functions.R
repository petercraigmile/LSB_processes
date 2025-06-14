
source("functions/LSB_AR2.R")
source("functions/basis_functions.R")
source("functions/LSB_AR_SDF.R")

LSB.AR2.exact.est.one.rep <- function (beta.tt1, beta.tt2, X, sigma.sq,
                                       exp.gen = FALSE) {
  ## ======================================================================
  ## Simulate single LSB-AR(2) process and estimate the basis parameters
  ## using exact likelihood method
  ## ======================================================================
  
  if(exp.gen){
    y <- LSB.AR2.sim.exp(beta.tt1, beta.tt2, sigma.sq, X)
  } else{
    y <- LSB.AR2.sim(beta.tt1, beta.tt2, sigma.sq, X)
  }
  
  optim(
    rep(1, 2 * ncol(X) + 1),
    LSB.AR2.exact.cond.lik,
    y = y,
    X = X,
    method = "BFGS"
  )$par
}

LSB.AR2.whittle.est.one.rep <- function (beta.tt1, beta.tt2, X, sigma.sq, 
                                         N, S, exp.gen = FALSE) {
  ## ======================================================================
  ## Simulate single LSB-AR(2) process and estimate the basis parameters
  ## using whittle likelihood method
  ## ======================================================================
  
  if(exp.gen){
    y <- LSB.AR2.sim.exp(beta.tt1, beta.tt2, sigma.sq, X)
  } else{
    y <- LSB.AR2.sim(beta.tt1, beta.tt2, sigma.sq, X)
  }
  
  obj <- LSB.AR2.whit.lik.setup(y, N, S)
  
  optim(
    rep(1, 2 * ncol(X) + 1),
    LSB.AR2.whit.lik,
    obj = obj,
    X = X,
    method = "BFGS"
  )$par
}

LSB.AR2.method.comp.one.rep <- function (beta.tt1, beta.tt2, 
                                         X, X.lsts,
                                         sigma.sq, N, S) {
  ## ======================================================================
  ## Simulate single LSB-AR(2) process and estimate the basis parameters
  ## using exact likelihood, whittle likelihood and LSTS package
  ## ======================================================================
  
  require(LSTS)
  
  y <- LSB.AR2.sim(beta.tt1, beta.tt2, sigma.sq, X)
  
  
  lsbar2.exact.est <- optim(
    rep(0.01, 2 * ncol(X) + 1),
    LSB.AR2.exact.cond.lik,
    y = y,
    X = X,
    method = "BFGS"
  )$par
  
  obj <- LSB.AR2.whit.lik.setup(y, N=NULL, S=NULL)
  
  lsbar2.whittle.est <- optim(
    rep(0.01, 2 * ncol(X) + 1),
    LSB.AR2.whit.lik,
    obj = obj,
    X = X,
    method = "BFGS"
  )$par
  
  lsts.ar2.est <- LS.whittle(
    series = y, start = rep(0.01, 2 * ncol(X.lsts) + 1), 
    order = c(p = 2, q = 0),
    ar.order = c(ncol(X.lsts)-1, ncol(X.lsts)-1), sd.order = 0
  )$coef
  
  list(lsbar2.exact.est = lsbar2.exact.est,
       lsbar2.whittle.est = lsbar2.whittle.est,
       lsts.ar2.est = lsts.ar2.est)
}


LSB.AR2.est.combine.reps <- function (n, beta.tt1, beta.tt2, sigma.sq, 
                                      X, lik.method = "whittle",
                                      N = NULL, S = NULL,
                                      REPS,
                                      exp.gen = FALSE,
                                      use.parallel = TRUE) {
  
  ## ======================================================================
  ## Calculate beta estimates for REPS replications for LSB-AR(2)
  ## simulations of length n
  ## ======================================================================
  
  if (lik.method == "whittle") {
    if (is.null(N)) {
      N = trunc(n ^ 0.6)
    }
    if (is.null(S)) {
      S = trunc(0.35 * N)
    }
  }
  
  
  if (use.parallel) {
    
    require(parallel)
  
    require(pbmcapply)
    
    sims <- pbmcapply::pbmclapply(1:REPS, function (the.rep) {
      if (lik.method == "whittle") {
        
        LSB.AR2.whittle.est.one.rep(beta.tt1, beta.tt2, X, 
                                    sigma.sq, N, S, exp.gen)
        
      } else if (lik.method == "exact") {
        
        LSB.AR2.exact.est.one.rep(beta.tt1, beta.tt2, X,
                                  sigma.sq, exp.gen)
        
      } else{
        
        stop("lik.method has to be either whittle or exact")
        
      }
    }, mc.cores = parallel::detectCores() - 2)
    
  } else {
    sims <-
      lapply(1:REPS, function (the.rep) {
        
        if (the.rep %% 250 == 0)
          cat(the.rep, " ")
        
        if (lik.method == "whittle") {
          
          LSB.AR2.whittle.est.one.rep(beta.tt1, beta.tt2, X,
                                      sigma.sq, N, S, exp.gen)
          
        } else if (lik.method == "exact") {
          
          LSB.AR2.exact.est.one.rep(beta.tt1, beta.tt2, X, 
                                    sigma.sq, exp.gen)
          
        } else{
          
          stop("lik.method has to be either whittle or exact")
          
        }
        
      })
  }
  
  sapply(sims, function (x) x)
  
}


LSB.AR2.est.comp.combine.reps <- function (n, beta.tt1, beta.tt2, sigma.sq, 
                                           X, X.lsts, N = NULL, S = NULL, REPS,
                                           use.parallel = TRUE) {
  
  ## ======================================================================
  ## Calculate beta estimates for REPS replications for LSB-AR(2)
  ## simulations of length n
  ## ======================================================================
  

  if (is.null(N)) {
    N = trunc(n ^ 0.6)
  }
  if (is.null(S)) {
    S = trunc(0.35 * N)
  }
  
  if (use.parallel) {
    
    require(parallel)
    
    require(pbmcapply)
    
    sims <- pbmcapply::pbmclapply(1:REPS, function (the.rep) {
        
      LSB.AR2.method.comp.one.rep(beta.tt1, beta.tt2, X, X.lsts, sigma.sq, N, S)

    }, mc.cores = parallel::detectCores() - 2)
    
  } else {
    
    sims <-
      lapply(1:REPS, function (the.rep) {
        
        if (the.rep %% 250 == 0)
          cat(the.rep, " ")
        
        LSB.AR2.method.comp.one.rep(beta.tt1, beta.tt2, X, X.lsts, sigma.sq, N, S)
          
      })
  }
  
  exact.ests <- sapply(1:REPS, function(x) sims[[x]]$lsbar2.exact.est)
  whittle.ests <- sapply(1:REPS, function(x) sims[[x]]$lsbar2.whittle.est)
  lsts.ests <- sapply(1:REPS, function(x) sims[[x]]$lsts.ar2.est)
  
  list(exact.ests = exact.ests,
       whittle.ests = whittle.ests,
       lsts.ests = lsts.ests)
  
}
  
