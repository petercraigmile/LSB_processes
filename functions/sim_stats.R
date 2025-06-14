source("functions/check_causal_AR2.R")


calc.stats <- function(n, samp, theta, digits = 4) {
  ## ======================================================================
  ## Calculate statistics from bootstrapped samples 
  ## ======================================================================
  
  the.par <- round(rowMeans(samp), digits)
  
  the.bias <-
    round(rowMeans(samp) - theta, digits)
  
  the.sd <- round(apply(samp, 1, sd), digits)
  
  the.rmse <-
    round(sqrt(rowMeans((samp - theta) ^ 2)),
          digits)
  
  cbind(
    par.est = the.par,
    bias = the.bias,
    SD = the.sd,
    RMSE = the.rmse,
    sqrtn.RMSE = round(sqrt(n) * the.rmse, digits)
  )
}

calc.sim.irmses <- function(true.param, sims, X,
                            X.lsts = NULL, fs, REPS, lsts = FALSE){
  ## ======================================================================
  ## Calculate IRMSEs for the AR2 parameter curves phi21 and phi22 and 
  ## log(SDF) based on LSB AR2 sims and true.param
  ## ======================================================================
  
  require(parallel)
  require(pbmcapply)
  
  # Calculate true phis and log.sdf 
  true.phis <- LSB.AR2.PARtoAR(true.param, X)
  
  true.sigma.sq <- eta.to.sigma2(true.param[length(true.param)])

  true.log.sdf <- LSB.AR2.log.sdf(true.phis$phi.21, 
                                  true.phis$phi.22, 
                                  true.sigma.sq, 
                                  fs)
  
  # Calculate IRMSE for phis
  if(lsts == FALSE){
    calc.ar.phis <- pbmclapply(1:REPS, function (k) {
      LSB.AR2.PARtoAR(sims[, k], X)
    }, mc.cores = parallel::detectCores() - 2)
  } else {
    calc.ar.phis <- pbmclapply(1:REPS, function (k) {
      LSTS.AR2.beta.to.phi(sims[, k], X.lsts)
    }, mc.cores = parallel::detectCores() - 2)
  }

  
  phi.21.ests <- sapply(1:REPS, function(x) calc.ar.phis[[x]]$phi.21)
  phi.22.ests <- sapply(1:REPS, function(x) calc.ar.phis[[x]]$phi.22)
  
  IRMSE.phi.21 <- sqrt(colMeans((phi.21.ests - as.numeric(true.phis$phi.21))^2))
  IRMSE.phi.22 <- sqrt(colMeans((phi.22.ests - as.numeric(true.phis$phi.22))^2))
  
  if(lsts == FALSE) {
    calc.log.sdf <- pbmclapply(1:REPS, function (k) {
      
      LSB.AR2.log.sdf(phi.21.ests[,k], phi.22.ests[,k], 
                      eta.to.sigma2(sims[nrow(sims), k]),
                      fs)
      
    }, mc.cores = parallel::detectCores() - 2)
  } else {
    calc.log.sdf <- pbmclapply(1:REPS, function (k) {
      
      LSB.AR2.log.sdf(phi.21.ests[,k], phi.22.ests[,k], 
                      sims[nrow(sims), k], fs)
      
    }, mc.cores = parallel::detectCores() - 2)
  }
  
  IRMSE.log.sdf <-
    sapply(1:REPS, 
           function(x) sqrt(mean((calc.log.sdf[[x]] - true.log.sdf) ^ 2)))
  
  list(
    ar1.ests = phi.21.ests,
    ar2.ests = phi.22.ests,
    log.sdf.ests = calc.log.sdf,
    IRMSEs = cbind(
      IRMSE.phi.21 = IRMSE.phi.21,
      IRMSE.phi.22 = IRMSE.phi.22,
      IRMSE.log.sdf = IRMSE.log.sdf
    )
  )

  
}

boot.ci <- function(asymp.stats, boot.stats) {
  ## ======================================================================
  ## Calculate bias, rmse, sqrtn.rmse and 95% CIs from bootstrapped samples 
  ## ======================================================================
  
  the.bias <- sapply(1:length(asymp.stats), 
                       function(x) asymp.stats[[x]][, "bias"])
  the.rmse <- sapply(1:length(asymp.stats), 
                       function(x) asymp.stats[[x]][, "RMSE"])
  the.sqrtn.rmse <- sapply(1:length(asymp.stats), 
                       function(x) asymp.stats[[x]][, "sqrtn.RMSE"])
  
  the.bias.ucb <- sapply(1:length(boot.stats), function(y) {
    apply(sapply(1:length(boot.stats[[y]]), function(x) 
      boot.stats[[y]][[x]][, "bias"]), 1, quantile, 0.975)
  })
  
  the.bias.lcb <- sapply(1:length(boot.stats), function(y) {
    apply(sapply(1:length(boot.stats[[y]]), function(x) 
      boot.stats[[y]][[x]][, "bias"]), 1, quantile, 0.025)
  })
  
  the.rmse.ucb <- sapply(1:length(boot.stats), function(y) {
    apply(sapply(1:length(boot.stats[[y]]), function(x) 
      boot.stats[[y]][[x]][, "RMSE"]), 1, quantile, 0.975)
  })
  
  the.rmse.lcb <- sapply(1:length(boot.stats), function(y) {
    apply(sapply(1:length(boot.stats[[y]]), function(x) 
      boot.stats[[y]][[x]][, "RMSE"]), 1, quantile, 0.025)
  })
  
  the.sqrtn.rmse.ucb <- sapply(1:length(boot.stats), function(y) {
    apply(sapply(1:length(boot.stats[[y]]), function(x) 
      boot.stats[[y]][[x]][, "sqrtn.RMSE"]), 1, quantile, 0.975)
  })
  
  the.sqrtn.rmse.lcb <- sapply(1:length(boot.stats), function(y) {
    apply(sapply(1:length(boot.stats[[y]]), function(x) 
      boot.stats[[y]][[x]][, "sqrtn.RMSE"]), 1, quantile, 0.025)
  })
  
  list(
    bias = the.bias,
    bias.lcb = the.bias.lcb,
    bias.ucb = the.bias.ucb,
    RMSE = the.rmse,
    RMSE.lcb = the.rmse.lcb,
    RMSE.ucb = the.rmse.ucb,
    sqrtn.RMSE = the.sqrtn.rmse,
    sqrtn.RMSE.lcb = the.sqrtn.rmse.lcb,
    sqrtn.RMSE.ucb = the.sqrtn.rmse.ucb
  )
  
}

calc.causal.percent <- function(method, n.ind) {
  ## ======================================================================
  ## Calculate pecentage of AR curve estimates which are causal
  ## of PFC check.causal function
  ## ======================================================================
  
  if(method == "exact"){
    mean(sapply(1:REPS,
                function(the.rep) {
                  mean(check.causal.AR2(
                    irmse.stats.exact[[n.ind]]$ar1.ests[, the.rep],
                    irmse.stats.exact[[n.ind]]$ar2.ests[, the.rep]
                  ))
                }))
  } else if (method == "whittle") {
    mean(sapply(1:REPS,
                function(the.rep) {
                  mean(check.causal.AR2(
                    irmse.stats.whittle[[n.ind]]$ar1.ests[, the.rep],
                    irmse.stats.whittle[[n.ind]]$ar2.ests[, the.rep]
                  ))
                }))
  } else if (method == "lsts"){
    mean(sapply(1:REPS,
                function(the.rep) {
                  mean(check.causal.AR2(
                    irmse.stats.lsts[[n.ind]]$ar1.ests[, the.rep],
                    irmse.stats.lsts[[n.ind]]$ar2.ests[, the.rep]
                  ))
                }))
  }
  
}

NS.rmse.whittle <- function(n, N, S, X, REPS) {
  
  ## ======================================================================
  ## Calculate mean RMSE corresponding to whittle estimate sims for given
  ## N and S. 
  ## ======================================================================
  
  NS.stat <- LSB.AR2.est.combine.reps(
    n = n,
    beta.tt1 = beta_tt1,
    beta.tt2 = beta_tt2,
    sigma.sq = sigma.sq,
    X = X,
    lik.method = "whittle",
    N = N,
    S = S,
    REPS = REPS
  )
  
  cbind(N, S,
        mean.RMSE = round(mean(calc.stats(
          n, NS.stat, theta.tt
        )[, "RMSE"]), 4))
  
}

NS.list.calc <- function(n, logN.logT.seq, S.N.seq){
  
  ## ======================================================================
  ## Calculate expanded grid on N and S combinations given length of 
  ## time series n and sequences for log(N)/log(T) and S/N
  ## ======================================================================
  
  Ns = trunc(n^logN.logT.seq)
  
  Ss = data.frame(sapply(Ns, function(N) trunc(S.N.seq*N)))
  
  names(Ss) <- Ns
  
  NS.list <- tidyr::pivot_longer(Ss, everything(),
                                 names_to = "Ns", values_to = "Ss")
  
  NS.list$Ns <- as.numeric(NS.list$Ns)
  
  NS.list = NS.list[order(NS.list$Ns),]
  
}

