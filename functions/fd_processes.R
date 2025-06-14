## ======================================================================
## Copyright 1998--, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================

## ======================================================================
## File     : fd_processes.R
## Purpose  : R code for fractionally differenced (FD) processes
##            (also known as ARFIMA(0,d,0) processes).
## Version  : 1.0
## Requires : DaviesHarte.R
## Updated  : pfc@stat.osu.edu, May 2014.
##!
## References :
##
## 1. J. Beran, Statistics for Long Memory Processes. New York:
##    Chapman and Hall, 1994.
## ======================================================================


fd.acvs <- function (lag.max=NULL, d, sigma2=1, lags) {
  ## ======================================================================
  ## Purpose : Calculate the autocovariance sequence (ACVS) for an
  ##           FD('d','sigma2') process, from lags zero up to 'lag.max'.
  ##           if lag.max is NULL you can calculate the ACVS at lags 'lags'.
  ## Assumes : lag.max >= 0
  ## ======================================================================
  
  if (!is.null(lag.max)) {
    
    acvs0 <- sigma2 * exp(lgamma(1-2*d)-2*lgamma(1-d))
    if (lag.max>0) {
      ks <- 1:lag.max
      cumprod(c(acvs0, (ks-1+d)/(ks-d)))
    }
    else acvs0
  } else {

    sigma2 * exp(lgamma(1-2*d) + lgamma(lags+d) - lgamma(1-d) - lgamma(d) - lgamma(1+lags-d))
  }  
}





fd.acf <- fd.acs <- function (lag.max, d) {
  ## ======================================================================
  ## Purpose : Calculate the auto-correlation sequence for an
  ##           FD('d') process, from lags zero up to 'lag.max'. 
  ## Assumes : lag.max >= 0
  ## ======================================================================

  if (lag.max>0)
  {
    ks <- 1:lag.max
    cumprod(c(1, (ks-1+d)/(ks-d)))
  }
  else 1
}



fd.pacs <- fd.pacf <- function (lag.max, d) {
  ## ======================================================================
  ## Purpose : Calculate the partial auto-correlation sequence for an
  ##           FD('d') process, from lags one up to 'lag.max'.
  ##           See Hosking (1981)
  ## Assumes : lag.max >= 0
  ## ======================================================================

  d/(seq(lag.max)-d)
}



## fd.LD <- function (k, d) {
##   ## ======================================================================
##   ## Calculate the kth order Levinson-Durbin 
##   ## See page 168 of Hosking (1981).
##   ## ======================================================================
  
##   js <- 1:(k-1)

##     if (k==1) {
        
##         d/(1-d)
        
##     } else {
        
##     c(-exp(lgamma(k+1) + lgamma(js-d) + lgamma(k-d-js+1) -
##       lgamma(js+1) - lgamma(k-js+1) - lgamma(-d) - lgamma(k-d+1))),
##       d/(k-d))


## #      c(-gamma(k+1) gamma(js-d) * gamma(k-d-js+1) /
## #      (gamma(js+1) * gamma(k-js+1) * gamma(-d) * gamma(k-d+1)),
## #      d/(k-d))
## }
## }


fd.LD <- function (k, d) {
  ## ======================================================================
    ## Calculate the kth order Levinson-Durbin using log-gamma
    ## See page 168 of Hosking (1981).
    ## ======================================================================
    
    js <- 1:(k-1)
    
    if (k==1) {
        
        d/(1-d)
    } else {
        
        c(-exp(lgamma(k-d-js+1) + lgamma(k+1) + lgamma(js-d) - lgamma(js+1) - lgamma(k-js+1) - lgamma(k-d+1)) / gamma(-d),
          d/(k-d))
    }
}




fd.sdf <- function (fs, d, sigma2=1, deltat=1, deriv=NULL) {
  ## ======================================================================
  ## Purpose : calculates the spectrum density function of a FD(d, sigma2)
  ##           process at the frequencies 'fs', assuming a sampling rate
  ##           of 'deltat'.
  ##           - if 'deriv' is not equal to NULL calculate instead
  ##             the deriv'th derivative of this sdf with respect to d.
  ## ======================================================================
  
  tsy  <- abs(2 * sin(pi * fs * deltat))

  if (is.null(deriv)) {

    deltat * sigma2 * tsy^(-2*d)
  } else {
    
    deltat * sigma2 * tsy^(-2*d) * (-2*log(tsy))^deriv
  }
}



fd.fract.int.d <- function (d) {
  ## ======================================================================
  ## Purpose : calculates the integer and fractional differencing
  ##           components for a difference parameter d.
  ## ======================================================================
  
  int   <- floor(d+0.5)
  fract <- d-int
  list(d=d, int=int, fract=fract)
}





fd.sim.setup <- function (N, d, sigma2=1) {
  ## ===========================================================================
  ## Purpose : Setup for simulating a FD(d,sigma2) process of length N using the
  ##           Davies-Harte algorithm due to Chan and Wood (1994).
  ## Assumes : DaviesHarte.R is loaded.
  ## =========================================================================== 

  ds <- fd.fract.int.d(d)

  ## if fractional part of d close to close to zero
  if (abs(ds$fract)<1e-8)
    list(N=N, ds=ds, sigma=sqrt(sigma2), DH.obj=NULL)
  else
    list(ds=ds, sigma=sqrt(sigma2),
         DH.obj=Davies.Harte.sim.setup(N, fd.acvs, d=ds$fract, sigma2=sigma2))
}





fd.sim <- function (N, d, sigma2 = 1, setup = NULL) { ##, zs = NULL)
  ## ===========================================================================
  ## Purpose : Simulate a FD(d,sigma2) process of length N using the
  ##           Davies-Harte algorithm due to Chan and Wood (1994).
  ## Assumes : DaviesHarte.R is loaded.
  ## To do   : Add the capability to include innovations, 'zs'.           
  ## ===========================================================================
  
  if (is.null(setup$DH.obj))
    setup <- fd.sim.setup(N, d, sigma2)

#  if (is.null(zs))
    zs <- rnorm(setup$DH.obj$M)

  Davies.Harte.sim(DH.obj = setup$DH.obj, csum = setup$ds$int) ##, zs = zs)
}



