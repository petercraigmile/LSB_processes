
## ======================================================================
## Script: application2_fit_model.R
## Purpose: Fit chosen LSB-AR(p) model to EEG data and calculate time
## varying SDF. Additionally, compute the test statistic for the test of
## stationarity in Section 4 for EEG series
## ======================================================================

library(spectral)
library(fields)
library(mvnfast)


source("functions/basis_functions.R")
source("functions/LSB_ARp.R")
source("functions/LSB_AR_SDF.R")

load("results/models_orders.RData")

eeg <- scan("EEG_data/eeg_F3.txt")

T. <- length(eeg)

ts <- seq(from=0, len=length(eeg), by=1/256)

fs <- seq(0, 1/2, length=256)

## Fit selected LSBAR(p) model based on NIC

p = nic.model$ar.order

b = nic.model$basis.order

s = nic.var.basis.order

beta_ord = c(rep(b, p), s)

# Build orthonormal spline basis 
X <- lapply(1:length(beta_ord), function(i)
  build.X(
    n = T.,
    X.deg = beta_ord[i],
    X.method = "spline"
  ))

# Model fit
fit <- optim(
  rep(0.01, sum(beta_ord + 1)),
  LSB.ARp.lik,
  y = eeg,
  X = X,
  method = "L-BFGS-B",
  hessian = T
)

# Estimated cov matrix
cov.mat <- solve(fit$hessian)

## Calculate the AR curve estimates and time varying variance curve

the.beta <- theta.to.beta.list(X, fit$par, p)

the.phi <- beta.list.to.phi.list(X, the.beta, p)

the.sigma.sq <- eta.to.sigma2(calc.eta(X[[p + 1]], the.beta[[p + 1]]))

## Calculate time varying SDF

ar.sdfs <-
  lapply(1:(length(eeg) - p), function(k)
    dB(ar.sdf(
      freqs = fs,
      phi = as.vector(the.phi[[p]][1:p, k]),
      sigma2 = the.sigma.sq[k]
    )))


sp.matrix <- t(sapply(ar.sdfs, function (x) x))

pdf(file = "figures/eeg_8_2_2.pdf", height = 3.5, width = 6.5)

par(mfrow=c(1,1), cex=0.75, mar=c(3.1,3.1,1,0.5), 
    mgp=c(1.8,0.5,0), bty="L", oma = c(1,0,1,0))

image.plot(
  ts[2:(length(eeg) - nic.model$ar.order + 1)],
  256 * fs,
  sp.matrix,
  zlim = c(10, 70),
  xlab = "time",
  ylab = "Frequency",
  ylim = c(0, 31)
)
abline(h = c(3.5, 7.5, 14.5, 31), lwd = 2)
text(12, 2, expression(delta))
text(12, 6, expression(theta))
text(12, 12, expression(alpha))
text(12, 24, expression(beta))

dev.off()



## Calculating the confidence intervals

par.sim <- mvnfast::rmvn(
  1000,
  mu = fit$par,
  sigma = cov.mat,
  ncores = parallel::detectCores() - 2
)

sdfs <- mclapply(1:nrow(par.sim),
                 function(i)
                   LSB.ARp.sdf(par = par.sim[i,], X = X, fs = fs),
                 mc.cores = parallel::detectCores() - 1)

# Calculate sd of time varying SDF
spec.sd <-
  array(sapply(1:(nrow(sdfs[[1]]) * ncol(sdfs[[1]])),
               function(x) {
                 sd(unlist(lapply(sdfs, '[[', x)))
               }), dim(sdfs[[1]]))


image.plot(
  ts[2:(length(eeg) - p + 1)],
  256 * fs,
  spec.sd,
  xlab = "Time(seconds)",
  ylab = "Frequency",
  ylim = c(0, 31),
  cex = 0.75,
  zlim = c(0, 3.5)
)
abline(h = c(3.5, 7.5, 14.5, 31), lwd = 2)
text(12, 2, expression(delta))
text(12, 6, expression(theta))
text(12, 12, expression(alpha))
text(12, 24, expression(beta))

save(p, sp.matrix, spec.sd, 
     file = "results/sdf_data_fig8.RData")

# Fit stationary model 

p = nic.model$ar.order

beta_ord.stat = c(rep(0, p), 0)

X.stat <- lapply(1:length(beta_ord.stat), function(i)
  build.X(n = T., X.method = "spline", X.deg = beta_ord.stat[i]))


fit.stat <- optim(
  rep(0.01, sum(beta_ord.stat + 1)),
  LSB.ARp.lik,
  y = eeg,
  X = X.stat,
  method = "L-BFGS-B",
  hessian = T
)

# Test statistic for testing stationarity
test.stat <- -2*(fit$value - fit.stat$value)
