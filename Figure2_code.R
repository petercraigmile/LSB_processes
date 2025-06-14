
## ======================================================================
## Script: Figure2_code.R
## Purpose: Calculate RMSEs for N, S combinations of Whittle likelihood 
## method to select ideal N and S for given length of time series
## ======================================================================

source("functions/LSB_AR2_sim_functions.R")
source("functions/basis_functions.R")
source("functions/sim_stats.R")

library(fields)

# Basis parameters for LSB-AR(2)
beta_tt1 <- c(0.61, 1.71, -1.27)
beta_tt2 <- c(-3.52, 5.5, -3)

sigma.sq = 1

theta.tt <- c(beta_tt1, beta_tt2, log(sigma.sq))

# Number of sims 
REPS <- 100

# Seq of N (logN/logT) for which we calculate RMSE
logN.logT.seq <- seq(from = .3, to = .8, by = .05)

# Seq of S (S/N) for which we calculate RMSE
S.N.seq <- seq(from = 0.3, to = 0.79, by = 0.07)

# Calculate expanded grid of N and S
the.NS.list.512 <- NS.list.calc(
  n = 512,
  logN.logT.seq = seq(from = .3, to = .8, by = .05),
  S.N.seq = seq(from = 0.3, to = 0.79, by = 0.07)
)

the.NS.list.1024 <- NS.list.calc(
  n = 1024,
  logN.logT.seq = seq(from = .3, to = .8, by = .05),
  S.N.seq = seq(from = 0.3, to = 0.79, by = 0.07)
)

# We use orthogonal polynomial basis 
X.512 <- poly.basis.for.sim(N = 512)
X.1024 <- poly.basis.for.sim(N = 1024)

# Calculate RMSEs for N, S combos
NS.512 <- mapply(function(i, j)
  NS.rmse.whittle(n = 512, i, j, X.512, REPS),
  the.NS.list.512$Ns,
  the.NS.list.512$Ss)

NS.1024 <- mapply(function(i, j)
  NS.rmse.whittle(n = 1024, i, j, X.1024, REPS),
  the.NS.list.1024$Ns,
  the.NS.list.1024$Ss)

# Extract RMSEs
NS.512.mse <- matrix(
  NS.512[3, ],
  nrow = length(logN.logT.seq),
  ncol = length(S.N.seq),
  byrow = TRUE
)

NS.1024.mse <- matrix(
  NS.1024[3, ],
  nrow = length(logN.logT.seq),
  ncol = length(S.N.seq),
  byrow = TRUE
)

save(NS.512.mse, NS.1024.mse,
     file = "results/NS_selection_data.RData")

## ======================================================================
## Plot Figure 2
## ======================================================================

zlim.min <- round(min(min(NS.512.mse), min(NS.1024.mse)) - 0.1, 1)
zlim.max <- round(max(max(NS.512.mse), max(NS.1024.mse)) + 0.1, 1)

pdf("figures/Figure2_NS_selection.pdf", height = 3.5, width = 6.5)

par(mfrow=c(1,2), cex=0.75, mar=c(3,3,3,2), mgp=c(1.8,0.5,0), bty="L")

fields::image.plot(
  x = logN.logT.seq,
  y = S.N.seq,
  z = NS.512.mse,
  xlab = "log(N)/log(T)",
  ylab = "S/N",
  zlim = c(zlim.min, zlim.max),
  main = "T = 512",
  col = gray(seq(0.1, 0.9, len = 200))
)

par(cex = 0.75)

fields::image.plot(
  x = logN.logT.seq,
  y = S.N.seq,
  z = NS.1024.mse,
  xlab = "log(N)/log(T)",
  ylab = "S/N",
  zlim = c(zlim.min, zlim.max),
  main = "T = 1024",
  col = gray(seq(0.1, 0.9, len = 200))
)

dev.off()

