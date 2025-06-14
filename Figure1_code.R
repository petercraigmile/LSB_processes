
## ======================================================================
## Script: Figure1_code.R
## Purpose: Generate 3d plot of spectral density of LSB-AR and 
## LSB-FD processes
## ======================================================================


library(spectral)
library(fields)

source("functions/LSB_links.R")
source("functions/basis_functions.R")
source("functions/LSB_AR2.R")
source("functions/fd_processes.R")

## ======================================================================
## LSB AR 
## ======================================================================

# Specify length of series
N <- 500

# Genrate rescaled time points
us <- 0:(N-1)/N

# Select an orthogonal polynomial basis of order 2
ar.X <- poly.basis.for.sim(N)

# Basis parameters for LSB-AR(2)
ar.beta.1 <- c(0.61, 1.71, -1.27)
ar.beta.2 <- c(-3.52, 5.5, -3)

ar.sig.sq <- 1

theta.tt <- c(ar.beta.1, ar.beta.2, log(ar.sig.sq))

# Generate phis
AR2.phis <- LSB.AR2.PARtoAR(theta = theta.tt, X = ar.X)

phi.21 <- as.numeric(AR2.phis$phi.21)
phi.22 <- as.numeric(AR2.phis$phi.22)

# Specify frequencies
fs <- seq(0, 1/2, length=50)

# Time indices
ts <- 3:N

# Time stamps where sdf is calculated
tt.ar <- seq(ts[1], ts[N-2], length.out = 50)

# Calculate time varying sdfs at times tt.ar 
ar2.sdfs <- lapply(tt.ar,
                   function(k)
                     dB(ar.sdf(fs,
                               c(phi.21[k], phi.22[k]),
                               ar.sig.sq)))

# Organize into a matrix
sp.matrix.ar2 <- t(sapply(ar2.sdfs, function (x) x))

# Specify color palette 

palette.col = gray.colors(10, 0.8, 0.9)

jet.colors = colorRampPalette(palette.col)
nbcol = 100
color = jet.colors(nbcol)
z = sp.matrix.ar2
nrz = nrow(sp.matrix.ar2)
ncz = ncol(sp.matrix.ar2)
zfacet = z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol = cut(t(zfacet), nbcol)


## ======================================================================
## LSB FD 
## ======================================================================

# Select polynomial basis of order 1
fd.X <- poly.basis(N, orthogonal = FALSE, normalized = FALSE)

# Specify basis parameters
fd.beta <- c(0.1,2)

# Specify variance of the process
fd.sig.sq <- 1

# Calculate the time varying FD curve as basis function
the.ds <- calc.delta(fd.X, fd.beta)

# Frequencies
fs <- seq(0, 1/2, length=50)

# Time indices
tt.fd <- seq(1, N, length.out = 50)

# Calculate SDF for FD process at tt.fd
fd.sdfs <-
  lapply(tt.fd, function(k)
    dB(fd.sdf(fs[-1], the.ds[k], fd.sig.sq)))

# Convert to matrix for ease of plotting 
sp.matrix.fd <- t(sapply(fd.sdfs, function (x) x))


# Set up parameters for 3d plot of SDF for LSB FD example

z.fd = sp.matrix.fd
nrz.fd = nrow(sp.matrix.fd)
ncz.fd = ncol(sp.matrix.fd)
zfacet.fd = z.fd[-1, -1] + z.fd[-1, -ncz.fd] + 
  z.fd[-nrz.fd, -1] + z.fd[-nrz.fd, -ncz.fd]
facetcol.fd = cut(t(zfacet.fd), nbcol)

## ======================================================================
## Plot Figure 1
## ======================================================================

pdf("figures/Figure1_ar2fdsim.pdf", height = 5, width = 8)

par(mfrow=c(2,2), mar=c(3.1,3.1,1,0.5), mgp=c(1.8,0.5,0),
    bty="L", oma = c(0,0,1,0))

persp(
  2 * pi * fs,
  tt.ar / N,
  t(sp.matrix.ar2),
  theta = 45,
  phi = 10,
  xlab = "Frequency",
  ylab = "Rescaled time",
  zlab = "Spectral density",
  expand = 0.5,
  ticktype = "detailed",
  nticks = 6,
  main = "(a)",
  cex = 0.5,
  col = color[facetcol]
)

persp(
  2 * pi * fs[-1],
  tt.fd / N,
  t(sp.matrix.fd),
  theta = 45,
  phi = 10,
  xlab = "Frequency",
  ylab = "Rescaled time",
  zlab = "Spectral density",
  expand = 0.5,
  ticktype = "detailed",
  nticks = 6,
  col = color[facetcol.fd],
  main = "(b)",
  cex = 0.5
)

plot(
  x = us[-c(1, 2)],
  y = phi.21,
  type = "l",
  xlab = "Rescaled time",
  ylab = "",
  main = "(c)",
  bty = "L",
  cex = 0.5,
  ylim = c(-1, 1),
  lwd = 2
)
lines(x = us[-c(1, 2)],
      y = phi.22,
      col = "gray",
      lwd = 2)
abline(h = 0, lty = 2)

plot(
  x = us,
  y = the.ds,
  ylim = c(0, .5),
  xlab = "Rescaled time",
  ylab = ""
  #, xlab = "Rescaled time", ylab = "LRD curve"
  ,
  type = "l",
  main = "(d)",
  bty = "L",
  cex = 0.5,
  lwd = 2
)

dev.off()

