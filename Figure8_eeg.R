
## ======================================================================
## Script: Figure8_code.R
## Purpose: Generate naive windowed estimate of time varying SDF, 
## LSB-AR(p) estimates from model chosen using NIC and its corresponding 
## time varying SDs for the eeg series 
## ======================================================================

library(spectral)
library(fields)


source("functions/basis_functions.R")
source("functions/windowed_stats.R")
source("functions/LSB_ARp.R")
source("functions/LSB_AR_SDF.R")

load("results/sdf_data_fig8.RData")

eeg <- scan("EEG_data/eeg_F3.txt")

T. <- length(eeg)

ts <- seq(from=0, len=length(eeg), by=1/256)

fs <- seq(0, 1/2, length=256)

# Calculate naive windowed SDF estimates of eeg series
sps.naive <- wstats(eeg, function (x) pgram(x, delta.t=1/256), B=512)

sp.matrix.naive <- t(dB(sapply(sps.naive$stats, function (x) x$spec)))

## ======================================================================
## Plot Figure 8
## ======================================================================

pdf("figures/Figure8_windowed_ar_gray.pdf", height = 3.5, width = 9.5)

par(mfrow=c(1,3), cex=0.75, mar=c(3.1,3.1,2.1,3.1),
    mgp=c(1.8,0.5,0), bty="L", oma = c(0.5,1,0.5,1))

image.plot(
  ts[sps.naive$indexes],
  sps.naive$stats[[1]]$freq,
  sp.matrix.naive,
  ylim = c(0, 31),
  zlim = c(10, 70),
  xlab = "Time(secs)",
  ylab = "Frequency(dB)",
  cex = 0.75,
  main = "(a)",
  col = gray.colors(200, start = 0.01, end = 0.99)
)
abline(h = c(3.5, 7.5, 14.5, 31), lwd = 2)
text(12, 2, expression(delta))
text(12, 6, expression(theta))
text(12, 12, expression(alpha))
text(12, 24, expression(beta))

image.plot(
  ts[2:(length(eeg) - p + 1)],
  256 * fs,
  sp.matrix,
  ylim = c(0, 31),
  zlim = c(10, 70),
  xlab = "Time(secs)",
  ylab = "Frequency(dB)",
  cex = 0.75,
  main = "(b)",
  col = gray.colors(200, start = 0.01, end = 0.99)
)
abline(h = c(3.5, 7.5, 14.5, 31), lwd = 2)
text(12, 2, expression(delta))
text(12, 6, expression(theta))
text(12, 12, expression(alpha))
text(12, 24, expression(beta))

image.plot(
  ts[2:(length(eeg) - p + 1)],
  256 * fs,
  spec.sd,
  xlab = "Time(secs)",
  ylab = "Frequency(dB)",
  ylim = c(0, 31),
  zlim = c(0, 3.5),
  cex = 0.75,
  main = "(c)",
  col = gray.colors(200, start = 0.01, end = 0.99)
)
abline(h = c(3.5, 7.5, 14.5, 31), lwd = 2)
text(12, 2, expression(delta))
text(12, 6, expression(theta))
text(12, 12, expression(alpha))
text(12, 24, expression(beta))

dev.off()


## ======================================================================
## Plot Figure 8 for supplement
## ======================================================================


pdf("figures/Figure8_windowed_ar_color.pdf", height = 3.5, width = 9.5)

par(mfrow=c(1,3), cex=0.75, mar=c(3.1,3.1,2.1,3.1),
    mgp=c(1.8,0.5,0), bty="L", oma = c(0.5,1,0.5,1))

image.plot(
  ts[sps.naive$indexes],
  sps.naive$stats[[1]]$freq,
  sp.matrix.naive,
  ylim = c(0, 31),
  zlim = c(10, 70),
  xlab = "Time(secs)",
  ylab = "Frequency(dB)",
  cex = 0.75,
  main = "(a)"
)
abline(h = c(3.5, 7.5, 14.5, 31), lwd = 2)
text(12, 2, expression(delta))
text(12, 6, expression(theta))
text(12, 12, expression(alpha))
text(12, 24, expression(beta))

image.plot(
  ts[2:(length(eeg) - p + 1)],
  256 * fs,
  sp.matrix,
  zlim = c(10, 70),
  xlab = "Time(secs)",
  ylab = "Frequency(dB)",
  ylim = c(0, 31),
  cex = 0.75,
  main = "(b)"
)
abline(h = c(3.5, 7.5, 14.5, 31), lwd = 2)
text(12, 2, expression(delta))
text(12, 6, expression(theta))
text(12, 12, expression(alpha))
text(12, 24, expression(beta))

image.plot(
  ts[2:(length(eeg) - p + 1)],
  256 * fs,
  spec.sd,
  xlab = "Time(secs)",
  ylab = "Frequency(dB)",
  ylim = c(0, 31),
  cex = 0.75,
  zlim = c(0, 3.5),
  main = "(c)"
)
abline(h = c(3.5, 7.5, 14.5, 31), lwd = 2)
text(12, 2, expression(delta))
text(12, 6, expression(theta))
text(12, 12, expression(alpha))
text(12, 24, expression(beta))

dev.off()
