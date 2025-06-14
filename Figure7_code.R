
## ======================================================================
## Script: Figure7_code.R
## Purpose: Plot eeg series
## ======================================================================

eeg <- scan("EEG_data/eeg_F3.txt")

T. <- length(eeg)

ts <- seq(from=0, len=length(eeg), by=1/256)


pdf("figures/Figure7_eeg_naive.pdf", height = 1.5, width = 6.5)

par(mfrow=c(1,1), cex=0.75, mar=c(3.1,3.1,1,0.5), mgp=c(1.8,0.5,0), bty="L")

plot(ts, eeg, type="l", xlab = "Time (secs)", 
     ylab = "EEG (millivolts)")

dev.off()
