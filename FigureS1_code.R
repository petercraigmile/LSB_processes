## ======================================================================
## Script: FigureS1_code.R
## Purpose: Plot the PACF for the 8 sections of the eeg data to get 
## an idea of the time varying nature of the ACF. 
## ======================================================================

eeg <- scan("EEG_data/eeg_F3.txt")

T. <- length(eeg)

ts <- seq(from=0, len=length(eeg), by=1/256)

# Separate the EEG series into 8 equal sections 
the.sec <- sapply(1:8, function(i) c((i-1)*(T./8)+1,(i*(T./8))))

pdf("figures/FigureS1_eeg_pacf.pdf", height = 4, width = 6.5)

par(
  mfrow = c(2, 4),
  cex = 0.75,
  mar = c(3.1, 2.5, 3.1, 0.5),
  mgp = c(1.3, 0.5, 0),
  bty = "n",
  oma = c(1, 0, 1, 0)
)

for (i in 1:8) {
  eeg.pacfs <- pacf(eeg[((i - 1) * (T. / 8) + 1):(i * (T. / 8))], 
                    plot = FALSE)
  plot(
    eeg.pacfs,
    ylim = c(-0.7, 1),
    main = paste0("t = ", the.sec[1, i], " - ", the.sec[2, i])
  )
}

dev.off()