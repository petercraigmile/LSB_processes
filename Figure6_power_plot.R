
## ======================================================================
## Script: Figure6_power_plot.R
##
## Purpose: Produce Figure 6 comparing the power curves for our test
## of stationarity with the Paparoditis (2010) test for stationarity
## in the case of a LSB-AR(1) process, as the sample length T (N in
## this code) and slope beta_{11} is varied.
## ======================================================================


source("functions/load_save.R")
source("functions/Paparoditis_statistics.R")



load.it("Q300.0")
load.it("Q300.Z")
load.it("Q300")

load.it("Q500.0")
load.it("Q500")

load.it("Q1000.0")
load.it("Q1000")

load.it("Q2000.0")
load.it("Q2000")



beta1s <- seq(0, 1.5, 0.3)

pow300  <- sapply(Q300,  calc.Q.pow, Q0=Q300.0)
pow500  <- sapply(Q500,  calc.Q.pow, Q0=Q500.0)
pow1000 <- sapply(Q1000, calc.Q.pow, Q0=Q1000.0)
pow2000 <- sapply(Q2000, calc.Q.pow, Q0=Q2000.0)

pow <- rbind(pow300, pow500, pow1000, pow2000)



pdf(file="Figure6_power.pdf", width=6.3, height=2.5)
par(mfrow=c(1,2), cex=0.65, mar=c(3,3,1,0.5), mgp=c(1.8,0.5,0), bty="L")

x <- read.csv("power.csv")

Ns <- c(300, 500, 1000, 2000)

beta11s <- seq(0, 1.5, 0.3)

plot(beta11s, x[1,7:12], lty=4, lwd=2, ylim=c(0, 1),
     xlab=expression(beta[1][1]), ylab="Power", type="n")

abline(h=0.05, col="gray", lty=1, lwd=2)

for (j in 1:length(Ns)) {

    lines(beta11s, x[j, 7:12], type="l", lty=5-j, lwd=2)
}

mtext(side=3, "(a) Our test", line=0, cex=0.8)

plot(beta11s, pow[1,], lty=4, lwd=2, ylim=c(0, 1),
     xlab=expression(beta[1][1]), ylab="Power", type="n")

abline(h=0.05, col="gray", lty=1, lwd=2)

for (j in 1:length(Ns)) {

    lines(beta11s, pow[j,], type="l", lty=5-j, lwd=2)
}

mtext(side=3, "(b) Paparoditis test", line=0, cex=0.8)

legend(0.05, 1, paste("T = ", rev(Ns)), lty=1:4, lwd=2,  box.col="gray70")

dev.off()
