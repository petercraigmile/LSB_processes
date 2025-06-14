
## ======================================================================
## Script: Figure3_code.R
## Purpose: Plot choices of N for different T
## ======================================================================

pdf(file="figures/Figure3_NT_plot.pdf", width=5, height=2.5)

par(mfrow=c(1,1), cex=0.65, mar=c(3,3,1,0.5), mgp=c(1.7,0.5,0), bty="L")

n = seq(1, 100000000, 100000)
plot(
  n,
  n ^ .4,
  type = "l",
  xlab = "T",
  ylab = "N",
  lty = 2,
  xaxt = "n",
  col = "gray50"
)
lines(n, n ^ 0.35, lty = 3, col = "gray50")
lines(n, n ^ 0.3, lty = 4, col = "gray50")
lines(n, n ^ 0.25, col = "black", lwd = 2)
lines(n, n ^ 0.5 / log(n), col = "black", lwd = 2)

legend(
  1,
  1500,
  legend = c("N=T^0.4", "N=T^0.35", "N=T^0.3"),
  lty = 2:4,
  cex = 0.95,
  col = "gray50"
)

axis(
  1,
  c(0, 20000000, 40000000, 60000000, 80000000, 100000000),
  c(
    "0",
    paste0("2x", expression(10 ^ 7)),
    paste0("4x", expression(10 ^ 7)),
    paste0("6x", expression(10 ^ 7)),
    paste0("8x", expression(10 ^ 7)),
    paste0("1x", expression(10 ^ 8))
  ),
  cex = 0.65
)


dev.off()
