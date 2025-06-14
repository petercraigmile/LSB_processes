
## ======================================================================
## Script: Figure4_code.R
## Purpose: Plot asymptotic results for Exact and
## Whittle likelihood results
## ======================================================================

## ======================================================================
## Plot Figure 4
## ======================================================================

load("results/asysd.RData")
load("results/boot_cis.RData")

true.var <- Re(asysd)

graph.names <- expression(beta[10], beta[11], beta[12], beta[20], beta[21], 
                          beta[22], beta[30])

log10.ns <- log10(ns)

pdf("figures/Figure4_tvar2_asymp.pdf", height = 11, width = 8)

par(mfcol=c(7,3), cex=0.7, mar=c(2,1.8,0.3,0.5), 
    mgp=c(1.8,0.5,0), bty="L", oma=c(1.2,1.4,1.2,0))

for(i in 1:length(true.param)){
  
  plot(
    x = log10.ns,
    y = boot.ci.exact$bias[i, ],
    type = "n",
    ylim = c(
      min(boot.ci.exact$bias.lcb[i, ], boot.ci.whittle$bias.lcb[i, ]),
      max(boot.ci.exact$bias.ucb[i, ], boot.ci.whittle$bias.ucb[i, ])
    ),
    xlab = "",
    ylab = "",
    pch = 3,
    bty = "l",
    cex = 1.2,
    cex.main = 1.6,
    col = "black"
  )
  
  segments(log10.ns-0.03, boot.ci.whittle$bias[i,], 
           log10.ns+0.03, boot.ci.whittle$bias[i,], 
           lwd=3, col="gray")
  
  segments(log10.ns, boot.ci.whittle$bias.lcb[i,], 
           log10.ns, boot.ci.whittle$bias.ucb[i,], 
           lwd=3, col="gray")
  
  segments(log10.ns - 0.03, boot.ci.exact$bias[i,], 
           log10.ns+0.03, boot.ci.exact$bias[i,], 
           lwd=3)
  
  segments(log10.ns, boot.ci.exact$bias.lcb[i,],
           log10.ns, boot.ci.exact$bias.ucb[i,], 
           lwd=3)
  
  
  abline(h = 0, lty = 2)
  
}

for(i in 1:length(true.param)){
  
  plot(
    x = log10.ns,
    y = boot.ci.exact$RMSE[i, ],
    type = "n",
    ylim = c(min(c(boot.ci.whittle$RMSE.lcb[i, ]), boot.ci.exact$RMSE.lcb[i, ]),
             max(c(boot.ci.whittle$RMSE.ucb[i, ]), boot.ci.exact$RMSE.ucb[i, ])),
    xlab = "",
    ylab = "",
    pch = 3,
    bty = "l",
    cex = 1.2,
    cex.main = 1.6,
    col = "black"
  )
  
  segments(log10.ns-0.03, boot.ci.whittle$RMSE[i,], 
           log10.ns+0.03, boot.ci.whittle$RMSE[i,], 
           lwd=3, col="gray")
  
  segments(log10.ns, boot.ci.whittle$RMSE.lcb[i,], 
           log10.ns, boot.ci.whittle$RMSE.ucb[i,], 
           lwd=3, col="gray")
  
  segments(log10.ns-0.03, boot.ci.exact$RMSE[i,], 
           log10.ns+0.03, boot.ci.exact$RMSE[i,], 
           lwd=3)
  
  segments(log10.ns, boot.ci.exact$RMSE.lcb[i,], 
           log10.ns, boot.ci.exact$RMSE.ucb[i,], 
           lwd=3)
  
}



for(i in 1:length(true.param)){
  
  plot(
    x = log10.ns,
    y = boot.ci.exact$sqrtn.RMSE[i, ],
    type = "n",
    ylim = c(min(c(boot.ci.whittle$sqrtn.RMSE.lcb[i, ]),
                 boot.ci.exact$sqrtn.RMSE.lcb[i, ]),
             max(c(boot.ci.whittle$sqrtn.RMSE.ucb[i, ]), 
                 boot.ci.exact$sqrtn.RMSE.ucb[i, ])),
    xlab = "",
    ylab = "",
    pch = 3,
    bty = "l",
    cex = 1.2,
    cex.main = 1.6,
    col = "black"
  )
  
  segments(log10.ns-0.03, boot.ci.whittle$sqrtn.RMSE[i,], 
           log10.ns+0.03, boot.ci.whittle$sqrtn.RMSE[i,], 
           lwd=3, col="gray")
  
  segments(log10.ns, boot.ci.whittle$sqrtn.RMSE.lcb[i,], 
           log10.ns, boot.ci.whittle$sqrtn.RMSE.ucb[i,], 
           lwd=3, col="gray")
  
  segments(log10.ns-0.03, boot.ci.exact$sqrtn.RMSE[i,],
           log10.ns+0.03, boot.ci.exact$sqrtn.RMSE[i,], 
           lwd=3)
  
  segments(log10.ns, boot.ci.exact$sqrtn.RMSE.lcb[i,], 
           log10.ns, boot.ci.exact$sqrtn.RMSE.ucb[i,], 
           lwd=3)
  
  abline(h = true.var[i], lty = 2)
  
}

mtext("Bias", side=3, at=0.19, outer=TRUE, cex=0.85)
mtext("RMSE", side=3, at=0.53, outer=TRUE, cex=0.85)
mtext(expression(sqrt(T) %*% RMSE), side=3, at=0.85, outer=TRUE, cex=0.85)

mtext(expression(log[10](T)), side=1, at=0.19, outer=TRUE, cex=0.85)
mtext(expression(log[10](T)), side=1, at=0.53, outer=TRUE, cex=0.85)
mtext(expression(log[10](T)), side=1, at=0.85, outer=TRUE, cex=0.85)

for(i in 1:length(true.param)){
  
  mtext(graph.names[i], side=2, at=1.09-0.144*i, outer=TRUE, cex=0.75, line=0)
  
}

dev.off()
