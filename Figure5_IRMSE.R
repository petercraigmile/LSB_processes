
## ======================================================================
## Script: Figure5_IRMSE.R
##
## Purpose: Produce Figure 5 comparing the integrated root mean square
## errors (IRMSEs) for parameter curves and the log spectral density
## function for the exact maximum likelihood estimates, our block
## Whittle estimator, and the Whittle estimator due to Palma and Olea
## (2010), when estimating an LSB-AR(2) process.
## ======================================================================



load("results/IRMSE_data.RData")

Ns <- 2^(7:12)

pchs <- 15:17

cols <- c("black", "gray40", "gray70")



one.panel <- function (exact, whittle, lsts,
                       ylim, title) {
    
    plot(log10(Ns), colMeans(exact),
         xlab=expression(log[10](T)), ylab="IRMSE", ylim=ylim,
         xlim=c(2,3.7), col=cols[1], pch=pchs[1])
    
    points(log10(Ns), colMeans(whittle), pch=pchs[2], col=cols[2])
    
    points(log10(Ns), colMeans(lsts), pch=pchs[3], col=cols[3])
    
    mtext(title, side=3, line=0, cex=0.7)

    max.SEs <- cbind(apply(exact, 2, sd) / sqrt(nrow(exact)),
                     apply(whittle, 2, sd) / sqrt(nrow(whittle)),
                     apply(lsts, 2, sd) / sqrt(nrow(lsts)))
    round(max(max.SEs), 3)
}


pdf(file="Figure5_IRMSE.pdf", width=6.3, height=2.5)
par(mfrow=c(1,3), cex=0.7, mar=c(3,3,1,0.5), mgp=c(1.8,0.5,0), bty="L")

one.panel(phi.21.irmse.exact,
          phi.21.irmse.whittle,          
          phi.21.irmse.lsts,
          c(0, 0.3),         
          "(a) AR(1) parameter curve")

one.panel(phi.22.irmse.exact,
          phi.22.irmse.whittle,          
          phi.22.irmse.lsts,
          c(0, 0.3),         
          "(b) AR(2) parameter curve")

one.panel(log.sdf.irmse.exact,
          log.sdf.irmse.whittle,          
          log.sdf.irmse.lsts,
          c(0.05, 0.6),
          "(c) Time varying log SDF")

legend(2.54, 0.6,
       rev(c("Maximum likelihood",
         "Block Whittle",
         "Palma and Olea")),
       col=rev(cols),
       pch=rev(pchs), cex=0.8, box.col="gray80")
       
dev.off()
