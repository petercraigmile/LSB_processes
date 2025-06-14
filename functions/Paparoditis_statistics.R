


u.range <- function (N, m, len) {

    L <- ceiling(m/2)/N

    seq(L, 1-L, length=len)    
}



Q.stat <- function (us, y, m, h,
                    lw.y = lw.pgram(y, b=h),
                    b=(m/length(y))^1/5 * h) {

    Q <- sapply(us, function (u) {
        
        lsp <- local.pgram(u, y, m)

        K <- length(lsp$freq)
        ks <- 2:(K-1)

        f.nznn.inds <- find.freqs(lsp$freq[ks], lw.y$freq)
        
        Y <- lsp$spec[ks] / lw.y$spec[f.nznn.inds] - 1
        
        lsp.sm <- ksmooth(lsp$freq[ks], Y, "normal", bandwidth=b)
        
        mean(lsp.sm$y^2)
    })

    out <- list(Q=Q, A=max(Q))

    return (out)
}
    


sim.Q.stats <- function (us, m, h, b, sim.fun, REPS=2000, CORES=6, every=100) {
    
    sims <- mclapply(1:REPS, function (k) {
        
        if (k%%every==0) cat(k, " ")
        
        y <- sim.fun()
        
        Q.stat(us, y, m, h, b=b)
        
    }, mc.cores=CORES)

    sims
}


calc.Q.crit <- function (Q0) {

    quantile(sapply(Q0$sims, function (x) x$A), 0.95)
}


calc.Q.pow <- function (Q, Q0) {

    crit <- calc.Q.crit(Q0)
    
    mean(sapply(Q$sims, function (x) x$A) >= crit)
}
