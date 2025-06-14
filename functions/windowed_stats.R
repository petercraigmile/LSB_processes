
my.acf <- function (x, lag, type = "correlation") {

    acf(x, plot=FALSE, lag.max=lag, type = type)$acf[lag+1]
}


wstats <- function (x, FUN, B, ...) {

    N <- length(x)

    bs <- 1:B
    ks <- 0:(N-B-1)

    med.B <- median(bs)

    structure(list(indexes=ks+med.B,
                   stats=lapply(ks, function (k) FUN(x[k+bs], ...)),
                   x=x, N=N),
              class="wstats")
}


plot.wstats <- function (obj, ts, ...) {

    if (missing(ts)) {
        
        plot(obj$indexes, unlist(obj$stats), xlim=range(1:obj$N), type="l", ...)            
        
    } else {
        
        plot(ts[obj$indexes], unlist(obj$stats), xlim=range(ts), type="l", ...)
    }
}
