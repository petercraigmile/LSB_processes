

LSB.AR1.sim <- function (N, phi, sigma2=1) {

    if (length(phi)) {

        phi <- rep(phi, N)
    }

    if (N==1) {

        the.sds <- sqrt(sigma2 / (1.0 - phi[1]^2))

    } else {
        
        the.sds <- sqrt(sigma2 * c(1.0/(1.0-phi[1]^2), rep(1.0, N-1)))
    }

    z <- rnorm(N, 0, the.sds)

    x <- rep(NA, N)

    x[1] <- z[1]

    if (N > 1) {

        for (t in 2:N) {
            
            x[t] <- phi[t] * x[t-1] + z[t]
        }
    }
    
    return(x)
}



LSB.AR1.m2l <- function (beta, y, X.phi, X.sigma2) {
    
    K <- ncol(X.phi)
    
    phi.beta <- beta[1:K]
    phi      <- calc.phi(X.phi, phi.beta)

    sigma2.beta <- beta[-(1:K)]
    sigma2.eta  <- drop(X.sigma2 %*% sigma2.beta)
    sigma2      <- exp(sigma2.eta)

    N <- length(y)
    
    ts <- 2:N

    the.sds <- sqrt(sigma2 * c(1.0/(1.0-phi[1]^2), rep(1, N-1)))

    -2 * dnorm(y[1], 0, the.sds[1], log=TRUE) -
        2 * sum(dnorm(y[ts], phi[ts] * y[ts-1], the.sds[ts], log=TRUE))
}





LSB.AR1.ML <- function (y, X.phi, X.sigma2) {

    opt <- optim(rep(1e-3, ncol(X.phi)+ncol(X.sigma2)), LSB.AR1.m2l,
                 method="BFGS",
                 y=y, X.phi=X.phi, X.sigma2=X.sigma2)
    
    list(beta.hat=opt$par,
         m2l=opt$value)
}







LSB.AR1.sdf <- function (beta, fs, X.phi, X.sigma2, deltat=1) {

    K <- ncol(X.phi)
    
    phi.beta <- beta[1:K]
    phi      <- calc.phi(X.phi, phi.beta)

    sigma2.beta <- beta[-(1:K)]
    sigma2.eta  <- drop(X.sigma2 %*% sigma2.beta)
    sigma2      <- exp(sigma2.eta)

    t(sapply(1:length(phi), function (k)
        ar.sdf(fs, phi[k], sigma[k], deltat)))
}
