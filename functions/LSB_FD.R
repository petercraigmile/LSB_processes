
LSB.FD.sim.setup <- function (N, delta) {
    ## ======================================================================
    ## Setup code to simulate LSB-FB processes of length 'N' with
    ## time-varying fractional parameter 'delta'
    ## ======================================================================
    
    ## 'Npower2' is the next power of two greater than or equal to 'N'.
    Npower2 <- 2^ceiling(log2(N))

    K <- length(delta)
    M <- 2 * Npower2

    if (K < Npower2) {

        delta.ext <- c(delta, rep(delta[K], Npower2-K))
    } else {
        
        delta.ext <- delta
    }

    acvs.list <- lapply(1:Npower2, function (k)
        fd.acvs(Npower2, delta.ext[k], sigma2=1))

    Sk <- lapply(acvs.list, function (acvs)
        Re(fft(c(acvs, acvs[Npower2:2]), inverse=FALSE)))
    
    if (any(unlist(Sk)<=0)) stop("Some values in Davies.Harte.sim.setup are <= 0!!")  
    
    sqrt.Sk <- lapply(Sk, function (xx) sqrt(xx))

    list(N=N, Npower2=Npower2, M=M, sqrt.Sk=sqrt.Sk)
}



LSB.FD.sim <- function (N, delta, sigma2=1, setup=LSB.FD.sim.setup(N, delta)) {
    ## ======================================================================
    ## Simulate LSB-FB processes of length 'N' with
    ## time-varying fractional parameter 'delta' and variance 'sigma2'
    ##
    ## The argument 'setup' can be used for repeated simulations.
    ## ======================================================================
  
    ## Simulate the common set of random normals
    zs <- rnorm(setup$M)

    ks <- 2:setup$Npower2

    A <- sqrt(0.5) * complex(real=zs[2*ks-2], imag=zs[2*ks-1])
    
    Y <- lapply(setup$sqrt.Sk, function (ssk) c(ssk[1] * zs[1],
                                                ssk[ks] * A,
                                                ssk[setup$Npower2+1] * zs[setup$M],
                                                Conj(rev(ssk[ks] * A))))

    ## Now use the cut and paste scheme
    sqrt(sigma2) * sapply(1:setup$N, function (j)
        Re(fft(Y[[j]], inverse=TRUE))[j])/sqrt(setup$M)
}



LSB.FD.m2l <- function (beta, y, X.delta, X.sigma2) {

    K <- ncol(X.delta)
    
    delta.beta <- beta[1:K]
    delta.eta  <- drop(X.delta %*% delta.beta)
    delta.exp.eta <- exp(delta.eta)
    delta <- (delta.exp.eta - 1) / (2 * (1 + delta.exp.eta))
    ##delta      <- eta.to.delta(delta.eta)

    sigma2.beta <- beta[-(1:K)]
    sigma2.eta  <- drop(X.sigma2 %*% sigma2.beta)
    sigma2      <- exp(sigma2.eta)
    ##sigma2      <- eta.to.sigma2(sigma2.eta)

    gamma0 <- fd.acvs(0, delta, sigma2)

    N <- length(y)

    ll <- dnorm(y[1], 0, sqrt(gamma0[1]), log=TRUE)

    if (N > 1) {

        for (k in 2:N) {

            LD.k   <- fd.LD(k-1, delta[k])
            pacf.k <- fd.pacf(k-1, delta[k])            
            v.k    <- gamma0[k] * prod(1.0 - pacf.k^2)

            cond.mean.k <- sum(y[(k-1):1] * LD.k)

            ll <- ll + dnorm(y[k], cond.mean.k, sqrt(v.k), log=TRUE)
        }
    }

    -2 * ll
}





LSB.FD.ML <- function (y, X.delta, X.sigma2) {

    opt <- optim(rep(1e-3, ncol(X.delta)+ncol(X.sigma2)), LSB.FD.m2l,
                 method="L-BFGS-B",
                 y=y, X.delta=X.delta, X.sigma2=X.sigma2)

    list(beta.hat=opt$par,
         m2l=opt$value)
}




LSB.FD.sdf <- function (beta, fs, X.delta, X.sigma2, deltat=1) {

    p1 <- ncol(X.delta)
    p2 <- ncol(X.sigma2)

    delta.beta <- beta[1:p1]  
    delta.eta  <- drop(X.delta %*% delta.beta)
    delta.exp.eta <- exp(delta.eta)
    delta <- (delta.exp.eta - 1) / (2 * (1 + delta.exp.eta))

    sigma2.beta <- beta[(p1+1):(p1+p2)]
    sigma2.eta  <- drop(X.sigma2 %*% sigma2.beta)
    sigma2      <- exp(sigma2.eta)

    t(sapply(1:length(delta), function (k)
        fd.sdf(fs, delta[k], sigma[k], deltat)))
}

