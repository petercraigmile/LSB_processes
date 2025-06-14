
check.causal.AR2 <- function (phi1, phi2) {

    m <- length(phi1)

    sapply(1:m, function (k)
        all(Mod(polyroot(c(1, -phi1[k], -phi2[k]))) > 1))
}


if (FALSE) {
    
REPS <- 50000
phi1 <- runif(REPS, -3, 3)
phi2 <- runif(REPS, -3, 3)

causal <- check.causal.AR2(phi1, phi2)

plot(phi1, phi2, col=ifelse(causal, "blue", "red"), pch=16)


}
