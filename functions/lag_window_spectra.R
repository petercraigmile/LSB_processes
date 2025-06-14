

lw.pgram <- function (y, b=0.5 * length(y)^(-1/3)) {
    
    sp  <- pgram(y)
    
    lw <- ksmooth(sp$freq, sp$spec, "normal", bandwidth=b)

    sp$spec <- lw$y

    sp
}
