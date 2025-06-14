
library(spectral)

tap.prdgram <- function(y, taper = NULL) {
  ## ======================================================================  
  ## Purpose: To get the tapered periodogram as given in Dahlhaus 1997
  ## with default data taper h = 0.5 * [1 - cos(2*pi*x)] for the time series y
  ## ======================================================================
  
  y <- y - mean(y)
  n <- length(y)
  
  # Data taper
  if (is.null(taper)) {
    
    taper <- 0.5 * (1 - cos(2 * pi * (0:(n - 1)) / n))
    
  }
  
  y <- taper * y
  
  D <- Mod(fft(y)) ^ 2
  
  H <- sum(taper ^ 2)
  
  spectral::spect(D / (2 * pi * H), n)
  
}