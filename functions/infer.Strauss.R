# Inference procedure for Strauss process incl. assessing interaction range
# returns parameters beta, gamma, R, respectively, as per standard notation,

# pattern - point pattern object using simWind specification

infer.Strauss <- function(pattern){
  estimates <- rep(NA, 3)
  xrange <- c(pattern[[2]][1,1], pattern[[2]][3,1])
  yrange <- c(pattern[[2]][2,2], pattern[[2]][3,2])
  L <- Lest(as.ppp(pattern$pattern, owin(xrange=xrange, yrange=yrange)))
  ix <- which(L$theo - L$trans == max (L$theo - L$trans))
  estimates[3] <- L$r[ix]
  estimates[1:2] <- exp(ppm(as.ppp(pattern$pattern, owin(xrange=xrange, yrange=yrange)) ~ 1, Strauss(estimates[3]))$coef)
  return(estimates)
}



