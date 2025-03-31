# INFERENCE FOR THOMAS PROCESS



# NEW

infer.Thomas <- function(pattern, startpar=c(kappa=10,scale=0.1)){
  estimates <- rep(NA, 3)
  xrange <- c(pattern[[2]][1,1], pattern[[2]][3,1])
  yrange <- c(pattern[[2]][2,2], pattern[[2]][3,2])
  estimates[1:2] <- thomas.estpcf(as.ppp(pattern$pattern, owin(xrange=xrange, yrange=yrange)), startpar=startpar)$par
  estimates[2] <- sqrt(estimates[2])
  estimates[3] <- nrow(pattern$pattern)/estimates[1]/(diff(xrange)*diff(yrange))
  return(estimates)
}


# OLD

infer.Thomas <- function(pattern, startpar=c(kappa=10,scale=0.1)){
  estimates <- rep(NA, 3)
  xrange <- c(pattern[[2]][1,1], pattern[[2]][3,1])
  yrange <- c(pattern[[2]][2,2], pattern[[2]][3,2])
  estimates[1:2] <- thomas.estK(as.ppp(pattern$pattern, owin(xrange=xrange, yrange=yrange)), startpar=startpar)$par
  estimates[2] <- sqrt(estimates[2])
  estimates[3] <- nrow(pattern$pattern)/estimates[1]/(diff(xrange)*diff(yrange))
  return(estimates)
}
# parent intensity, sigma, avg no of offsprings

# CORRECT