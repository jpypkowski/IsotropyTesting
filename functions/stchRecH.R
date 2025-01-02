# STOCHASTIC RECONSTRUCTION ALGORITHM WITH H-FUNCTION
# spatPat - point pattern specified using simWind
# rmax    - optional; a pre-set upper integration bound for the H-function
# rlength - length of the sequence of distances up to rmax; has to be used if rmax is specified
# iter    - the number of iterations for which the algorithm is run



stchRecH <- function(spatPat, rmax=NULL, rlength=NULL, iter=1e4){
  
  xrange <- c(spatPat[[2]][1,1], spatPat[[2]][3,1])
  yrange <- c(spatPat[[2]][2,2], spatPat[[2]][3,2])
  
  # SET UP INTEGRATION GRID FOR R
  if(!is.null(rmax) & !is.null(rlength)){
    interval <- rmax/rlength
    rs <- seq(0, rmax, by=interval)
    
    originH <- Hest(X=as.ppp(spatPat[[1]], W=owin(xrange=xrange, yrange=yrange)), r = rs, correction = "raw", conditional=FALSE)$raw
  } else {
    # set automatically by spatstat so that H(rmax) = 1
    primH <- Hest(X=as.ppp(spatPat[[1]], W=owin(xrange=xrange, yrange=yrange)), correction = "raw", conditional=FALSE)
    originH <- primH$raw
    rs <- primH$r
    interval <- mean(diff(rs))
  }
  N <- nrow(spatPat[[1]])
  
  
  
  # INITIALISATION - BINOMIAL PROCESS with the same number of points as the original pattern
  curPat <- rHomPois(xdim=xrange, ydim=yrange, numPoints = N)
  curH <- Hest(X=as.ppp(curPat[[1]], W=owin(xrange=xrange, yrange=yrange)), r = rs, correction = "raw", conditional=FALSE)$raw
  curD <- sum(interval*(curH-originH)^2)
  for(i in 1:iter){
    move <- sample(1:N, 1) # sampled index of a point to be moved
    propPat <- curPat
    propPat[[1]][move,] <- c(runif(1, xrange[1], xrange[2]), runif(1, yrange[1], yrange[2])) # random new location
    propH <- Hest(X=as.ppp(propPat[[1]], W=owin(xrange=xrange, yrange=yrange)), r = rs, correction = "raw", conditional=FALSE)$raw
    propD <- sum(interval*(propH-originH)^2) # deviation for the H-function
    if(curD >= propD){
      accept <- 1 # accept if deviance improved
    } else {
      TS <- (0.1- 0.1*(i/(iter+1))^(1/10)) #low-temperature scaling temperatures
      Paccept <- exp((propD-curD)/TS)/100 # probability of acceptance if deviance isn't improved by the proposal
      accept <- runif(1, 0, 1)<Paccept
    }
    if(accept == 1){
      curPat <- propPat
      curH <- propH
      curD <- propD
    }
  }
  return(curPat)
}
