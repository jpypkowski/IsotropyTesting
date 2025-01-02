# Function for estimating parameters of LGCP with Matern
# returns a vector of estimates variance, scale, and mu, respectively
# spatPat    - point pattern specified using simWind
# covariance - assumed covariance function
# startpar   - initial parameters for the search 
# ...        - further arguments to control covariance model



infer.LGCP <- function(spatPat, covariance, startpar=c(var=1, scale=0.8), ...){
  estimates <- rep(NA, 3) # variance, scale, mu
  lambda <- nrow(spatPat$pattern)/(spatPat$xDiff*spatPat$yDiff)
  estimates[1:2] <- lgcp.estpcf(as.ppp(spatPat$pattern, 
                                       W=owin(c(spatPat$winBor[1,1], spatPat$winBor[2,1]), 
                                              c(spatPat$winBor[1,2], spatPat$winBor[3,2]))),
                                startpar=startpar, covmodel = list(model=covariance, ...), lambda=lambda)$par
  estimates[3] <- log(lambda)-estimates[1]/2
  return(estimates)
  
}