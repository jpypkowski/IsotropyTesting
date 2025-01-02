# BIRTH DEATH MOVE ALGORITHM FOR SIMULAITON OF GIBBS PROCESSS
# interactionF  - pairwise interaction function, see fAnisoLennardJones for an example
# p             - probability that a new birth will be proposed at one step of the algorithm
# q             - probability that a move will be proposed at one step of the algorithm
# M             - number of iterations for which the algorithm is run
# xdim, ydim    - bounds of the observation window in x and y dimensions
# adddMargin    - a distance by which the observation window is extended for the simulation
# ...           - additional arguments passed to interactionF
# initN         - pre-specified number of points of the initial pattern simulated from a binomial process
# initIntensity - intensity of the homogeneous Poisson process used to simulate the initial pattern
# initPat       - initial pattern specified directly

BirthDeathMove <- function(interactionF, p, q, M, xdim=c(0,1), ydim=c(0,1), addMargin = 0.05, ..., initn=NULL, initIntensity=NULL, initPat=NULL){
  
  
  xdim_old <- xdim
  ydim_old <- ydim
  xdim <- xdim+c(-addMargin, addMargin)
  ydim <- ydim+c(-addMargin, addMargin)
  
  if(is.null(initn) & is.null(initIntensity) & is.null(initPat)) stop("Define initn, initIntensity, or initPat")
  # simulate the initial pattern if not specified
  if(!is.null(initPat)){
    curPat <- initPat
  } else  if(is.null(initn)) {
    curPat <- rHomPois(intensity = initIntensity, xdim=xdim, ydim=ydim)
  } else {
    curPat <- rHomPois(numPoints = initn, xdim=xdim, ydim=ydim)
  }
  densities <- rep(NA, M+1)
  densities[1] <- interactionF(spatPat=curPat, log=TRUE, ...)
  
  for(m in 1:M){
    print(m)
    Rm3 <- runif(1)
    n <- nrow(curPat$pattern)
    if (Rm3 <= q && n>0) {
      # MOVE
      I <- sample(1:n, 1) # index of point to be moved
      Rm <- runif(1)
      ksi <- c(runif(1, xdim[1], xdim[2]), runif(1, ydim[1], ydim[2])) # proposed new location
      propPat <- curPat
      propPat$pattern[I,] <- ksi # pattern with the proposed moved point
      r <- exp(interactionF(spatPat=propPat, log=TRUE, ...)-interactionF(spatPat=curPat, log=TRUE, ...)) # acceptance probability
      if(Rm <= r){
        curPat <- propPat
      } 
    } else {
      Rm1 <- runif(1)
      Rm2 <- runif(1)
      ########### ADAPT
      propPat <- curPat
      if(Rm1 <= p){
        # BIRTH
        ksi <- c(runif(1, xdim[1], xdim[2]), runif(1, ydim[1], ydim[2])) # proposed new location
        propPat$pattern <- rbind(propPat$pattern, ksi) # pattern with the proposed new point
        r <- exp(interactionF(spatPat=propPat, log=TRUE, ...)-interactionF(spatPat=curPat, log=TRUE, ...)) # acceptance probability
        if(Rm <= r){
        if(Rm2 <= r){
          curPat <- propPat
        } 
      } else {
        # DEATH
        if (n>0) {
          eta <- sample(1:n, 1)
          if(n==2){
            propPat$pattern <- matrix(propPat$pattern[-eta,], nrow=1)
          } else if(n==1){
            propPat$pattern <- propPat$pattern[-eta,]
          } else
            propPat$pattern <- propPat$pattern[-eta,]
          r <- exp(interactionF(spatPat=propPat, log=TRUE, ...)-interactionF(spatPat=curPat, log=TRUE,...))
          if(Rm2 <= r){
            curPat <- propPat
          } 
        } 
        
      }
      
    }
    densities[m+1] <- interactionF(spatPat=curPat, log=TRUE, ...) # tracking densities 
    
    # uncomment if to print patterns as the algorithm progresses (set to every 500 iterations)
    # if(m%%500 == 0) {
    #   plotSpatPat(curPat, ylab=nrow(curPat[[1]])/(diff(xdim)*diff(ydim)))
    # }
  }
  window <- simWind(xdim_old, ydim_old)
  finalPat <- list(pattern=curPat$pattern[(curPat$pattern[,1]>=xdim_old[1] & curPat$pattern[,1]<=xdim_old[2] & curPat$pattern[,2]>=ydim_old[1] & curPat$pattern[,2]<=ydim_old[2]),],
                   winBor=window[[2]], lim=window[[3]],  xDiff=window[[4]], yDiff=window[[5]])
  return(list(spatPat = finalPat, densities = densities))
}

