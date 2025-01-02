# THOMAS PROCESS SIMULATION FOR BOTH ISOTROPIC AND ANISOTROPIC CASE
# parentIntensity - intensity of the underlying homogeneous Poisson process generating unobserved parents
# meanClsuterSize - average number of offsprings generated from parent points
# sigma           - standard deviation of the (isotropic) normal distribution describing the spread of the offsprings 
#                   around their parent
# theta           - rotation angle for the geometric transformation
# xScalingFactor  - stretch/compression factor along the angle theta, setting it to 1 yields an isotropic point pattern
# xdim, ydim      - observation window bounds (see simWind function)
# addMargin       - used for extending the simulation window by addMargin*sigma in each direction relative to the observation window
#                   this is done to avoid distortions near observation window borders


rThomas <- function(parentIntensity, meanClusterSize, 
                    sigma, theta=0, xScalingFactor=1, xdim, ydim, addMargin=4){
  
  xdim2 <- xdim+c(-addMargin*sigma, addMargin*sigma)
  ydim2 <- ydim+c(-addMargin*sigma, addMargin*sigma)
  trueWindow <- simWind(xdim, ydim)
  parents <- rHomPois(parentIntensity, xdim2, ydim2) # generate parent points
  clusterSizes <- rpois(nrow(parents[[1]]), meanClusterSize) # generate Poisson numbers of offsprings for each parent
  offsprings <- matrix(NA, nrow=sum(clusterSizes), ncol=2)
  if(xScalingFactor==1){  # isotropic case
    for(i in 1:nrow(parents[[1]])){
      stoppedAt <- max(c(0,which(!is.na(offsprings[,1]))))
      cs <- clusterSizes[i]
      if(cs>0) offsprings[(stoppedAt+1):(stoppedAt+cs),]  <- rmvnorm(cs, # generate offsprings' locations
                                                                     parents[[1]][i,], diag(sigma^2, nrow=2))
    }
  } else { # anisotropic case
    R <- matrix(c(cos(theta),-sin(theta), sin(theta), cos(theta)), 
                byrow=TRUE, ncol=2) # rotation matrix
    C <- diag(c(xScalingFactor, 1/xScalingFactor)) # stretch-compression matrix
    Trans <- R%*%C
    for(i in 1:nrow(parents[[1]])){
      stoppedAt <- max(c(0,which(!is.na(offsprings[,1]))))
      cs <- clusterSizes[i]
      if(cs>0) offsprings[(stoppedAt+1):(stoppedAt+cs),]  <- rmvnorm(cs, # generate offsprings' locations
                                                                     parents[[1]][i,], Trans%*%diag(sigma^2, nrow=2)%*%t(Trans))
    }
  }
  
  inWindow <- (offsprings[,1] >= xdim[1] & offsprings[,1] <= xdim[2] &
                 offsprings[,2] >= ydim[1] & offsprings[,2] <= ydim[2]) 
  observedOffsprings <- offsprings[inWindow,] # truncate to the observation window
  return(list(pattern=observedOffsprings, winBor=trueWindow$winBor, 
              lim=trueWindow$lim, xDiff=trueWindow$xDiff, yDiff=trueWindow$yDiff))
}

