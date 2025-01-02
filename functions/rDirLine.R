# POISSON LINE CLUSTER PROCESS SIMULATION

# lineIntensity - intensity of the underlying Poisson line process
# alongIntensity - intensity of the one-dimensional Poisson process on each latent line
# theta          - preferred direction, the location parameter of the von Mises distribution
# xdim, ydim     - observation windows bounds in x, y dimensions, respectively
# sigma          - standard deviation of the normal distribution according to which points are
#                  dispersed around the latent lines
# kappa          - concentration parameter, often expressed as KofA(a), where a is the anisotropy 
#                  parameter (see below)
# margin         - additional margin by which the observation window is extended for the simulation

rDirLine <- function(lineIntensity, alongIntensity, theta, xdim, ydim, sigma, kappa, margin=0.05){
  window <- simWind(xdim, ydim)
  
  xDiff <- window$xDiff+2*margin
  yDiff <- window$yDiff+2*margin
  r <- sqrt(xDiff^2+yDiff^2)/2
  noLines <- max(1, rpois(1, 2*r*lineIntensity)) # the number of lines on a circumscribed circle
  intersections <- runif(noLines, -r, r)
  thetas <- as.numeric(rvonmises(n=noLines, mu=theta, kappa=kappa)) # angles of the lines
  observedLine <- matrix(NA, nrow=1, ncol=2)
  for(i in 1:noLines){ # for each line
    lineNumPoints <- rpois(1, alongIntensity*2*r) #number of points on a line of maximum length
    linePointsAcross <- rnorm(lineNumPoints, 0, sigma) # distance of points from the line
    linePointsAlong <- runif(lineNumPoints, -r, r)     # location along the lines
    linePoints <- cbind(linePointsAlong, linePointsAcross+intersections[i]) # shifting the line with points
    R <- matrix(c(cos(thetas[i]),-sin(thetas[i]), sin(thetas[i]), cos(thetas[i])), 
                byrow=TRUE, ncol=2)
    rotated <- t(R%*%t(linePoints)) # rotating the shifted line with points
    observedLine <- rbind(observedLine, cbind(rotated[,1], rotated[,2]))
  }
  if(nrow(observedLine)>2){
    observedLine <- observedLine[-1,]
  } else {
    observedLine <- matrix(observedLine[-1,], ncol=2)
  }
  
  inWindow <- (observedLine[,1] >= xdim[1] & observedLine[,1] <= xdim[2] &
                 observedLine[,2] >= ydim[1] & observedLine[,2] <= ydim[2])
  observedPoints <- observedLine[inWindow,] # trimming to observation window
  
  
  if(sum(inWindow==1)){
    return(list(pattern=matrix(observedPoints, ncol=2), winBor=window$winBor, 
                lim=window$lim,  xDiff=window$xDiff, yDiff=window$yDiff, intersections = intersections, theta =theta))
  } else {
    return(list(pattern=observedPoints, winBor=window$winBor, 
                lim=window$lim,  xDiff=window$xDiff, yDiff=window$yDiff, intersections = intersections, theta =theta))
  }
  
}

# AUXILIARY FUNCTION THAT CONVERTS ANISOTROPY PARAMETER a TO kappa

KofA <- function(a) 5*(1-exp(1-1/a))
