# HOMOGENOUS POISSON POINT PROCES SIMULATION
# intensity  - intensity of the process (an expected number of points per unit area)
# xdim, ydim - bounds of the observation window in x and y dimensions
# numPoints  - a pre-specified number of points; if set, the function simulates a binomial process
rHomPois <- function(intensity=NA, xdim, ydim, numPoints=NA){
  window <-simWind(xdim, ydim)
  winSize <- window[[1]]
  if(is.na(numPoints)) numPoints <- max(1, rpois(1, intensity*winSize))
  x <- runif(numPoints, min=xdim[1], max=xdim[2])
  y <- runif(numPoints, min=ydim[1], max=ydim[2])
  return(list(pattern = cbind(x, y), winBor=window[[2]],
              lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff))
}
