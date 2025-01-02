# Function to apply a transformation leading to geometric anisotropy (used for LGCP)
# theta          - direction of stretch
# xScalingFactor - a factor $c$ by which the process is stretched along theta (or compressed if 0<=c<1)
# process        - a function simulating an anisotropic point process to be transformed
# xdim, ydim     - observation window bounds in x and y dimensions
# ...            - further arguments ot be passed to process


rAnisotrop <- function(theta, xScalingFactor, process, xdim, ydim, ...){
  R <- matrix(c(cos(theta),-sin(theta), sin(theta), cos(theta)), 
              byrow=TRUE, ncol=2) # rotation matrix
  C <- diag(c(xScalingFactor, 1/xScalingFactor)) # compression-stretch matrix
  Trans <- R%*%C
  window <- simWind(xdim, ydim)
  #extending a simulation window, so that after transformation it will include the requested observation window:
  extrMat <- solve(Trans) %*% t(window$winBor)
  xdim2 <- c(min(extrMat[1,]), max(extrMat[1,]))
  ydim2 <- c(min(extrMat[2,]), max(extrMat[2,]))
  isoPattern <- process(xdim=xdim2, ydim=ydim2, ...)[[1]] # simulate isotropic
  anisoPattern <- t(Trans%*%t(isoPattern)) # apply transofrmation
  inWindow <- (anisoPattern[,1] >= xdim[1] & anisoPattern[,1] <= xdim[2] &
                 anisoPattern[,2] >= ydim[1] & anisoPattern[,2] <= ydim[2]) # truncate to the observation window
  observedAnisoPattern <- anisoPattern[inWindow,]
  return(list(pattern=observedAnisoPattern, winBor=window$winBor, 
              lim=window$lim, xDiff=window$xDiff, yDiff=window$yDiff))
}
