# PERIODOGRAM 
# spatPat      - point pattern specified using simWind
# p            - integers series $p_1$ and $p_2$ (assumed to be equal)
# contribution - returns a contribution of each point to the total periodogram
#              - keep FALSE for a faster computation, set TRUE only if implementing Marked Point Method

periodogram <- function(spatPat, p, contribution=FALSE){
  p1 <- p
  p2 <- p
  omega1 <- 2*pi*p1/spatPat$xDiff
  omega2 <- 2*pi*p2/spatPat$yDiff
  grid <- as.matrix(expand.grid(omega1, omega2))
  area <- spatPat$xDiff*spatPat$yDiff
  if(contribution){
    # compute contributions of each point
    contrib <- matrix(NA, ncol=nrow(grid), nrow=nrow(spatPat$pattern))
    contrib.og <- exp(-1i*grid%*%t(spatPat[[1]]))/sqrt(area)
    contrib.cc <- Re(contrib.og)-1i*Im(contrib.og)
    for(i in 1:nrow(grid)){
      out <- outer(contrib.og[i,], contrib.cc[i,], "*")
      contrib[,i] <- Re(apply(out,1,sum)+apply(out,2,sum))/2
    }
    sdf <- apply(contrib, 2, sum, na.rm=TRUE)
    
  } else {
    # compute for an entire pattern at once
    DFT <- rowSums(exp(-1i*grid%*%t(spatPat[[1]])))/sqrt(area)
    sdf <- Re(DFT)^2 + Im(DFT)^2
  }
  # exclude zero
  zero1 <- which(omega1==0)
  zero2 <- which(omega2==0)
  matrixSDF <- matrix(NA, nrow=length(omega1), ncol=length(omega2))
  angles <- matrix(NA, nrow=length(omega1), ncol=length(omega2))
  ang <- atan((grid[,2]/grid[,1]))
  for(i in 1:length(omega1)){
    matrixSDF[i,] <- sdf[((i-1)*length(omega2)+1):((i)*length(omega2))] # put in a matrix form instead of vector
    angles[i,] <- ang[((i-1)*length(omega2)+1):((i)*length(omega2))] # compute angles corresponding to different frequencies
  }
  matrixSDF[zero1,zero2] <- NA
  angles[zero1,zero1] <- NA
  angles <- angles %% pi
  
  
  if(contribution){
    return(list(sdf = matrixSDF, omega1=omega1, omega2=omega2, p1, p2, angles=angles,
                TContrib = contrib))
  } else {
    return(list(sdf = matrixSDF, omega1=omega1, omega2=omega2, p1, p2, angles=angles))
  }
 
}
