# CYLINDRICAL K-FUNCTION WITH A FIXED ASPECT RATIO
# spatPat      - point pattern specified using simWind
# alpha        - an angle for which the statistic is computed
# aspect       - aspect of the rectangle; parameter $\zeta$
# r            - vector of distances for which the statistic is computed
# contribution - returns a contribution of each point to the final value of statistic
#              - set TRUE if implementing Marked Point Method

cylK <- function(spatPat, alpha, aspect=NA, r=seq(0,1, length=100),                 
                 contribution=FALSE){
  R <- matrix(c(cos(-alpha), -sin(-alpha), sin(-alpha), cos(-alpha)),
              nrow=2, byrow=TRUE) # rotation matrix by angle -'alpha
  rotPattern <- t(R%*%t(spatPat[[1]]))
  
  lamHat <- nrow(spatPat[[1]])/(spatPat$xDiff*spatPat$yDiff) # intensity estimate
  
  
  xDist <- spatPat$xDiff-as.matrix(dist(cbind(spatPat[[1]][,1],0), "euclidean")) # differences of points and boundaries
  yDist <- spatPat$yDiff-as.matrix(dist(cbind(spatPat[[1]][,2],0), "euclidean"))
  xRotDist <- as.matrix(dist(cbind(rotPattern[,1],0), "euclidean")) # differences between points of the pattern
  yRotDist <- as.matrix(dist(cbind(rotPattern[,2],0), "euclidean"))
  diag(xRotDist) <- NA
  diag(xRotDist) <- NA
  windowSize <- xDist*yDist
  contributions <- matrix(NA, ncol=length(r), nrow= nrow(spatPat[[1]]))
  for(i in 1:length(r)){
    # computes contributions of each point
    contributions[,i] <- apply((xRotDist<=r[i] & yRotDist <=(aspect*r[i]))/windowSize, 1, sum, na.rm=TRUE)/lamHat^2
  }
  # aggregates the contributions
  K <- apply(contributions, 2, sum, na.rm=TRUE)
  if(contribution){
    return(list(K=K, KContrib = contributions))
  } else {
    return(K=K)
  }
  
}
