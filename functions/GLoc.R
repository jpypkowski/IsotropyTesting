# LOCAL DIRECTIONAL NEAREST NEIGHBOUR DISTANCE DISTRIBUTION
# spatPat       - point pattern specified using simWind
# alpha         - an angle for which the statistic is computed
# epsilon       - paramter $\epsilon$ of the statistic
# r             - vector of distances for which the statistic is computed
# distance      - optional, input a matrix of distances between points if already available
# angles        - optional, input a matrix of angles between points if already available
# retMatrices   - include matrices distance and angles in the function's output
# contribution  - returns a contribution of each point to the final value of statistic
#               - set TRUE if implementing Marked Point Method


GLoc<- function(spatPat, alpha, epsilon, r=seq(0,1, length=100), 
                distance=NA, angles=NA, retMatrices=FALSE, contribution=FALSE) {
  if(!is.matrix(distance)){
    distance <- (as.matrix(dist(x=spatPat[[1]], method="euclidean"))) # creates a matrix of distances between points
  }
  borders <- cbind(spatPat[[1]][,1]-spatPat$winBor[1,1], spatPat$winBor[2,1]-spatPat[[1]][,1],
                   spatPat[[1]][,2]-spatPat$winBor[1,2], spatPat$winBor[3,2]-spatPat[[1]][,2]) # difference between points and the boundaries of the observation window
  distBorder <- apply(borders, 1, min)
  
  distQual <- matrix(NA, nrow=nrow(distance), ncol=ncol(distance))
  if(!is.matrix(angles)){
    angles <- angleMatrix(spatPat[[1]]) # creates a matrix of angles at which difference vectors between points are directed
  }
  for(i in 1:nrow(spatPat[[1]])){
    TF <- ((angles[i,] > alpha-epsilon & angles[i,] < alpha+epsilon ) |
             (angles[i,] > alpha+pi-epsilon & angles[i,] < alpha+pi+epsilon )) # check if eahc pair of points falls into the infinite double-cone
    distQual[i, TF] <- distance[i, TF]
  }
  diag(distQual) <- NA
  minDist <- apply(distQual, 1, function(vec){ifelse(all(is.na(vec)), NA, min(vec, na.rm=TRUE)) }) # find minimum distances
  
  if(alpha-epsilon>=0 && alpha+epsilon<=pi/2){ # for computation of the estimation correction
    xMar <- 1*cos(alpha-epsilon)
    yMar <- 1*cos(pi/2-alpha-epsilon)
  } else if(alpha-epsilon>=0 && alpha-epsilon<=pi/2  && alpha+epsilon>=pi/2) {
    xMar <- max(1*cos(alpha-epsilon), 1*cos(pi-alpha-epsilon))
    yMar <- 1
  } else if(alpha-epsilon>=pi/2  && alpha+epsilon<=pi){
    xMar <- 1*cos(pi-alpha-epsilon)
    yMar <- 2*cos(alpha-epsilon-pi/2)
  } else {
    xMar <- 1
    yMar <- max(abs(1*sin(alpha+epsilon)), 1*sin(pi-alpha-epsilon))
  }
  volumes <- (spatPat$xDiff-2*xMar*minDist)*(spatPat$yDiff-2*xMar*minDist) 
  NAs <- which(is.na(volumes)|is.na(minDist))
  if(length(NAs)>0){
    volumes2 <- volumes[-NAs]
    minDist2 <- minDist[-NAs]
    spatPat2 <- spatPat[[1]][-NAs,]
  } else{
    volumes2 <- volumes
    minDist2 <- minDist
    spatPat2 <- spatPat[[1]]
  }
  contributions <- matrix(NA, ncol=length(r), nrow=nrow(spatPat[[1]]))
  border <- which(spatPat[[1]][,1]>spatPat[[2]][1,1]+xMar*minDist & spatPat[[1]][,1]<spatPat[[2]][3,1]-xMar*minDist &
                    spatPat[[1]][,2]>spatPat[[2]][1,2]+yMar*minDist & spatPat[[1]][,2]<spatPat[[2]][3,2]-yMar*minDist) # 1[d_i < r]
  if(length(border)>0){
    for(i in 1:length(r)){
      # each point's contribution at each r
      contributions[border,i] <- (minDist[border]<r[i])/volumes[border] 
      contributions[-border,i] <- 0
      contributions[-NAs,i]
      
    }
    numerator <- apply(contributions, 2, sum, na.rm=TRUE) # aggregate
    denominatorContirb <- rep(0, nrow(spatPat[[1]]))
    denominatorContirb[border] <- (minDist[border]<sqrt(spatPat$xDiff^2+spatPat$yDiff^2))/volumes[border]
    denominator <- sum(denominatorContirb)
    if(retMatrices && contribution){
      return(list(Gloc=numerator/denominator, distance=distance, angles=angles, 
                  GContrib = list(numContrib = contributions, denContrib = denominatorContirb)))
    } else if(retMatrices){
      return(list(Gloc=numerator/denominator, distance=distance, angles=angles))
    } else if(contribution){
      return(list(Gloc=numerator/denominator, 
                  GContrib = list(numContrib = contributions, denContrib = denominatorContirb)))
    } else {
      return(numerator/denominator)
      
    }
  } else {
    if(retMatrices && contribution){
      return(list(Gloc=rep(0, length(r)), distance=distance, angles=angles, 
                  GContrib = list(numContrib = matrix(0, ncol=length(r), nrow=nrow(spatPat[[1]])), 
                                  denContrib = matrix(1, ncol=length(r), nrow=nrow(spatPat[[1]])))))
    } else if(retMatrices){
      return(list(Gloc=rep(0, length(r)), distance=distance, angles=angles))
    } else if(contribution){
      return(list(Gloc=rep(0, length(r)), 
                  GContrib = list(numContrib = matrix(0, ncol=length(r), nrow=nrow(spatPat[[1]])), 
                                  denContrib = matrix(1, ncol=length(r), nrow=nrow(spatPat[[1]])))))
    } else {
      return(rep(0, length(r)))
      
    }
  }
  
  
}
