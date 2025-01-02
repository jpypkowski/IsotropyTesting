# A function used to simulate bigger Gibbs patterns faster by splitting the observation window
# Uses BirthDeathMove in each of the subregions
# sqrt.sqlit      - square root of the number of subregions into which the observation window is split
# M_small         - number of iterations for which BirthDeathMove is run in each subergion
# M_big           - number of iterations for which BirthDeathMove is run for the full observation window
#                 - starting from the pattern consiting of combined subregions
# addMargin_small - a distance by which each subregion is extended for the simulation
# addMargin_big   - a distance by which the full pbservation window is extended for the simulation
# intn_small      - a vector of nunbers of points in the initial patterns for each extended subregion

BDM_BIG <- function(sqrt.split=2, interactionF, p, q, M_small, M_big, xdim=c(0,1), ydim=c(0,1), 
                    addMargin_small = 0.05, addMargin_big = 0.05, ..., initn_small){
  xdim_old <- xdim
  ydim_old <- ydim
  xdim <- xdim+c(-addMargin_big, addMargin_big)
  ydim <- ydim+c(-addMargin_big, addMargin_big)
  xDist <- diff(xdim)/2/sqrt.split
  xMids <- seq(xdim[1]+xDist, xdim[2]-xDist, length=sqrt.split)
  yDist <- diff(ydim)/2/sqrt.split
  yMids <- seq(ydim[1]+yDist, ydim[2]-yDist, length=sqrt.split)
  mids <- expand.grid(xMids, yMids)
  pats <- vector(mode="list", length=sqrt.split^2)
  lengths <- c(0, rep(NA, sqrt.split^2))
  for(i in 1:sqrt.split^2){
    pats[[i]] <- BirthDeathMove(interactionF, p, q, M_small, xdim=c(-xDist,xDist), ydim=c(-yDist,yDist), addMargin = addMargin_small, 
                                ..., initn=initn_small[[i]])[[1]][[1]]
    if(!is.null(nrow(pats[[i]]))){
      lengths[i+1] <- nrow(pats[[i]])
    } else lengths[i+1] <- 0
    
  }
  auxPat <- matrix(NA, ncol=2, nrow=sum(lengths))
  for(i in 1:sqrt.split^2){
    if(lengths[i+1]>0){
      auxPat[(1+sum(lengths[1:i])):sum(lengths[1:(i+1)]),] <- cbind(pats[[i]][,1]+mids[i,1], pats[[i]][,2]+mids[i,2])
      
    }
  }
  window <- simWind(xdim_old, ydim_old)
  
  finalPat <- BirthDeathMove(interactionF, p, q, M_big, xdim=xdim_old, ydim=ydim_old, 
                             addMargin = addMargin_big, ..., initPat=list(pattern=auxPat, winBor=window[[2]], lim=window[[3]],  
                                                                          xDiff=window[[4]], yDiff=window[[5]]))[[1]]
  return(finalPat)
}
