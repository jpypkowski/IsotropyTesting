# FUNCTION PRODUCING A SINGLE PATTERN REPLICATE USING TILING (for square observation windows)  
# spatPat - point pattern specified using simWind
# sqrt.n  - suqare root of the number of tiles $k$, s.t. $N_{tile}=k^2$

singleBoot <- function(spatPat, sqrt.n){
  n <- sqrt.n^2
  xDiff <- spatPat$xDiff
  yDiff <- spatPat$yDiff
  pattern <- rbind(t(spatPat$pattern), 1:nrow(spatPat$pattern))
  if(yDiff!=xDiff) stop("square pls")
  r <- xDiff*sqrt(2)/sqrt.n/2
  # set out the grid of centres of subregions of the reconstructed pattern
  x.centres <- seq(spatPat$winBor[1,1]+r, spatPat$winBor[3,1]-r, length=sqrt.n)
  y.centres <- seq(spatPat$winBor[1,2]+r, spatPat$winBor[3,2]-r, length=sqrt.n)
  grid.centres <- expand.grid(x.centres, y.centres)
  if(nrow(grid.centres)!=n) stop("num. of points in the centres (reconstrutcted) grid not correct")
  
  angles <- runif(n, min=0, max=2*pi)
  
  # set out the grid of centres of tiles (in the observed pattern)
  x.tiles <- seq(spatPat$winBor[1,1]+xDiff/sqrt.n/2, spatPat$winBor[3,1]-xDiff/sqrt.n/2, length=sqrt.n)
  y.tiles <- seq(spatPat$winBor[1,2]+yDiff/sqrt.n/2, spatPat$winBor[3,2]-yDiff/sqrt.n/2, length=sqrt.n)
  grid.tiles <- expand.grid(x.tiles, y.tiles)
  if(nrow(grid.tiles)!=n) stop("num. of points in the tiles (original) grid not correct")
  
  # SAVING TRANSLATED (but not rotated) TILES
  transTiles <- vector(mode="list", length = n)
  for(i in 1:n){
    aux.pattern <- pattern
    aux.pattern[1,] <- pattern[1,]-grid.centres[i,1]
    aux.pattern[2,] <- pattern[2,]-grid.centres[i,2] # translate the whole pattern
    distances <- sqrt(aux.pattern[1,]^2+aux.pattern[2,]^2) # find distances from the origin
    trans <- aux.pattern[,distances <= r] # select points within required distance
    if(is.matrix(trans)) {
      transTiles[[i]] <- trans
    } else {
      transTiles[[i]] <- as.matrix(trans)
    }
    
  }
  
  
  
  which.tiles <- sample(1:n, n, replace=TRUE)
  
  new.pattern <- vector(mode="list", length=n)
  sampled.points <- vector(mode="list", length=n)
  lengths <- rep(NA, n)
  for(i in 1:n){
    # put it all together
    tile <- which.tiles[i] # sample a tile
    aux.pattern <- transTiles[[tile]]
    
    angle <- angles[i] # sample a random angle
    R <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),ncol=2, byrow=TRUE)
    rotated <- rbind(R%*%(as.matrix(aux.pattern[-3,])), aux.pattern[3,]) # rotate
    rot.sq <- as.matrix(rotated[,(rotated[1,] <= xDiff/sqrt.n/2 & rotated[1,] >= -xDiff/sqrt.n/2 & 
                                    rotated[2,] <= xDiff/sqrt.n/2 & rotated[2,] >= -xDiff/sqrt.n/2)]) # truncate to a square
    new.pattern[[i]] <- cbind(rot.sq[1,]+grid.tiles[i,1], rot.sq[2,]+grid.tiles[i,2]) # translate by the centre of the relevant subregion
    sampled.points[[i]] <- rot.sq[3,] 
    lengths[i] <- nrow(new.pattern[[i]])
  }
  reconstr <- matrix(NA, nrow=sum(lengths), ncol=2)
  for(i in 1:n){
    # put all the tiles together
    if(i==1 && lengths[i]>0) {
      reconstr[1:lengths[1],] <- new.pattern[[1]]
    } else {
      if(lengths[i]>0) {
        reconstr[(1+sum(lengths[1:(i-1)])):sum(lengths[1:i]),] <- new.pattern[[i]]
      }
      
      
    }
    
    
  }
    
    

  return(list(pattern = reconstr, winBor=spatPat$winBor,
              lim=spatPat$lim, xDiff=xDiff, yDiff=yDiff))
}