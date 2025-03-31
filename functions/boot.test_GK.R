boot.test_GK <- function(spatPat, sqrt.n, reps, r.grid.W, 
                          prefDir, Gepsilon, Kaspect){
  
  
  
  # setting up tiles and subregions
  {
    n <- sqrt.n^2
    xDiff <- spatPat$xDiff
    yDiff <- spatPat$yDiff
    pattern <- t(spatPat$pattern)
    if(yDiff!=xDiff) stop("square pls")
    r <- xDiff*sqrt(2)/sqrt.n/2
    x.centres <- seq(spatPat$winBor[1,1]+r, spatPat$winBor[3,1]-r, length=sqrt.n)
    y.centres <- seq(spatPat$winBor[1,2]+r, spatPat$winBor[3,2]-r, length=sqrt.n)
    grid.centres <- expand.grid(x.centres, y.centres)
    if(nrow(grid.centres)!=n) stop("num. of points in the centres (reconstrutcted) grid not correct")
    x.tiles <- seq(spatPat$winBor[1,1]+xDiff/sqrt.n/2, spatPat$winBor[3,1]-xDiff/sqrt.n/2, length=sqrt.n)
    y.tiles <- seq(spatPat$winBor[1,2]+yDiff/sqrt.n/2, spatPat$winBor[3,2]-yDiff/sqrt.n/2, length=sqrt.n)
    grid.tiles <- expand.grid(x.tiles, y.tiles)
    if(nrow(grid.tiles)!=n) stop("num. of points in the tiles (original) grid not correct")
  }
  
  # SAVING TRANSLATED (but not rotated) TILES
  transTiles <- vector(mode="list", length = n)
  
  for(i in 1:n){
    aux.pattern <- pattern
    aux.pattern[1,] <- pattern[1,]-grid.centres[i,1]
    aux.pattern[2,] <- pattern[2,]-grid.centres[i,2]
    distances <- sqrt(aux.pattern[1,]^2+aux.pattern[2,]^2)
    trans <- aux.pattern[,distances <= r]
    if(is.matrix(trans)) {
      transTiles[[i]] <- trans
    } else {
      transTiles[[i]] <- as.matrix(trans)
    }
    
  }
  
  vGtile.W <- matrix(NA, nrow= length(r.grid.W), ncol= reps)
  vKtile.W <- matrix(NA, nrow= length(r.grid.W), ncol= reps)

  
  
  
 
  
  for(m in 1:reps){
    # sampling angles and tiles
    
    
    # extracting sampled subregions and points
    new.pattern <- vector(mode="list", length=n)
    lengths <- rep(0, n)
    while(sum(lengths) < 2){
      for(i in 1:n){
        angles <- runif(n, min=0, max=2*pi)
        which.tiles <- sample(1:n, n, replace=TRUE)
        tile <- which.tiles[i]
        aux.pattern <- transTiles[[tile]]
        
        angle <- angles[i]
        R <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),ncol=2, byrow=TRUE)
        rotated <- R%*%(as.matrix(aux.pattern))
        rot.sq <- as.matrix(rotated[,(rotated[1,] <= xDiff/sqrt.n/2 & rotated[1,] >= -xDiff/sqrt.n/2 & 
                                        rotated[2,] <= xDiff/sqrt.n/2 & rotated[2,] >= -xDiff/sqrt.n/2)])
        new.pattern[[i]] <- cbind(rot.sq[1,]+grid.tiles[i,1], rot.sq[2,]+grid.tiles[i,2])
        lengths[i] <- nrow(new.pattern[[i]])
      }
    }
    
    
    #putting everything together
    reconstr <- matrix(NA, nrow=sum(lengths), ncol=2)
    for(i in 1:n){
      if(i==1 && lengths[i]>0) {
        reconstr[1:lengths[1],] <- new.pattern[[1]]
      } else {
        if(lengths[i]>0) {
          reconstr[(1+sum(lengths[1:(i-1)])):sum(lengths[1:i]),] <- new.pattern[[i]]
        }
        
      }
      
      
    }
    recPat <- list(pattern = reconstr, winBor=spatPat$winBor,
                   lim=spatPat$lim, xDiff=xDiff, yDiff=yDiff)
    
    # GET VALUES OF VECTORS OF (contrasts of) TEST STATISTICS 
    
    secDir <- (prefDir+pi/2)%%pi
    
    # *** TILING *** #
    
    #G
    prefG <- GLoc(spatPat=recPat, alpha=prefDir, epsilon=Gepsilon, r=r.grid.W, retMatrices = TRUE)
    secG <- GLoc(spatPat=recPat, alpha=secDir, epsilon=Gepsilon, r=r.grid.W, 
                 distance=prefG$distance, angles=prefG$angles)
    vGtile.W[,m] <- prefG$Gloc - secG
    
    
    
    
    #K
    
    prefK <- cylK(spatPat=recPat, alpha=prefDir, aspect=Kaspect, r=r.grid.W)
    secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid.W)
    vKtile.W[,m] <- prefK-secK
   
    
    
  }
  
  return(list(vGtile.W=vGtile.W, vKtile.W=vKtile.W))
}
