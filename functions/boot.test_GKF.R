# VALUES FOR OBTAINING TEST STATISTICS USING TILING
# Includes G, K, and Theta spectrum
# This function was designed for the simulation study and produces all replicates for all DSSs
# in one go without parallelisation; it is inefficient for a single pattern and one DSS
# spatPat    - point pattern specified using simWind
# sqrt.n     - square root of the number of tiles, i.e. $N_{tile}=sqrt.n^2$
# reps       - number of replicates to be produced
# r.grid     - a series of distances for G and K taking distance as an argument
# theta.grid - a series of angles for Theta spectrum taking angle as an argument
# prefDir    - tested direction $\alpha_1$
# Gepsilon   - parameter $\epsiolon$ of the G-function
# Kaspect    - parameter $\zeta$ of the cylindrical K-function
# Tepsilon   - bandwidth $h$ of the Theta spectrum
# Tp         - series of in integers, $p_1$ and $p_2

boot.test_GKF <- function(spatPat, sqrt.n, reps, r.grid, theta.grid, 
                      prefDir, Gepsilon, Kaspect, Tepsilon, Tp){
  
  
  
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
  
  vGtile <- matrix(NA, nrow= length(r.grid), ncol= reps)
  vKtile <- matrix(NA, nrow= length(r.grid), ncol= reps)
  vTtile <- matrix(NA, nrow= length(theta.grid), ncol= reps)


  

  omega1 <- omega2 <- 2*pi*Tp/spatPat$xDiff
  grid <- as.matrix(expand.grid(omega1, omega2))
  zero <- which(grid[,1]==0 & grid[,2]==0)
  angles <- atan(expand.grid(Tp,Tp)[,2]/expand.grid(Tp,Tp)[,1])
  angles[zero] <- NA
  
  for(m in 1:reps){
    
    # sampling angles and tiles
    angles <- runif(n, min=0, max=2*pi)
    which.tiles <- sample(1:n, n, replace=TRUE)
    
    # extracting sampled subregions and points
    new.pattern <- vector(mode="list", length=n)
    lengths <- rep(NA, n)
    for(i in 1:n){
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
    prefG <- GLoc(spatPat=recPat, alpha=prefDir, epsilon=Gepsilon, r=r.grid, retMatrices = TRUE)
    secG <- GLoc(spatPat=recPat, alpha=secDir, epsilon=Gepsilon, r=r.grid, 
                 distance=prefG$distance, angles=prefG$angles)
    vGtile[,m] <- prefG$Gloc - secG
    
    #K
    
    prefK <- cylK(spatPat=recPat, alpha=prefDir, aspect=Kaspect, r=r.grid)
    secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid)
    vKtile[,m] <- prefK-secK
    
    #Theta spectrum
    vTtile[,m] <- thetaSpectrum(prdgrm = periodogram(spatPat=recPat, p=Tp), 
                                thetas = theta.grid, epsilon = Tepsilon)
  }
  
  return(list(vGtile=vGtile, vKtile=vKtile, vTtile=vTtile))
}


# THIS VERSION OF THE FUNCTION CAN BE USED IF ONE WISHES TO TREY TWO DIFFERENT SEQUENCES OF DISTANCESFOR 2 INTEGRATION RANGES:

boot.test_GKF_2 <- function(spatPat, sqrt.n, reps, r.grid.W, r.grid.L, theta.grid, 
                          prefDir, Gepsilon, Kaspect, Tepsilon, Tp){
  
  
  
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
  vGtile.L <- matrix(NA, nrow= length(r.grid.L), ncol= reps)
  vKtile.W <- matrix(NA, nrow= length(r.grid.W), ncol= reps)
  vKtile.L <- matrix(NA, nrow= length(r.grid.L), ncol= reps)
  vTtile <- matrix(NA, nrow= length(theta.grid), ncol= reps)
  
  
  
  
  omega1 <- omega2 <- 2*pi*Tp/spatPat$xDiff
  grid <- as.matrix(expand.grid(omega1, omega2))
  zero <- which(grid[,1]==0 & grid[,2]==0)
  angles <- atan(expand.grid(Tp,Tp)[,2]/expand.grid(Tp,Tp)[,1])
  angles[zero] <- NA
  
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
    
    prefG2 <- GLoc(spatPat=recPat, alpha=prefDir, epsilon=Gepsilon, r=r.grid.L, distance=prefG$distance, angles=prefG$angles)
    secG <- GLoc(spatPat=recPat, alpha=secDir, epsilon=Gepsilon, r=r.grid.L,  distance=prefG$distance, angles=prefG$angles)
    vGtile.L[,m] <- prefG2 - secG
    
    
    #K
    
    prefK <- cylK(spatPat=recPat, alpha=prefDir, aspect=Kaspect, r=r.grid.W)
    secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid.W)
    vKtile.W[,m] <- prefK-secK
    
    prefK <- cylK(spatPat=recPat, alpha=prefDir, aspect=Kaspect, r=r.grid.L)
    secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid.L)
    vKtile.L[,m] <- prefK-secK
    
    #Theta spectrum
    vTtile[,m] <- thetaSpectrum(prdgrm = periodogram(spatPat=recPat, p=Tp), 
                                thetas = theta.grid, epsilon = Tepsilon)
  }
  
  return(list(vGtile.W=vGtile.W, vGtile.L=vGtile.L, vKtile.W=vKtile.W, vKtile.L=vKtile.L, vTtile=vTtile))
}
