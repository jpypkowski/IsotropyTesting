library(doParallel)

# SIMULATION WINDOW GEOMETRIC - AUXILIARY FUNCTIONS
simWind <- function(xdim, ydim){
  xDiff <- xdim[2]-xdim[1]
  yDiff <- ydim[2]-ydim[1]
  winSize <- xDiff*yDiff
  if(xDiff < yDiff){
    xlim <- c(mean(xdim)-yDiff/2, mean(xdim)+yDiff/2)
    ylim <- ydim
  } else if(yDiff < xDiff){
    ylim <- c(mean(ydim)-xDiff/2, mean(ydim)+xDiff/2)
    xlim <- xdim
  } else {
    xlim <- xdim
    ylim <- ydim
  }
  return(list(winSize, winBor=cbind(c(xdim[1],xdim[2],xdim[2],xdim[1]),
                                    c(ydim[1],ydim[1],ydim[2],ydim[2])),
              lim=cbind(xlim, ylim), xDiff=xDiff, yDiff=yDiff))
}


data <- read.table("doi_10.5063_AA_connolly.206.1-DATA.data", header = TRUE)
uq <- unique(data[,1])
for(i in uq){
  ixs <- which(data[,1]==i)
  if(length(ixs) > 1) data <- data[-ixs[-1],]
}
window <- simWind(xdim=c(-50,50), ydim=c(-50,50))
Ambrosia <-list(pattern=cbind(as.numeric(data[,2])-50, as.numeric(data[,3])-50), winBor=window[[2]],
                lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff)


cylK <- function(spatPat, alpha, aspect=NA, r=seq(0,1, length=100),                 
                 contribution=FALSE){
  R <- matrix(c(cos(-alpha), -sin(-alpha), sin(-alpha), cos(-alpha)),
              nrow=2, byrow=TRUE)
  rotPattern <- t(R%*%t(spatPat[[1]]))
  
  lamHat <- nrow(spatPat[[1]])/(spatPat$xDiff*spatPat$yDiff)
  
  
  xDist <- spatPat$xDiff-as.matrix(dist(cbind(spatPat[[1]][,1],0), "euclidean"))
  yDist <- spatPat$yDiff-as.matrix(dist(cbind(spatPat[[1]][,2],0), "euclidean"))
  xRotDist <- as.matrix(dist(cbind(rotPattern[,1],0), "euclidean"))
  yRotDist <- as.matrix(dist(cbind(rotPattern[,2],0), "euclidean"))
  diag(xRotDist) <- NA
  diag(xRotDist) <- NA
  windowSize <- xDist*yDist
  contributions <- matrix(NA, ncol=length(r), nrow= nrow(spatPat[[1]]))
  #numerator <- rep(NA, length(r))
  for(i in 1:length(r)){
    contributions[,i] <- apply((xRotDist<=r[i] & yRotDist <=(aspect*r[i]))/windowSize, 1, sum, na.rm=TRUE)/lamHat^2
  }
  K <- apply(contributions, 2, sum, na.rm=TRUE)
  if(contribution){
    return(list(K=K, KContrib = contributions))
  } else {
    return(K=K)
  }
  
}




spatPat = Ambrosia 
sqrt.n = 5
reps = 1000 
r.grid.W = seq(0, 25, length=101)[-1]
prefDir=0
Kaspect=0.15
  
  
  
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
  
  
  
  
  
set.seed(3611)
seeds <- sample(1e6, 1e3, replace = FALSE)

registerDoParallel(cores=8)
vKtile <- foreach(l=1:reps, .combine=cbind, .export = c("recPat", "prefDir", "secDir", "spatPat", "yDiff", "seeds",
                                           "Kaspect", "r.grid.W", "cylK", "n", "transTiles", "xDiff", "grid.tiles"), 
        .multicombine = TRUE, .maxcombine = 8, .inorder=TRUE) %dopar% {
    # sampling angles and tiles
    set.seed(seeds[l])
    
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
    
   
    
            prefK <- cylK(spatPat=recPat, alpha=prefDir, aspect=Kaspect, r=r.grid.W)
            secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid.W)
            

    print(l)

    prefK-secK
  }
  
write.csv(vKtile, "vKtile_k5_025.csv")
