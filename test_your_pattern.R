# This scripts allows users to run a testing procedure for their pattern using tiling and stochastic reconstruction. 
# The script is written for square observation windows. Ajustments must be made for other observation window shapes.

# Input coordinates of points in yourr pattern:
x <- x_coordinates_of_your_points 
y <- y_coordinates_of_your_points 

# Input lower and upper bounds of the observation window in each dimension
xWin <- c(lower_x, upper_x)
yWin <- c(lower_y, upper_y)

# Select a test statistic
DSS <- "statistic" # input either "G", "K", or "T" for G-function, K-function or Theta spectrum

# if you selected G
Gepsilon <- NA # the epsilon parameter of the G-function, we recommend setting it to pi/4

# if you selected K
Kaspect <- NA # aspect ration of the rectangle

# if you selected T
Tp <- NA # a sequence of integeres, note that longer chains will increase the computation time, we recommend -15:15
Tepsilon <- NA # bandiwdth h, we recommend 7.5*pi/180

# for G and K, specify tested angles and a sequence of distances at wich the DSSs are esitmated

prefDir <- angle1
secDir <- angle2 # typically prefDir + pi/2
r.grid <- seq(0, rmax, length=51)[-1] # equidistant sequence of distances
  
# for T, specify a sequence of angles at which the DSS is esimated
theta.grid <- seq(0, pi, length=51)[-1] # only change the length argument to get finer or sparser sequence

# Select a replication method and number of replicates
repMethod <- "method_name" # set to either "tiling" or "SR" (for stochastic reconstruciton with spherical contact)
M <- num_repl # number of replicates; remember that stochastic reconstruction is much more computationally expensive than tiling

# if you selected tiling
sqrtN_tile <- 4 # a sqraure root of the number of tiles used in tiling

# if you selected stochastic reconstruction
numIter <-  5e3 # number of iterations for stochastic reconstruction algorithm, we recommend at least 5,000 iterations per 100 points.

# Select the number of cores you want to rune the parallelised code on
# We particularly suggest running parallelised code for big patterns or when stochastic reconstruction is used
cores <- 1 # Set to a number lower than the number of cores in your computer! Setting to 1 results in a non-parallelised code.

# Set a seed for reproducible results
seed <- NA

# Loading packages and functions
{
  library(spatstat)
  library(doParallel)
  
  angleMatrix <- function(points){
    angles <- matrix(data=NA, nrow=nrow(points), ncol=nrow(points))
    for(i in 1:nrow(points)){
      x <- points[i,1]-points[,1]
      y <- points[i,2]-points[,2]
      angles[i,] <- ifelse(x!=0, atan(y/x), 0)
      angles[i,angles[i,]<0] <- angles[i,angles[i,]<0]+pi
      angles[i,i] <- 0
    }
    return(angles)
  }
  
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
  
  rHomPois <- function(intensity=NA, xdim, ydim, numPoints=NA){
    window <-simWind(xdim, ydim)
    winSize <- window[[1]]
    if(is.na(numPoints)) numPoints <- max(1, rpois(1, intensity*winSize))
    x <- runif(numPoints, min=xdim[1], max=xdim[2])
    y <- runif(numPoints, min=ydim[1], max=ydim[2])
    return(list(pattern = cbind(x, y), winBor=window[[2]],
                lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff))
  }
  
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
  
  stchRecH <- function(spatPat, rmax=NULL, rlength=NULL, iter=1e4){
    
    xrange <- c(spatPat[[2]][1,1], spatPat[[2]][3,1])
    yrange <- c(spatPat[[2]][2,2], spatPat[[2]][3,2])
    
    # SET UP INTEGRATION GRID FOR R
    if(!is.null(rmax) & !is.null(rlength)){
      interval <- rmax/rlength
      rs <- seq(0, rmax, by=interval)
      
      originH <- Hest(X=as.ppp(spatPat[[1]], W=owin(xrange=xrange, yrange=yrange)), r = rs, correction = "raw", conditional=FALSE)$raw
    } else {
      # set automatically by spatstat so that H(rmax) = 1
      primH <- Hest(X=as.ppp(spatPat[[1]], W=owin(xrange=xrange, yrange=yrange)), correction = "raw", conditional=FALSE)
      originH <- primH$raw
      rs <- primH$r
      interval <- mean(diff(rs))
    }
    N <- nrow(spatPat[[1]])
    
    
    
    # INITIALISATION - BINOMIAL PROCESS with the same number of points as the original pattern
    curPat <- rHomPois(xdim=xrange, ydim=yrange, numPoints = N)
    curH <- Hest(X=as.ppp(curPat[[1]], W=owin(xrange=xrange, yrange=yrange)), r = rs, correction = "raw", conditional=FALSE)$raw
    curD <- sum(interval*(curH-originH)^2)
    for(i in 1:iter){
      move <- sample(1:N, 1) # sampled index of a point to be moved
      propPat <- curPat
      propPat[[1]][move,] <- c(runif(1, xrange[1], xrange[2]), runif(1, yrange[1], yrange[2])) # random new location
      propH <- Hest(X=as.ppp(propPat[[1]], W=owin(xrange=xrange, yrange=yrange)), r = rs, correction = "raw", conditional=FALSE)$raw
      propD <- sum(interval*(propH-originH)^2) # deviation for the H-function
      if(curD >= propD){
        accept <- 1 # accept if deviance improved
      } else {
        TS <- (0.1- 0.1*(i/(iter+1))^(1/10)) #low-temperature scaling temperatures
        Paccept <- exp((propD-curD)/TS)/100 # probability of acceptance if deviance isn't improved by the proposal
        accept <- runif(1, 0, 1)<Paccept
      }
      if(accept == 1){
        curPat <- propPat
        curH <- propH
        curD <- propD
      }
    }
    return(curPat)
  }
  
  testAniso <- function(OGv, Vs, DSS){
    if(!(DSS%in%c("G","K","T","P"))) stop("DSS must be either 'G', 'K', 'T', or 'P'.")
    reps <- ncol(Vs)
    mHat <- apply(Vs, 1, mean, na.rm=TRUE)
    vector <- OGv-mHat
    TstatRep <- rep(NA, reps)
    if(!(DSS=="G")){
      CHatDiag <- diag(1/apply(Vs, 1, var, na.rm=TRUE))
      CHatDiag[CHatDiag==Inf] <- 0 # if the variance is 0, a contribution at the corresponding distance 
      # will not be include to not overwhelm the test statistic
      Tstat <- t(vector)%*%CHatDiag%*%vector
      for(i in 1:reps){
        vecRep <- Vs[,i]-mHat
        TstatRep[i] <- t(vecRep)%*%CHatDiag%*%vecRep
      }
    } else {
      Tstat <- t(vector)%*%vector
      for(i in 1:reps){
        vecRep <- Vs[,i]-mHat
        TstatRep[i] <- t(vecRep)%*%vecRep
      }
    }
    p.val <- (1+sum(as.numeric(Tstat) < TstatRep, na.rm=TRUE))/(1+sum(!is.na(TstatRep)))
    return(list(p.val=p.val, Tstat=Tstat, TstatRep=TstatRep))
  }
  
  thetaSpectrum <- function(prdgrm, thetas, epsilon){ 
    spectrum <- rep(NA, length(thetas))
    for(i in 1:length(thetas)){
      consider <- (abs((prdgrm$angles-thetas[i])%%pi) < epsilon)
      spectrum[i] <- sum(prdgrm$sdf[consider],
                         na.rm=TRUE)/ sum(rowSums(consider, na.rm=TRUE))
    }
    return(spectrum)                 
  }
}

{
  if(diff(xWin)!=diff(xWin)) stop("the observation window is not a square!")
  # Setting up the point pattern object; re-centred at (0,0)
  xdim <- xWin - mean(xWin)
  ydim <- yWin - mean(xWin)
  window <- simWind(xdim, ydim)
  spatPat <- list(pattern = cbind(x-mean(xWin), y-mean(xWin)), winBor=window[[2]],
                  lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff)
  
  
  
  # sampling seeds for each of the replicates
  set.seed(seed)
  seeds <- sample(1:1e6, M, replace = FALSE)
  
  if(cores > 1){
    registerDoParallel(cores=cores)
    if(repMethod == "SR"){
      # PARALLELISED STOCHASTIC RECONSTRUCTION
      DDSvs <- foreach(l=1:M, .combine=cbind, .export = c("simWind", "angleMatrix", "Hest", "seeds",
                                                          "cylK", "GLoc", "periodogram", "thetaSpectrum", "stchRecH",
                                                          "owin", "xdim", "ydim", "DSS", "prefDir", "secDir",
                                                          "M", "Gepsilon", "Kaspect", "Tepsilon", "Tp", "theta.grid",
                                                          "as.ppp", "r.grid", "rHomPois", "numIter"),
                       .multicombine = TRUE, .inorder = TRUE, .maxcombine = M) %dopar% {
                         set.seed(seeds[l])
                         recPat <- stchRecH(spatPat, iter=numIter)
                         if(DSS == "G") {
                           prefG <- GLoc(spatPat=recPat, alpha=prefDir, epsilon=Gepsilon, r=r.grid, retMatrices = TRUE)
                           secG <- GLoc(spatPat=recPat, alpha=secDir, epsilon=Gepsilon, r=r.grid,
                                        distance=prefG$distance, angles=prefG$angles)
                           prefG$Gloc - secG
                         } else  if(DSS == "K") {
                           prefK <- cylK(spatPat=recPat, alpha=prefDir, aspect=Kaspect, r=r.grid)
                           secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid)
                           prefK-secK
                         } else {
                           thetaSpectrum(prdgrm = periodogram(spatPat=recPat, p=Tp),
                                         thetas = theta.grid, epsilon = Tepsilon)
                         }
                         
                       }
      
      
    } else {
      # PARALLELISED TILING
      # setting up tiles and subregions
      {
        n <- sqrtN_tile^2
        xDiff <- spatPat$xDiff
        yDiff <- spatPat$yDiff
        pattern <- t(spatPat$pattern)
        if(yDiff!=xDiff) stop("square pls")
        r <- xDiff*sqrt(2)/sqrtN_tile/2
        x.centres <- seq(spatPat$winBor[1,1]+r, spatPat$winBor[3,1]-r, length=sqrtN_tile)
        y.centres <- seq(spatPat$winBor[1,2]+r, spatPat$winBor[3,2]-r, length=sqrtN_tile)
        grid.centres <- expand.grid(x.centres, y.centres)
        if(nrow(grid.centres)!=n) stop("num. of points in the centres (reconstrutcted) grid not correct")
        x.tiles <- seq(spatPat$winBor[1,1]+xDiff/sqrtN_tile/2, spatPat$winBor[3,1]-xDiff/sqrtN_tile/2, length=sqrtN_tile)
        y.tiles <- seq(spatPat$winBor[1,2]+yDiff/sqrtN_tile/2, spatPat$winBor[3,2]-yDiff/sqrtN_tile/2, length=sqrtN_tile)
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
      
      DDSvs <- foreach(l=1:M, .combine=cbind, .export = c("simWind", "angleMatrix", "seeds", "sqrtN_tile", "n", "transTiles",
                                                          "cylK", "GLoc", "periodogram", "thetaSpectrum", "xDiff", "yDiff",
                                                          "DSS", "prefDir", "secDir", "pattern", "grid.tiles",
                                                          "M", "Gepsilon", "Kaspect", "Tepsilon", "Tp", "theta.grid",
                                                          "as.ppp", "r.grid"),
                       .multicombine = TRUE, .inorder = TRUE, .maxcombine = M) %dopar% {
                         set.seed(seeds[l])
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
                           rot.sq <- as.matrix(rotated[,(rotated[1,] <= xDiff/sqrtN_tile/2 & rotated[1,] >= -xDiff/sqrtN_tile/2 & 
                                                           rotated[2,] <= xDiff/sqrtN_tile/2 & rotated[2,] >= -xDiff/sqrtN_tile/2)])
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
                         if(DSS == "G") {
                           prefG <- GLoc(spatPat=recPat, alpha=prefDir, epsilon=Gepsilon, r=r.grid, retMatrices = TRUE)
                           secG <- GLoc(spatPat=recPat, alpha=secDir, epsilon=Gepsilon, r=r.grid,
                                        distance=prefG$distance, angles=prefG$angles)
                           prefG$Gloc - secG
                         } else  if(DSS == "K") {
                           prefK <- cylK(spatPat=recPat, alpha=prefDir, aspect=Kaspect, r=r.grid)
                           secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid)
                           prefK-secK
                         } else {
                           thetaSpectrum(prdgrm = periodogram(spatPat=recPat, p=Tp),
                                         thetas = theta.grid, epsilon = Tepsilon)
                         }
                         
                       }
      
    }
  } else {
    if(repMethod == "SR"){
      # STOCHASTIC RECONSTRUCTION WITHOUT PARALLELISATION
      if(DSS == "T") lt <- length(theta.grid) else lt <- length(r.grid)
      DDSvs <- matrix(NA, nrow = lt, ncol = M)
      for(i in 1:M){
        set.seed(seeds[i])
        recPat <- stchRecH(spatPat, iter=numIter)
        if(DSS == "G") {
          prefG <- GLoc(spatPat=recPat, alpha=prefDir, epsilon=Gepsilon, r=r.grid, retMatrices = TRUE)
          secG <- GLoc(spatPat=recPat, alpha=secDir, epsilon=Gepsilon, r=r.grid,
                       distance=prefG$distance, angles=prefG$angles)
          DDSvs[,i] <- prefG$Gloc - secG
        } else if(DSS == "K") {
          prefK <- cylK(spatPat=recPat, alpha=prefDir, aspect=Kaspect, r=r.grid)
          secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid)
          DDSvs[,i] <- prefK-secK
        } else {
          DDSvs[,i] <- thetaSpectrum(prdgrm = periodogram(spatPat=recPat, p=Tp),
                                     thetas = theta.grid, epsilon = Tepsilon)
        }
      }
    } else {
      # TILING WITHOUT PARALLELISATION
      # setting up tiles and subregions
      {
        n <- sqrtN_tile^2
        xDiff <- spatPat$xDiff
        yDiff <- spatPat$yDiff
        pattern <- t(spatPat$pattern)
        if(yDiff!=xDiff) stop("square pls")
        r <- xDiff*sqrt(2)/sqrtN_tile/2
        x.centres <- seq(spatPat$winBor[1,1]+r, spatPat$winBor[3,1]-r, length=sqrtN_tile)
        y.centres <- seq(spatPat$winBor[1,2]+r, spatPat$winBor[3,2]-r, length=sqrtN_tile)
        grid.centres <- expand.grid(x.centres, y.centres)
        if(nrow(grid.centres)!=n) stop("num. of points in the centres (reconstrutcted) grid not correct")
        x.tiles <- seq(spatPat$winBor[1,1]+xDiff/sqrtN_tile/2, spatPat$winBor[3,1]-xDiff/sqrtN_tile/2, length=sqrtN_tile)
        y.tiles <- seq(spatPat$winBor[1,2]+yDiff/sqrtN_tile/2, spatPat$winBor[3,2]-yDiff/sqrtN_tile/2, length=sqrtN_tile)
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
      
      if(DSS == "T") lt <- length(theta.grid) else lt <- length(r.grid)
      DDSvs <- matrix(NA, nrow = lt, ncol = M)
      
      # obtaining replicates
      for(m in 1:M){
        set.seed(seeds[m])
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
          rot.sq <- as.matrix(rotated[,(rotated[1,] <= xDiff/sqrtN_tile/2 & rotated[1,] >= -xDiff/sqrtN_tile/2 & 
                                          rotated[2,] <= xDiff/sqrtN_tile/2 & rotated[2,] >= -xDiff/sqrtN_tile/2)])
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
        
        if(DSS == "G") {
          prefG <- GLoc(spatPat=recPat, alpha=prefDir, epsilon=Gepsilon, r=r.grid, retMatrices = TRUE)
          secG <- GLoc(spatPat=recPat, alpha=secDir, epsilon=Gepsilon, r=r.grid,
                       distance=prefG$distance, angles=prefG$angles)
          DDSvs[,m] <- prefG$Gloc - secG
        } else if(DSS == "K") {
          prefK <- cylK(spatPat=recPat, alpha=prefDir, aspect=Kaspect, r=r.grid)
          secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid)
          DDSvs[,m] <- prefK-secK
        } else {
          DDSvs[,m] <- thetaSpectrum(prdgrm = periodogram(spatPat=recPat, p=Tp),
                                     thetas = theta.grid, epsilon = Tepsilon)
        }
        
        
      }
      
      
    }
  }
  
  
  # Computing the selected DSS for the observed pattern
  if(DSS == "G") {
    prefG <- GLoc(spatPat=spatPat, alpha=prefDir, epsilon=Gepsilon, r=r.grid, retMatrices = TRUE)
    secG <- GLoc(spatPat=spatPat, alpha=secDir, epsilon=Gepsilon, r=r.grid,
                 distance=prefG$distance, angles=prefG$angles)
    DSS_obs <- prefG$Gloc - secG
  } else if(DSS == "K") {
    prefK <- cylK(spatPat=spatPat, alpha=prefDir, aspect=Kaspect, r=r.grid)
    secK <- cylK(spatPat=spatPat, alpha=secDir, aspect=Kaspect, r=r.grid)
    DSS_obs <- prefK-secK
  } else {
    DSS_obs <- thetaSpectrum(prdgrm = periodogram(spatPat=spatPat, p=Tp),
                             thetas = theta.grid, epsilon = Tepsilon)
  }
  
  # Obtaining the test statistic
  results <- testAniso(OGv=DSS_obs, Vs=DDSvs, DSS=DSS)
  print(paste0("Your test's p-value is ", round(results$p.val, 4),
               ". Print object 'results' to view the test statistic and its replicates."))
}
