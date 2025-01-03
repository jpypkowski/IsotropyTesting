# PLOTTING PAPER
library(latex2exp)
library(pROC)
library(mvtnorm)
library(spatstat)
library(circular)
{
  
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
  
  # HOMOGENOUS POISSON POINT PROCES SIMULATION
  rHomPois <- function(intensity=NA, xdim, ydim, numPoints=NA){
    window <-simWind(xdim, ydim)
    winSize <- window[[1]]
    if(is.na(numPoints)) numPoints <- max(1, rpois(1, intensity*winSize))
    x <- runif(numPoints, min=xdim[1], max=xdim[2])
    y <- runif(numPoints, min=ydim[1], max=ydim[2])
    return(list(pattern = cbind(x, y), winBor=window[[2]],
                lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff))
  }
  
  
  # Adjusting spatstat's simulation of LGCP to our data format
  
  rLGCPtransformed <- function(model, mu, xdim, ydim, ...){
    pat <- rLGCP(model=model, mu=mu, win=owin(xrange=xdim, yrange=ydim), ...)
    window <-simWind(xdim, ydim)
    return(list(pattern = cbind(pat$x, pat$y), winBor=window[[2]],
                lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff))
  }
  
  # POISSON LINE CLUSTER PROCESS SIMULATION
  
  rDirLine <- function(lineIntensity, alongIntensity, theta, xdim, ydim, sigma, kappa, margin=0.05){
    window <- simWind(xdim, ydim)
    
    xDiff <- window$xDiff+2*margin
    yDiff <- window$xDiff+2*margin
    r <- sqrt(xDiff^2+yDiff^2)/2
    noLines <- max(1, rpois(1, 2*r*lineIntensity))
    intersections <- runif(noLines, -r, r)
    thetas <- as.numeric(rvonmises(n=noLines, mu=theta, kappa=kappa))
    
    
    # plot(1,1,type="n", xlim=c(-r,r), ylim=c(-r,r))
    # for(i in 1:length(intersections)){
    #   points(-intersections[i]*sin(thetas[i]), intersections[i]*cos(thetas[i]), col="blue", pch=16)
    # }
    
    
    observedLine <- matrix(NA, nrow=1, ncol=2)
    for(i in 1:noLines){
      lineNumPoints <- rpois(1, alongIntensity*2*r)
      linePointsAcross <- rnorm(lineNumPoints, 0, sigma) #x-coord
      linePointsAlong <- runif(lineNumPoints, -r, r)                    #y-coord
      linePoints <- cbind(linePointsAlong, linePointsAcross+intersections[i])
      
      # points(linePointsAlong, linePointsAcross)
      # points(linePointsAlong, linePointsAcross+intersections[i], col=2)
      R <- matrix(c(cos(thetas[i]),-sin(thetas[i]), sin(thetas[i]), cos(thetas[i])), 
                  byrow=TRUE, ncol=2)
      rotated <- t(R%*%t(linePoints))
      # points(rotated[,1], rotated[,2], col=3)
      observedLine <- rbind(observedLine, cbind(rotated[,1], rotated[,2]))
      # points(observedLine[,1], observedLine[,2])
    }
    if(nrow(observedLine)>2){
      observedLine <- observedLine[-1,]
    } else {
      observedLine <- matrix(observedLine[-1,], ncol=2)
    }
    
    inWindow <- (observedLine[,1] >= xdim[1] & observedLine[,1] <= xdim[2] &
                   observedLine[,2] >= ydim[1] & observedLine[,2] <= ydim[2])
    observedPoints <- observedLine[inWindow,]
    
    
    if(sum(inWindow==1)){
      return(list(pattern=matrix(observedPoints, ncol=2), winBor=window$winBor, 
                  lim=window$lim,  xDiff=window$xDiff, yDiff=window$yDiff, intersections = intersections, theta =theta))
    } else {
      return(list(pattern=observedPoints, winBor=window$winBor, 
                  lim=window$lim,  xDiff=window$xDiff, yDiff=window$yDiff, intersections = intersections, theta =theta))
    }
    
  }
  
  # AUXILIARY FUNCTION THAT CONVERTS ANISOTROPY PARAMETER a TO kappa
  
  KofA <- function(a) 5*(1-exp(1-1/a))
  
  # Function to apply a transformation leading to geometric anisotropy (used for LGCP)
  
  rAnisotrop <- function(theta, xScalingFactor, process, xdim, ydim, ...){
    R <- matrix(c(cos(theta),-sin(theta), sin(theta), cos(theta)), 
                byrow=TRUE, ncol=2)
    C <- diag(c(xScalingFactor, 1/xScalingFactor))
    Trans <- R%*%C
    window <- simWind(xdim, ydim)
    extrMat <- solve(Trans) %*% t(window$winBor)
    xdim2 <- c(min(extrMat[1,]), max(extrMat[1,]))
    ydim2 <- c(min(extrMat[2,]), max(extrMat[2,]))
    isoPattern <- process(xdim=xdim2, ydim=ydim2, ...)[[1]]
    anisoPattern <- t(Trans%*%t(isoPattern))
    inWindow <- (anisoPattern[,1] >= xdim[1] & anisoPattern[,1] <= xdim[2] &
                   anisoPattern[,2] >= ydim[1] & anisoPattern[,2] <= ydim[2])
    observedAnisoPattern <- anisoPattern[inWindow,]
    return(list(pattern=observedAnisoPattern, winBor=window$winBor, 
                lim=window$lim, xDiff=window$xDiff, yDiff=window$yDiff))
  }
  
  
  # PLOTTING A SPATIAL PATTERN ON A RECTANGULAR OBSERVATION WINDOW
  plotSpatPat <- function(spatPat, axes=FALSE, 
                          xlab="", ylab="", ...){
    plot(spatPat$pattern, frame.plot = FALSE, axes=axes,
         xlim=spatPat$lim[,1], ylim=spatPat$lim[,2], 
         xlab=xlab, ylab=ylab, ...)
    polygon(x=spatPat$winBor[,1], y=spatPat$winBor[,2], 
            lty=1, lwd=1, border="darkgrey")
  }
  
  
  # Anisotropic 6-12 Lennard-Jones energy function
  
  fAnisoLennardJones <- function(spatPat, theta, alpha, epsilon1, epsilon2, sigma1, sigma2, log=FALSE){
    n <- nrow(spatPat$pattern)
    if(n>1){
      distances <- as.matrix(dist(spatPat$pattern))
      diag(distances) <- NA
      angles <- angleMatrix(spatPat$pattern)
      diag(angles) <- NA
      
      if(theta < pi/4){
        angles <- (angles < (theta+pi/4) | angles > ((theta-pi/4)%%pi))
      } else if (theta > 3*pi/4){
        angles <- (angles < ((theta+pi/4)%%pi) | angles > (theta-pi/4))
      } else {
        angles <- (angles < (theta+pi/4) & angles > (theta-pi/4))
      }
      
      PHI <- 4*epsilon1* ( (sigma1/distances)^12 - (sigma1/distances)^6 ) * angles + 4*epsilon2* ( (sigma2/distances)^12 - (sigma2/distances)^6 ) * (!angles)
      if(log) {
        return(-(n*alpha+sum(apply(PHI, 2, sum, na.rm=TRUE))/2))
      } else{
        return(exp(-(n*alpha+sum(apply(PHI, 2, sum, na.rm=TRUE))/2)))
      }
    } else if(n==1) {
      if(log) {
        return(-alpha)
      } else{
        return(exp(-alpha))
      }
    } else {
      if(log) {
        return(0)
      } else{
        return(1)
      }
    }
  }
  
  # Single replicate of tiling
  
  singleBoot <- function(spatPat, sqrt.n){
    n <- sqrt.n^2
    xDiff <- spatPat$xDiff
    yDiff <- spatPat$yDiff
    pattern <- rbind(t(spatPat$pattern), 1:nrow(spatPat$pattern))
    if(yDiff!=xDiff) stop("square pls")
    r <- xDiff*sqrt(2)/sqrt.n/2
    x.centres <- seq(spatPat$winBor[1,1]+r, spatPat$winBor[3,1]-r, length=sqrt.n)
    y.centres <- seq(spatPat$winBor[1,2]+r, spatPat$winBor[3,2]-r, length=sqrt.n)
    grid.centres <- expand.grid(x.centres, y.centres)
    if(nrow(grid.centres)!=n) stop("num. of points in the centres (reconstrutcted) grid not correct")
    
    angles <- runif(n, min=0, max=2*pi)
    
    x.tiles <- seq(spatPat$winBor[1,1]+xDiff/sqrt.n/2, spatPat$winBor[3,1]-xDiff/sqrt.n/2, length=sqrt.n)
    y.tiles <- seq(spatPat$winBor[1,2]+yDiff/sqrt.n/2, spatPat$winBor[3,2]-yDiff/sqrt.n/2, length=sqrt.n)
    grid.tiles <- expand.grid(x.tiles, y.tiles)
    if(nrow(grid.tiles)!=n) stop("num. of points in the tiles (original) grid not correct")
    
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
    
    
    
    which.tiles <- sample(1:n, n, replace=TRUE)
    
    new.pattern <- vector(mode="list", length=n)
    sampled.points <- vector(mode="list", length=n)
    lengths <- rep(NA, n)
    for(i in 1:n){
      tile <- which.tiles[i]
      aux.pattern <- transTiles[[tile]]
      
      angle <- angles[i]
      R <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),ncol=2, byrow=TRUE)
      rotated <- rbind(R%*%(as.matrix(aux.pattern[-3,])), aux.pattern[3,])
      rot.sq <- as.matrix(rotated[,(rotated[1,] <= xDiff/sqrt.n/2 & rotated[1,] >= -xDiff/sqrt.n/2 & 
                                      rotated[2,] <= xDiff/sqrt.n/2 & rotated[2,] >= -xDiff/sqrt.n/2)])
      new.pattern[[i]] <- cbind(rot.sq[1,]+grid.tiles[i,1], rot.sq[2,]+grid.tiles[i,2])
      sampled.points[[i]] <- rot.sq[3,] 
      lengths[i] <- nrow(new.pattern[[i]])
    }
    reconstr <- matrix(NA, nrow=sum(lengths), ncol=2)
    MPM.points <- rep(NA, sum(lengths))
    for(i in 1:n){
      if(i==1 && lengths[i]>0) {
        reconstr[1:lengths[1],] <- new.pattern[[1]]
        MPM.points[1:lengths[1]] <- sampled.points[[1]]
      } else {
        if(lengths[i]>0) {
          reconstr[(1+sum(lengths[1:(i-1)])):sum(lengths[1:i]),] <- new.pattern[[i]]
          MPM.points[(1+sum(lengths[1:(i-1)])):sum(lengths[1:i])] <- sampled.points[[i]]
        }
        
        
      }
      
      
    }
    return(list(pattern = reconstr, winBor=spatPat$winBor,
                lim=spatPat$lim, xDiff=xDiff, yDiff=yDiff, MPM.points=MPM.points))
  }
  
  # STOCHASTIC RECONSTRUCTION ALGORITHM WITH H-FUNCTION
  
  stchRecH <- function(spatPat, rmax=NULL, rlength=NULL, iter=1e4){
    
    xrange <- c(spatPat[[2]][1,1], spatPat[[2]][3,1])
    yrange <- c(spatPat[[2]][2,2], spatPat[[2]][3,2])
    
    # SET UP INTEGRATION GRID FOR R
    if(!is.null(rmax) & !is.null(rlength)){
      interval <- rmax/rlength
      rs <- seq(0, rmax, by=interval)
      
      originH <- Hest(X=as.ppp(spatPat[[1]], W=owin(xrange=xrange, yrange=yrange)), r = rs, correction = "raw", conditional=FALSE)$raw
    } else {
      primH <- Hest(X=as.ppp(spatPat[[1]], W=owin(xrange=xrange, yrange=yrange)), correction = "raw", conditional=FALSE)
      originH <- primH$raw
      rs <- primH$r
      interval <- mean(diff(rs))
    }
    N <- nrow(spatPat[[1]])
    
    
    
    # INITIALISATION - BINOMIAL PROCESS, U, KMatrix
    curPat <- rHomPois(xdim=xrange, ydim=yrange, numPoints = N)
    curH <- Hest(X=as.ppp(curPat[[1]], W=owin(xrange=xrange, yrange=yrange)), r = rs, correction = "raw", conditional=FALSE)$raw
    curD <- sum(interval*(curH-originH)^2)
    for(i in 1:iter){
      move <- sample(1:N, 1)
      propPat <- curPat
      propPat[[1]][move,] <- c(runif(1, xrange[1], xrange[2]), runif(1, yrange[1], yrange[2]))
      propH <- Hest(X=as.ppp(propPat[[1]], W=owin(xrange=xrange, yrange=yrange)), r = rs, correction = "raw", conditional=FALSE)$raw
      propD <- sum(interval*(propH-originH)^2)
      if(curD >= propD){
        accept <- 1
      } else {
        TS <- (0.1- 0.1*(i/(iter+1))^(1/10))
        Paccept <- exp((propD-curD)/TS)/100
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
  
}
black <- function(alpha=255) rgb(0, 0, 0, max = 255, alpha = alpha)
teal <- function(alpha=255) rgb(0, 153, 136, max = 255, alpha = alpha)
orange <- function(alpha=255) rgb(238, 119, 51, max = 255, alpha = alpha)
magenta <- function(alpha=255) rgb(238, 51, 119, max = 255, alpha = alpha)
grey <- function(alpha=255) rgb(187, 187, 187, max = 255, alpha = alpha)
violet <- function(alpha=255) rgb(238, 130, 238, max = 255, alpha = alpha)




# PATTERN EXAMPLES

{

var=3
mu=log(400)-var/2

set.seed(645) #22
LGCP <- rAnisotrop(theta=pi/6, 0.6, process=rLGCPtransformed, 
                   xdim=c(0,1), ydim=c(0,1), model="exponential", mu=mu, var=var, scale=0.02)
plotSpatPat(LGCP)

patterns <- read.csv("IsotropyTesting\\GIBBS\\Gibbs06_big_2.csv", header=TRUE)[,-1]
ixs <- which(patterns[,1]=="pattern")
l=18 #18  #845
ix <- ixs[l:(l+1)]
window <- simWind(c(-0.5,0.5), c(-0.5,0.5))
np1 <- patterns[(ix[1]+1):(ix[2]-1),]
np1 <- cbind(as.numeric(np1[,1]), as.numeric(np1[,2]))
np <- np1[np1[,1]>(-0.5) & np1[,1]<(0.5) & np1[,2]>(-0.5) & np1[,2]<(0.5),]
GIBBS <- list(pattern = np, winBor=window[[2]],
              lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff)
plotSpatPat(GIBBS)


set.seed(786)   #78
PLCP <- rDirLine(lineIntensity=16, alongIntensity=25, theta=pi/6, xdim=c(-0.5,0.5), ydim=c(-0.5,0.5), sigma=0.015, kappa=KofA(0.6))


par(mfrow=c(1,3), mar=c(1,1,1,1))
plotSpatPat(LGCP, cex = 0.75)
plotSpatPat(GIBBS, cex = 0.75)
plotSpatPat(PLCP, cex = 0.75)

}


# TILING PLOT
{
  
  cols <- c(black, violet, orange, teal)
  
  var=3
  scale= 0.03
  mu=log(200)-var/2
  model="exponential"
  par(mfrow=c(1,2))
  set.seed(894) #222
  example <- rAnisotrop(theta=pi/6, 0.4, process=rLGCPtransformed, 
                        xdim=c(0,1), ydim=c(0,1), model=model, mu=mu, var=var, scale=scale)
  r <- 1/2*sqrt(2)/2
  
  
  x.centres <- c(r, r, 1-r, 1-r)
  y.centres <- c(r, 1-r, r, 1-r)
  
  
  for(i in 1:4){
    xs <- seq(x.centres[i]-r, x.centres[i]+r, length=100)
    vrtx <- sqrt(ifelse(r^2-(xs-x.centres[i])^2>0,r^2-(xs-x.centres[i])^2, 0 ))+y.centres[i]
    vrtx <- c(vrtx[100:1], -sqrt(ifelse(r^2-(xs-x.centres[i])^2>0,r^2-(xs-x.centres[i])^2,0)) +y.centres[i])
    polygon(c(xs[100:1], xs), vrtx)
  }
  
  sqr <-  matrix(c(-1/4, 1/4, 1/4, -1/4, -1/4, -1/4, 1/4, 1/4), ncol=2, byrow=TRUE)
  set.seed(66)
  angles <- runif(4, 0, 2*pi)
  set.seed(12)
  tiles <- sample(1:4, 4, replace=TRUE)
  
  x.tileCen <- c(0.25, 0.75, 0.25, 0.75)
  y.tileCen <- c(0.75, 0.75, 0.25, 0.25)
  
  basevrtx <- matrix(c(-0.25, 0.25, 0.25, -0.25, -0.25, -0.25, 0.25, 0.25), nrow=2, byrow=TRUE)
  
  
  
  
  par(mfrow=c(1,2), mar=c(1,1,1,1))
  plotSpatPat(example, cex=0.5)
  for(i in 1:4){
    xs <- seq(x.centres[i]-r, x.centres[i]+r, length=100)
    vrtx <- sqrt(ifelse(r^2-(xs-x.centres[i])^2>0,r^2-(xs-x.centres[i])^2, 0 ))+y.centres[i]
    vrtx <- c(vrtx[100:1], -sqrt(ifelse(r^2-(xs-x.centres[i])^2>0,r^2-(xs-x.centres[i])^2,0)) +y.centres[i])
    polygon(c(xs[100:1], xs), vrtx, lty=2, border="darkgrey")
    angle <- -angles[i]
    R <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),ncol=2, byrow=TRUE)
    sqvrtx <- R%*%basevrtx
    colour  <- cols[i][[1]]
    polygon(sqvrtx[1,]+ x.centres[tiles[i]], sqvrtx[2,]+ y.centres[tiles[i]],
            col=colour(alpha=50), border = colour(alpha=255))
  }
  plotSpatPat(example, type="n")
  for(i in 1:4){
    
    tile <- tiles[i]
    pat_trans <- cbind(example[[1]][,1]- x.centres[tile], example[[1]][,2]-y.centres[tile])
    distances <- sqrt(pat_trans[,1]^2+pat_trans[,2]^2)
    pat_trans_r <- pat_trans[distances <= r,]
    angle <- angles[i]
    R <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),ncol=2, byrow=TRUE)
    rotated <- t(R%*%t(pat_trans_r))
    rotSq <- rotated[(rotated[,1] <= 0.25 & rotated[,1] >= -0.25 & 
                        rotated[,2] <= 0.25 & rotated[,2] >= -0.25),]  
    rotated_trans <- cbind(rotSq[,1]+x.tileCen[i], rotSq[,2]+y.tileCen[i])
    colour  <- cols[i][[1]]
    points(rotated_trans, cex=0.5)#, col=colour(alpha=255), cex=0.5)
    
  }
  polygon(c(0, 0.5, 0.5, 0), c(0.5, 0.5, 1, 1),
          col=cols[1][[1]](alpha=50), border = cols[1][[1]](alpha=255))
  polygon(c(0.5, 1, 1, 0.5), c(0.5, 0.5, 1, 1),
          col=cols[2][[1]](alpha=50), border = cols[2][[1]](alpha=255))
  polygon(c(0, 0.5, 0.5, 0), c(0, 0, 0.5, 0.5),
          col=cols[3][[1]](alpha=50), border = cols[3][[1]](alpha=255))
  polygon(c(0.5, 1, 1, 0.5), c(0, 0, 0.5, 0.5),
          col=cols[4][[1]](alpha=50), border = cols[4][[1]](alpha=255))
  
}



{

singleBoot <- function(spatPat, sqrt.n){
  n <- sqrt.n^2
  xDiff <- spatPat$xDiff
  yDiff <- spatPat$yDiff
  pattern <- rbind(t(spatPat$pattern), 1:nrow(spatPat$pattern))
  if(yDiff!=xDiff) stop("square pls")
  r <- xDiff*sqrt(2)/sqrt.n/2
  x.centres <- seq(spatPat$winBor[1,1]+r, spatPat$winBor[3,1]-r, length=sqrt.n)
  y.centres <- seq(spatPat$winBor[1,2]+r, spatPat$winBor[3,2]-r, length=sqrt.n)
  grid.centres <- expand.grid(x.centres, y.centres)
  if(nrow(grid.centres)!=n) stop("num. of points in the centres (reconstrutcted) grid not correct")
  
  angles <- runif(n, min=0, max=2*pi)
  
  x.tiles <- seq(spatPat$winBor[1,1]+xDiff/sqrt.n/2, spatPat$winBor[3,1]-xDiff/sqrt.n/2, length=sqrt.n)
  y.tiles <- seq(spatPat$winBor[1,2]+yDiff/sqrt.n/2, spatPat$winBor[3,2]-yDiff/sqrt.n/2, length=sqrt.n)
  grid.tiles <- expand.grid(x.tiles, y.tiles)
  if(nrow(grid.tiles)!=n) stop("num. of points in the tiles (original) grid not correct")
  
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
  
  
  
  which.tiles <- sample(1:n, n, replace=TRUE)
  
  new.pattern <- vector(mode="list", length=n)
  sampled.points <- vector(mode="list", length=n)
  lengths <- rep(NA, n)
  for(i in 1:n){
    tile <- which.tiles[i]
    aux.pattern <- transTiles[[tile]]
    
    angle <- angles[i]
    R <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),ncol=2, byrow=TRUE)
    rotated <- rbind(R%*%(as.matrix(aux.pattern[-3,])), aux.pattern[3,])
    rot.sq <- as.matrix(rotated[,(rotated[1,] <= xDiff/sqrt.n/2 & rotated[1,] >= -xDiff/sqrt.n/2 & 
                                    rotated[2,] <= xDiff/sqrt.n/2 & rotated[2,] >= -xDiff/sqrt.n/2)])
    new.pattern[[i]] <- cbind(rot.sq[1,]+grid.tiles[i,1], rot.sq[2,]+grid.tiles[i,2])
    sampled.points[[i]] <- rot.sq[3,] 
    lengths[i] <- nrow(new.pattern[[i]])
  }
  reconstr <- matrix(NA, nrow=sum(lengths), ncol=2)
  MPM.points <- rep(NA, sum(lengths))
  for(i in 1:n){
    if(i==1 && lengths[i]>0) {
      reconstr[1:lengths[1],] <- new.pattern[[1]]
      MPM.points[1:lengths[1]] <- sampled.points[[1]]
    } else {
      if(lengths[i]>0) {
        reconstr[(1+sum(lengths[1:(i-1)])):sum(lengths[1:i]),] <- new.pattern[[i]]
        MPM.points[(1+sum(lengths[1:(i-1)])):sum(lengths[1:i])] <- sampled.points[[i]]
      }
      
      
    }
    
    
  }
  return(list(pattern = reconstr, winBor=spatPat$winBor,
              lim=spatPat$lim, xDiff=xDiff, yDiff=yDiff, MPM.points=MPM.points))
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
    primH <- Hest(X=as.ppp(spatPat[[1]], W=owin(xrange=xrange, yrange=yrange)), correction = "raw", conditional=FALSE)
    originH <- primH$raw
    rs <- primH$r
    interval <- mean(diff(rs))
  }
  N <- nrow(spatPat[[1]])
  
  
  
  # INITIALISATION - BINOMIAL PROCESS, U, KMatrix
  curPat <- rHomPois(xdim=xrange, ydim=yrange, numPoints = N)
  curH <- Hest(X=as.ppp(curPat[[1]], W=owin(xrange=xrange, yrange=yrange)), r = rs, correction = "raw", conditional=FALSE)$raw
  curD <- sum(interval*(curH-originH)^2)
  for(i in 1:iter){
    move <- sample(1:N, 1)
    propPat <- curPat
    propPat[[1]][move,] <- c(runif(1, xrange[1], xrange[2]), runif(1, yrange[1], yrange[2]))
    propH <- Hest(X=as.ppp(propPat[[1]], W=owin(xrange=xrange, yrange=yrange)), r = rs, correction = "raw", conditional=FALSE)$raw
    propD <- sum(interval*(propH-originH)^2)
    if(curD >= propD){
      accept <- 1
    } else {
      TS <- (0.1- 0.1*(i/(iter+1))^(1/10))
      Paccept <- exp((propD-curD)/TS)/100
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

infer.LGCP <- function(spatPat, covariance, startpar=c(var=1, scale=0.8), ...){
  estimates <- rep(NA, 3) # variance, scale, mu
  lambda <- nrow(spatPat$pattern)/(spatPat$xDiff*spatPat$yDiff)
  estimates[1:2] <- lgcp.estpcf(as.ppp(spatPat$pattern, 
                                       W=owin(c(spatPat$winBor[1,1], spatPat$winBor[2,1]), 
                                              c(spatPat$winBor[1,2], spatPat$winBor[3,2]))),
                                startpar=startpar, covmodel = list(model=covariance, ...), lambda=lambda)$par
  estimates[3] <- log(lambda)-estimates[1]/2
  return(estimates)
  
}


# Adjusting spatstat's simulation of LGCP to our data format

rLGCPtransformed <- function(model, mu, xdim, ydim, ...){
  pat <- rLGCP(model=model, mu=mu, win=owin(xrange=xdim, yrange=ydim), ...)
  window <-simWind(xdim, ydim)
  return(list(pattern = cbind(pat$x, pat$y), winBor=window[[2]],
              lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff))
}

infer.Thomas <- function(pattern, startpar=c(kappa=10,scale=0.1)){
  estimates <- rep(NA, 3)
  xrange <- c(pattern[[2]][1,1], pattern[[2]][3,1])
  yrange <- c(pattern[[2]][2,2], pattern[[2]][3,2])
  estimates[1:2] <- thomas.estK(as.ppp(pattern$pattern, owin(xrange=xrange, yrange=yrange)), startpar=startpar)$par
  estimates[2] <- sqrt(estimates[2])
  estimates[3] <- nrow(pattern$pattern)/estimates[1]/(diff(xrange)*diff(yrange))
  return(estimates)
}

rThomas <- function(parentIntensity, meanClusterSize, 
                    sigma, theta=0, xScalingFactor=1, xdim, ydim, addMargin=4){
  
  xdim2 <- xdim+c(-addMargin*sigma, addMargin*sigma)
  ydim2 <- ydim+c(-addMargin*sigma, addMargin*sigma)
  trueWindow <- simWind(xdim, ydim)
  parents <- rHomPois(parentIntensity, xdim2, ydim2)
  clusterSizes <- rpois(nrow(parents[[1]]), meanClusterSize)
  offsprings <- matrix(NA, nrow=sum(clusterSizes), ncol=2)
  if(xScalingFactor==1){
    for(i in 1:nrow(parents[[1]])){
      stoppedAt <- max(c(0,which(!is.na(offsprings[,1]))))
      cs <- clusterSizes[i]
      if(cs>0) offsprings[(stoppedAt+1):(stoppedAt+cs),]  <- rmvnorm(cs, 
                                                                     parents[[1]][i,], diag(sigma^2, nrow=2))
    }
  } else {
    R <- matrix(c(cos(theta),-sin(theta), sin(theta), cos(theta)), 
                byrow=TRUE, ncol=2)
    C <- diag(c(xScalingFactor, 1/xScalingFactor))
    Trans <- R%*%C
    for(i in 1:nrow(parents[[1]])){
      stoppedAt <- max(c(0,which(!is.na(offsprings[,1]))))
      cs <- clusterSizes[i]
      if(cs>0) offsprings[(stoppedAt+1):(stoppedAt+cs),]  <- rmvnorm(cs, 
                                                                     parents[[1]][i,], Trans%*%diag(sigma^2, nrow=2)%*%t(Trans))
    }
  }
  
  inWindow <- (offsprings[,1] >= xdim[1] & offsprings[,1] <= xdim[2] &
                 offsprings[,2] >= ydim[1] & offsprings[,2] <= ydim[2])
  observedOffsprings <- offsprings[inWindow,]
  return(list(pattern=observedOffsprings, winBor=trueWindow$winBor, 
              lim=trueWindow$lim, xDiff=trueWindow$xDiff, yDiff=trueWindow$yDiff))
}

infer.Strauss <- function(pattern){
  estimates <- rep(NA, 3)
  xrange <- c(pattern[[2]][1,1], pattern[[2]][3,1])
  yrange <- c(pattern[[2]][2,2], pattern[[2]][3,2])
  L <- Lest(as.ppp(pattern$pattern, owin(xrange=xrange, yrange=yrange)))
  ix <-  which(((cumsum(L$theo-L$trans))/L$r)[-1]==max(((cumsum(L$theo-L$trans))/L$r)[-1]))
  estimates[3] <- L$r[1+ix]
  estimates[1:2] <- exp(ppm(as.ppp(pattern$pattern, owin(xrange=xrange, yrange=yrange)) ~ 1, Strauss(estimates[3]))$coef)
  return(estimates)
}

#STRAUSS SIMULATION USING METROPOLIS-HASTINGS

rmhStrauss_rep <- function(beta, gamma, R, xdim, ydim){
  pat <- rmh(list(cif="strauss",par=list(beta=beta,gamma=gamma,r=R),
                  w=c(xdim, ydim)), nrep=1e6, verbose=FALSE)
  window <-simWind(xdim, ydim)
  return(list(pattern=cbind(pat$x, pat$y), winBor=window$winBor, 
              lim=window$lim,  xDiff=window$xDiff, yDiff=window$yDiff))
}

# BAYESIAN INFERENCE FOR PLCP + ALL AUXILIARY FUNCTIONS

int_sum <- function(points, lines_new, sigma_sq, full = TRUE){
  #lines_new - matrix for many lines and full == TRUE
  #            or vector for only one new line
  if(full){
    vector <- rep(NA, nrow(lines_new))
    for(j in 1:nrow(lines_new)){
      point_diff <- cbind(points[,1] - lines_new[j,1], points[,2])
      projection <- point_diff - t(apply(point_diff, 1, function(x) sum(x*lines_new[j, 2:3])*lines_new[j, 2:3]))
      proj_dist <- apply(projection, 1, function(x) sqrt(x[1]^2+x[2]^2))
      vector[j] <- sum(dnorm(proj_dist, mean = 0, sd = sqrt(sigma_sq)))
    }
  } else {
    point_diff <- cbind(points[,1] - lines_new[1], points[,2])
    projection <- point_diff - t(apply(point_diff, 1, function(x) sum(x*lines_new[2:3])*lines_new[2:3]))
    proj_dist <- apply(projection, 1, function(x) sqrt(x[1]^2+x[2]^2))
    vector <- sum(dnorm(proj_dist, mean = 0, sd = sqrt(sigma_sq))) # in this case vector is a scalar!
  }
  
  return(vector)
}

sum_log_matrix <- function(points, lines_new, sigma_sq, full = TRUE){
  if(full){
    func_vals <- matrix(NA, nrow=nrow(points), ncol=nrow(lines_new))
    for(j in 1:nrow(lines_new)){
      point_diff <- cbind(points[,1] - lines_new[j,1], points[,2])
      projection <- point_diff - t(apply(point_diff, 1, function(x) sum(x*lines_new[j, 2:3])*lines_new[j, 2:3]))
      proj_dist <- apply(projection, 1, function(x) sqrt(x[1]^2+x[2]^2))
      func_vals[,j] <- dnorm(proj_dist, mean = 0, sd = sqrt(sigma_sq))
    }
  } else {
    point_diff <- cbind(points[,1] - lines_new[1], points[,2])
    projection <- point_diff - t(apply(point_diff, 1, function(x) sum(x*lines_new[2:3])*lines_new[2:3]))
    proj_dist <- apply(projection, 1, function(x) sqrt(x[1]^2+x[2]^2))
    func_vals <- dnorm(proj_dist, mean = 0, sd = sqrt(sigma_sq)) # in this case func_vals is a vector, not a matrix!
  }
  
  return(func_vals)
}

Rbirth <- function(l_prop, points, lines, lambda_J, rho_L, alpha, sigma_sq, MC_points, dens_mat, point_area, log=FALSE){
  log_1st_bit <- log(rho_L * lambda_J * abs(l_prop[3]) / (nrow(lines) + 1))
  int_prop_l <- int_sum(MC_points, lines=matrix(l_prop, nrow=1), sigma_sq = sigma_sq, full = FALSE) # M_integral
  log_2nd_bit <- -alpha*int_prop_l*point_area
  
  dens_prep_l <- sum_log_matrix(points, l_prop, sigma_sq, full = FALSE)
  
  if(is.matrix(dens_mat)){
    lll <- log(1 + dens_prep_l/apply(dens_mat, MARGIN = 1, sum))
    lll[lll==Inf] <- NaN
    sum_log <- sum(lll, na.rm=TRUE)
  } else {
    lll <- log(1 + dens_prep_l/sum(dens_mat))
    lll[lll=Inf] <- NaN
    sum_log <- sum(log(1 + dens_prep_l/sum(dens_mat)))
  }
  
  if(log){
    return(list(logRb = log_1st_bit + log_2nd_bit + sum_log, int = int_prop_l, dens = dens_prep_l))
  } else {
    return(list(Rb = exp(log_1st_bit + log_2nd_bit  + sum_log), int = int_prop_l, dens = dens_prep_l))
  }
}


Rdeath <- function(rm, l_rm, lambda_J, rho_L, alpha, dens_mat, int_vals, point_area, k, log=FALSE){
  log_1st_bit <- log(rho_L * lambda_J * abs(l_rm[3]) / (k + 1))
  log_2nd_bit <- -alpha*int_vals[rm]*point_area
  if(sum(dens_mat[,-rm]) > 0){
    if(ncol(dens_mat)==2){
      sum_log <- sum(log(1 + dens_mat[,rm]/sum(dens_mat[,-rm])), na.rm=TRUE)
    } else {
      sum_log <-sum(log(1 + dens_mat[,rm]/apply(dens_mat[,-rm], MARGIN = 1, sum)), na.rm=TRUE)
    }
    if(log){
      
      return(- log_1st_bit - log_2nd_bit - sum_log)
    } else {
      return(1/exp(log_1st_bit + log_2nd_bit  + sum_log))
    }
    
  } else{
    if(log){
      return(-Inf)
    } else {
      return(0)
    }
  }
}

infer.PLCP <- function(spatPat, a1, b1, a2, b2, M, M_integral=33^2, prop_sigma, max_sigma_sq, W_margin,
                       sigma_sq_init, lines, birth_prob, death_prob){
  
  a <- spatPat$xDiff/2+W_margin
  
  accept <- rep(NA, M)
  sigma_sq_sample <- c(sigma_sq_init, rep(NA, M))
  a1_sample <- c(a1, rep(NA,M))
  b1_sample <- c(b1, rep(NA, M))
  a2_sample <- c(a2, rep(NA, M))
  alpha_sample <- rho_L_sample <- rep(NA, M)
  
  # simulated points for MC approximation of an integral
  # MC_points <- rHomPois(xdim=spatPat$winBor[1:2,1], ydim=spatPat$winBor[2:3,2], numPoints=M_integral)$pattern
  xs <- seq(spatPat$winBor[1,1]+spatPat$xDiff/(2*sqrt(M_integral)), spatPat$winBor[2,1]-spatPat$xDiff/(2*sqrt(M_integral)), length=sqrt(M_integral))
  ys <- seq(spatPat$winBor[1,2]+spatPat$xDiff/(2*sqrt(M_integral)), spatPat$winBor[3,2]-spatPat$xDiff/(2*sqrt(M_integral)), length=sqrt(M_integral))
  MC_points<- expand.grid(xs, ys)
  point_area <- spatPat$xDiff*spatPat$yDiff/M_integral
  # function for integrating
  
  
  # *** UPDATING *** #
  
  # plotSpatPat(spatPat)
  
  # Gibbs step - outside of the loop
  
  a1_post <- a1 + nrow(spatPat$pattern)
  b2_post <- b2 + 8*a/pi
  
  ints_cur <- int_sum(MC_points, lines, sigma_sq_sample[1])
  dens_cur <- sum_log_matrix(spatPat$pattern, lines, sigma_sq_sample[1])
  
  # inside the loop 
  
  # BEGIN LOOP
  for(iter in 1:M){
    
    # Gibbs step - continuation
    
    # a1_sample[iter+1] <- a1_sample[1] + nrow(spatPat$pattern)/nrow(lines)
    b1_sample[iter+1] <- b1_sample[1] + sum(ints_cur) * point_area
    a2_sample[iter+1] <- a2_sample[1] + nrow(lines)
    
    
    # plot(seq(0,2, length=1000), dgamma(seq(0,2, length=1000), shape=a1_post, rate=b1_sample[iter+1]))
    # plot(seq(0,4, length=1000), dgamma(seq(0,4, length=1000), shape=a2_sample[iter+1], rate=b2_post))
    
    
    
    alpha_sample[iter] <- rgamma(1, shape=a1_post, rate=b1_sample[iter+1])
    rho_L_sample[iter] <- rgamma(1, shape=a2_sample[iter+1], rate=b2_post)
    
    
    # Metropolis-Hastings updating
    # Updating sigma_sq
    
    sigma_sq_prop <- runif(1, max(0, sigma_sq_sample[iter]-prop_sigma), min(sigma_sq_sample[iter]+prop_sigma, max_sigma_sq))
    
    log_p_ratio <- dunif(sigma_sq_prop, max(0, sigma_sq_sample[iter]-prop_sigma), min(sigma_sq_sample[iter]+prop_sigma, max_sigma_sq), log=TRUE) -
      dunif(sigma_sq_sample[iter], max(0, sigma_sq_prop-prop_sigma), min(sigma_sq_sample[iter]+prop_sigma, max_sigma_sq), log=TRUE)
    
    ints_prop <- int_sum(MC_points, lines, sigma_sq_prop)
    log_int_ratio <- alpha_sample[iter] *( sum(ints_cur) - sum(ints_prop) ) * point_area
    
    dens_prop <- sum_log_matrix(spatPat$pattern, lines, sigma_sq_prop)
    if(ncol(dens_cur)==2 || !is.matrix(dens_cur)){
      prod_sum_ratio <- sum(log( apply(dens_prop, MARGIN = 1, sum) / sum(dens_cur) ), na.rm = TRUE)
    } else {
      prod_sum_ratio <- sum(log( apply(dens_prop, MARGIN = 1, sum) / apply(dens_cur, MARGIN = 1, sum) ), na.rm = TRUE)
    }
    
    
    logR <- log_p_ratio + log_int_ratio + prod_sum_ratio
    
    threshold <- runif(1, 0, 1)
    
    if(log(threshold) < logR) {
      sigma_sq_sample[iter+1] <- sigma_sq_prop
      accept[iter] <- TRUE
      ints_cur <- ints_prop
      dens_cur <- dens_prop
    } else {
      sigma_sq_sample[iter+1] <- sigma_sq_sample[iter]
      accept[iter] <- FALSE
    }
    
    # Birth-Death-Move step
    
    threshold_BDM <- runif(1, 0, 1)
    # plotSpatPat(spatPat)
    # for(i in 1:nrow(lines)){
    #   abline(-lines[i,1]*tan(lines[i,4]), tan(lines[i,4]))
    # }
    if(threshold_BDM < birth_prob){
      # BIRTH
      u_angle <- runif(1, 0, 2*pi)
      u_prop <- c(cos(u_angle), sin(u_angle))
      
      if((u_angle > 0 && u_angle <= pi/2) || (u_angle >= pi && u_angle <= 3*pi/2)){
        y_prop <- runif(1, -a/tan(u_angle) - a, a/tan(u_angle)+a)
        lambda_J <- 2*a + 2*a/tan(u_angle)
        if(y_prop > 0){
          psi <- u_angle %% pi + pi
          x <- y_prop
        } else {
          psi <- (pi-u_angle) %% pi + pi
          x <- -y_prop
        }
        if((a-x)*tan(psi)<=a) {
          continue <- TRUE
        } else {
          continue <- FALSE
        }
        
      } else {
        y_prop <- runif(1, a/tan(u_angle) - a, a-a/tan(u_angle)+2)
        lambda_J <- 2*a - 2*a/tan(u_angle)
        if(y_prop > 0){
          psi <- u_angle %% pi
          x <- y_prop
        } else {
          psi <- (pi-u_angle) %% pi 
          x <- -y_prop
        }
        if((x-a)*tan(psi)<=a) {
          continue <- TRUE
        } else {
          continue <- FALSE
        }
      }
      
      
      # plot(1,1, type="n", xlim=c(-0.6, 0.6), ylim=c(--0.6,0.6))
      # abline(-y_prop*tan(u_angle), tan(u_angle))
      if(continue){
        l_prop <- c(y_prop, u_prop, u_angle)
        logRb <- Rbirth(l_prop, points=spatPat$pattern, lines, lambda_J, rho_L = rho_L_sample[iter], 
                        alpha = alpha_sample[iter], sigma_sq = sigma_sq_sample[iter+1], MC_points, dens_mat = dens_cur, point_area=point_area,
                        log=TRUE)
        
        threshold <- runif(1,0,1)
        if(log(threshold) < logRb[[1]]) {
          lines <- rbind(lines, l_prop)
          ints_cur <- c(ints_cur, logRb$int)
          dens_cur <- cbind(dens_cur, logRb$dens)
        }
      }
      
      
    } else if(threshold_BDM < birth_prob+death_prob && nrow(lines) >= 2){
      # DEATH
      rm <- sample(1:nrow(lines),1)
      
      l_rm <- lines[rm,]
      if((l_rm[4] > 0 && l_rm[4] <= pi/2) || (l_rm[4] >= pi && l_rm[4] <= 3*pi/2)){
        lambda_J <- 2*a + 2*a/tan(l_rm[4])
      } else {
        lambda_J <- 2*a - 2*a/tan(l_rm[4])
      }
      if(nrow(lines)==2) {
        l_less <- matrix(lines[-rm,], nrow=1)
      } else {
        l_less <- lines[-rm,]
      }
      logRd <- Rdeath(rm, l_rm, lambda_J, rho_L = rho_L_sample[iter], alpha = alpha_sample[iter], 
                      dens_mat = dens_cur, int_vals = ints_cur, point_area=point_area, k = nrow(lines), log=TRUE)
      threshold <- runif(1,0,1)
      if(log(threshold) < logRd) {
        lines <- l_less
        ints_cur <- ints_cur[-rm]
        dens_cur <- dens_cur[,-rm]
      }
    } else if(nrow(lines) >= 2){
      #MOVE
      
      #birth part
      u_angle <- runif(1, 0, 2*pi)
      u_prop <- c(cos(u_angle), sin(u_angle))
      
      if((u_angle > 0 && u_angle <= pi/2) || (u_angle >= pi && u_angle <= 3*pi/2)){
        y_prop <- runif(1, -a/tan(u_angle) - a, a/tan(u_angle)+a)
        lambda_J <- 2*a + 2*a/tan(u_angle)
        if(y_prop > 0){
          psi <- u_angle %% pi + pi
          x <- y_prop
        } else {
          psi <- (pi-u_angle) %% pi + pi
          x <- -y_prop
        }
        if((a-x)*tan(psi)<=a) {
          continue <- TRUE
        } else {
          continue <- FALSE
        }
        
      } else {
        y_prop <- runif(1, a/tan(u_angle) - a, a-a/tan(u_angle)+2)
        lambda_J <- 2*a - 2*a/tan(u_angle)
        if(y_prop > 0){
          psi <- u_angle %% pi
          x <- y_prop
        } else {
          psi <- (pi-u_angle) %% pi 
          x <- -y_prop
        }
        if((x-a)*tan(psi)<=a) {
          continue <- TRUE
        } else {
          continue <- FALSE
        }
      }
      if(continue){
        l_prop <- c(y_prop, u_prop, u_angle)
        log_num <- Rbirth(l_prop, points=spatPat$pattern, lines, lambda_J, rho_L = rho_L_sample[iter], 
                          alpha = alpha_sample[iter], sigma_sq = sigma_sq_sample[iter+1], MC_points, dens_mat = dens_cur, point_area=point_area,
                          log=TRUE)
        #death part
        rm <- sample(1:nrow(lines),1)
        l_rm <- lines[rm,]
        if((l_rm[4] > 0 && l_rm[4] <= pi/2) || (l_rm[4] >= pi && l_rm[4] <= 3*pi/2)){
          lambda_J <- 2*a + 2*a/tan(l_rm[4])
        } else {
          lambda_J <- 2*a - 2*a/tan(l_rm[4])
        }
        if(nrow(lines)==2) {
          l_less <- matrix(lines[-rm,], nrow=1)
        } else {
          l_less <- lines[-rm,]
        }
        log_denom <- Rdeath(rm, l_rm, lambda_J, rho_L = rho_L_sample[iter], alpha = alpha_sample[iter], 
                            dens_mat = dens_cur, int_vals = ints_cur, point_area=point_area, k = nrow(lines), log=TRUE)
        
        
        logRm <- log_num[[1]]+log_denom
        if(log(threshold) < logRm) {
          lines <- rbind(l_less, l_prop)
          ints_cur <- c(ints_cur[-rm], log_num$int)#### ADD THE NEW ONE AS IN BIRTH
          dens_cur <- cbind(dens_cur[,-rm], log_num$dens)
        }
      }
      
    }
    
    print(iter)
    
  } # END LOOP
  return(list(sigma_sq_sample=sigma_sq_sample, b1_sample=b1_sample, a2_sample=a2_sample,
              alpha_sample=alpha_sample, rho_L_sample=rho_L_sample, acccept=accept))
}

lines <- matrix(data=c(0.6, cos(pi*0.3), sin(pi*0.3), pi*0.3,
                       0.5, cos(11*pi/6), sin(11*pi/6), 11*pi/6,
                       -0.5, cos(2*pi/11), sin(2*pi/11), 2*pi/11,
                       1.5, cos(pi/6), sin(pi/6), pi/6,
                       -1.5, cos(pi/6), sin(pi/6), pi/6,
                       0.5, cos(pi*0.8), sin(pi*0.8), pi*0.8,
                       -2.4, cos(pi*0.9), sin(pi*0.9), pi*0.9,
                       -5, cos(pi/90), sin(pi/90), pi/90,
                       0.6, cos(pi*0.3), sin(pi*0.3), pi*0.3
),ncol=4, byrow=TRUE)

  

}


# REPLICATES

{
set.seed(4894)
srPLCP <- stchRecH(PLCP, iter=2e4)
set.seed(484)
srGIBBS <- stchRecH(GIBBS, iter=2e4)
set.seed(485655)
srLGCP <- stchRecH(LGCP, iter=2e4)
set.seed(134)
tilePLCP <- singleBoot(PLCP, sqrt.n=6)
set.seed(146)
tileGIBBS <- singleBoot(GIBBS, sqrt.n=6)
set.seed(4894)
tileLGCP <- singleBoot(LGCP, sqrt.n=6)

set.seed(363)
est_LGCP_LGCP <- infer.LGCP(LGCP, "exponential")
LGCP_LGCP <- rLGCPtransformed(var = est_LGCP_LGCP[1], scale = est_LGCP_LGCP[2], model="exponential",
                              mu = est_LGCP_LGCP[3], xdim = c(-0.5, 0.5), ydim = c(-0.5, 0.5))
set.seed(46231)
est_LGCP_Thomas <- infer.Thomas(LGCP)
LGCP_Thomas <- rThomas(parentIntensity = est_LGCP_Thomas[1], meanClusterSize = est_LGCP_Thomas[3],
                       sigma = est_LGCP_Thomas[2], xdim = c(-0.5, 0.5), ydim = c(-0.5, 0.5), 
                       xScalingFactor=1, theta=theta)
set.seed(1515)
est_PLCP_Thomas <- infer.Thomas(PLCP)
PLCP_Thomas <- rThomas(parentIntensity = est_PLCP_Thomas[1], meanClusterSize = est_PLCP_Thomas[3],
                       sigma = est_PLCP_Thomas[2], xdim = c(-0.5, 0.5), ydim = c(-0.5, 0.5), 
                       xScalingFactor=1, theta=theta)
set.seed(2787)
est_Gibbs_Strauss <- infer.Strauss(GIBBS)
Gibbs_Strauss <- rmhStrauss_rep(beta = est_Gibbs_Strauss[1], gamma = est_Gibbs_Strauss[2],
                                R = est_Gibbs_Strauss[3], xdim = c(-0.5, 0.5), ydim = c(-0.5, 0.5))

set.seed(54)
est_PLCP_PLCP <- infer.PLCP(spatPat=PLCP, a1=1+1e-6, b1=1e-12, a2=1+1e-6, b2=1e-12, M=7241, M_integral=33^2,
                            prop_sigma=1e-3, max_sigma_sq=0.3, W_margin=0.05,
                            sigma_sq_init=0.0225^2, lines=lines, birth_prob=1/3, death_prob=1/3)
set.seed(75)
PLCP_PLCP <- rDirLine(lineIntensity = est_PLCP_PLCP$rho_L[7241], alongIntensity = est_PLCP_PLCP$alpha[7241],
                      sigma = sqrt(est_PLCP_PLCP$sigma_sq_sample[7241]), 
                      xdim = c(-0.5, 0.5), ydim = c(-0.5, 0.5), kappa=KofA(1), theta=pi/6)
par(mfrow=c(1,1))
plotSpatPat(PLCP_PLCP, cex=0.75)


}




{
where=-2
layout(matrix(1:16, ncol=4, byrow=TRUE), widths=c(0.55, 5,5,5))
par(mar=c(0,0,0,0))
plot.new()
title(ylab=TeX("MC w/ correct model"), line=where, cex.lab=2)
par(mar=c(1,1,1,1))

plotSpatPat(LGCP_LGCP, cex=0.75)
plot.new()
plotSpatPat(PLCP_PLCP, cex=0.75)

par(mar=c(0,0,0,0))
plot.new()
title(ylab=TeX("MC w/ misspecified model"), line=where, cex.lab=2)
par(mar=c(1,1,1,1))

plotSpatPat(LGCP_Thomas, cex=0.75)
plotSpatPat(Gibbs_Strauss, cex=0.75)
plotSpatPat(PLCP_Thomas, cex=0.75)

par(mar=c(0,0,0,0))
plot.new()
title(ylab=TeX("Nonparametric: tiling"), line=where, cex.lab=2)
par(mar=c(1,1,1,1))

plotSpatPat(tileLGCP, cex=0.75)
plotSpatPat(tileGIBBS, cex=0.75)
plotSpatPat(tilePLCP, cex=0.75)

par(mar=c(0,0,0,0))
plot.new()
title(ylab=TeX("Nonparametric: stoch. rec."), line=where+0.4, cex.lab=2)
par(mar=c(1,1,1,1))

plotSpatPat(srLGCP, cex=0.75)
plotSpatPat(srGIBBS, cex=0.75)
plotSpatPat(srPLCP, cex=0.75)

}


# CUSTOM POINTS

{

points_custom <- function(x, y, shape, col = 'black', cex = 1, ...) {
  
  if(missing(shape)) {
    points(x, y, col = col, cex = cex, ...) 
  } 
  else {
    shape <- lapply(shape, function(z) z * cex)
    Map(function(x_i, y_i) {
      a <- grconvertX(grconvertX(x_i, 'user', 'inches') + shape$x, 'inches', 'user')
      b <- grconvertY(grconvertY(y_i, 'user', 'inches') + shape$y, 'inches', 'user')
      polygon(a, b, col = col, border = col, ...)
    }, x_i = x, y_i = y)
  }
  invisible(NULL)
}

my_shape1 <- list(x = c(cos(18/180*pi), cos(54/180*pi)/3, cos(pi/2), cos(126/180*pi)/3, cos(162/180*pi),
                        cos(198/180*pi)/3, cos(234/180*pi), cos(270/180*pi)/3, cos(306/180*pi), cos(342/180*pi)/3)/20, 
                  y = c(sin(18/180*pi), sin(54/180*pi)/3, sin(pi/2), sin(126/180*pi)/3, sin(162/180*pi),
                        sin(198/180*pi)/3, sin(234/180*pi), sin(270/180*pi)/3, sin(306/180*pi), sin(342/180*pi)/3)/20)


}




# LGCP

LGCP_B <- read.csv("POW_LB.csv")
LGCP_S <- read.csv("POW_LS.csv")
POW_LOS <- LGCP_S[1:3,-1]
POW_LLS <- LGCP_S[4:6,-1]
POW_LTS <- LGCP_S[7:9,-1]
POW_LT2S <- LGCP_S[10:12,-1]
POW_LT3S <- LGCP_S[13:15,-1]
POW_LT4S <- LGCP_S[10:18,-1]
POW_LT5S <- LGCP_S[19:21,-1]
POW_LSS <- LGCP_S[22:24,-1]

POW_LOB <- LGCP_B[1:3,-1]
POW_LLB <- LGCP_B[4:6,-1]
POW_LTB <- LGCP_B[7:9,-1]
POW_LT4B <- LGCP_B[10:12,-1]
POW_LT5B <- LGCP_B[13:15,-1]
POW_LT6B <- LGCP_B[16:18,-1]
POW_LT8B <- LGCP_B[19:21,-1]
POW_LSB <- LGCP_B[22:24,-1]



{
  # layout(matrix(c(1,2,3,4,5,6), ncol=1, byrow=FALSE), heights=c(3,3,3,3,3, 2))
  layout(matrix(c(1:13, 13, 13), ncol=3, byrow=TRUE), heights=c(5,5,5,1,1), widths = c(1, 7,7))
  par(mar=c(1,2.25,1.5,1))
  
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$G_{loc,\\alpha, \\epsilon}$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_LOS[1,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_LLS[1,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_LTS[1,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_LT2S[1,], 3:0, pch=16, col=blue(255), cex=1.5)
  # points(POW_LT3S[1,], 3:0, pch=16, col=cyan(255), cex=1.5)
  # points(POW_LT4S[1,], 3:0, pch=16, col=teal(255), cex=1.5)
  points(POW_LT5S[1,], 3:0-0.03, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_LSS[1,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  title(main =TeX("$W=[-0.25, 0.25]^2$"), line=0.5, cex.main=1.5)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_LOB[1,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_LLB[1,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_LTB[1,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_LT4B[1,], 3:0, pch=16, col=blue(255), cex=1.5)
  # points(POW_LT6B[1,], 3:0, pch=16, col=cyan(255), cex=1.5)
  points(POW_LT8B[1,], 3:0-0.03, pch=16, col=teal(255), cex=1.5)
  # points(POW_LT5B[1,], 3:0, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_LSB[1,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  title(main =TeX("$W=[-0.5, 0.5]^2$"), line=0.5, cex.main=1.5)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$K_{cyl,\\alpha, \\zeta}$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_LOS[2,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_LLS[2,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_LTS[2,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_LT2S[2,], 3:0, pch=16, col=blue(255), cex=1.5)
  points(POW_LT3S[2,], 3:0-0.03, pch=16, col=teal(255), cex=1.5)
  # points(POW_LT4S[2,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_LT5S[2,], 3:0, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_LSS[2,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_LOB[2,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_LLB[2,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_LTB[2,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_LT8B[2,], 3:0, pch=16, col=blue(255), cex=1.5)
  # points(POW_LT6B[2,], 3:0, pch=16, col=cyan(255), cex=1.5)
  points(POW_LT4B[2,], 3:0-0.03, pch=16, col=teal(255), cex=1.5)
  # points(POW_LT5B[2,], 3:0, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_LSB[2,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$\\Theta-spectrum$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_LOS[3,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_LLS[3,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_LTS[3,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_LT2S[3,], 3:0, pch=16, col=blue(255), cex=1.5)
  # points(POW_LT3S[3,], 3:0, pch=16, col=cyan(255), cex=1.5)
  # points(POW_LT4S[3,], 3:0, pch=16, col=teal(255), cex=1.5)
  points(POW_LT5S[3,], 3:0-0.03, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_LSS[3,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_LOB[3,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_LLB[3,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_LTB[3,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_LT4B[3,], 3:0, pch=16, col=blue(255), cex=1.5)
  # points(POW_LT5B[3,], 3:0, pch=16, col=cyan(255), cex=1.5)
  # points(POW_LT6B[3,], 3:0, pch=16, col=teal(255), cex=1.5)
  points(POW_LT8B[3,], 3:0-0.03, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_LSB[3,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot.new()
  plot.new()
  title(xlab="Rejection rate", line=-1, cex.lab=1.25)
  plot.new()
  title(xlab="Rejection rate", line=-1, cex.lab=1.25)
  
  
  par(mar=c(0,2.25,0,1))
  plot.new()
  
  l = 0.12
  legend(x=-0.03+l, y=1, c("Oracle",  "MC w/ LGCP", "MC w/ Thomas", "Tiling", "Stoch. rec."), pch=c(NA, 16,16,16,16), pt.cex=1.5, 
         col=c(NA, orange(255),  magenta(255), teal(255), grey(255)), ncol=5, cex=1.5)
  points_custom(-0.011+l, 0.66, shape=my_shape1, col=black(255), cex=1.5)
  
}



# GIBBS

GIBBS_B <- read.csv("POW_GB.csv")
GIBBS_S <- read.csv("POW_GS.csv")
POW_GOS <- GIBBS_S[1:3,-1]
POW_GSS <- GIBBS_S[4:6,-1]
POW_GT2S <- GIBBS_S[7:9,-1]
POW_GT3S <- GIBBS_S[10:12,-1]
POW_GT4S <- GIBBS_S[13:15,-1]
POW_GT5S <- GIBBS_S[10:18,-1]
POW_GSRS <- GIBBS_S[19:21,-1]

POW_GOB <- GIBBS_B[1:3,-1]
POW_GSB <- GIBBS_B[4:6,-1]
POW_GT4B <- GIBBS_B[7:9,-1]
POW_GT5B <- GIBBS_B[10:12,-1]
POW_GT6B <- GIBBS_B[13:15,-1]
POW_GT8B <- GIBBS_B[16:18,-1]
POW_GSRB <- GIBBS_B[19:21,-1]

{
  # layout(matrix(c(1,2,3,4,5,6), ncol=1, byrow=FAGSRE), heights=c(3,3,3,3,3, 2))
  layout(matrix(c(1:13, 13, 13), ncol=3, byrow=TRUE), heights=c(5,5,5,1,1), widths = c(1, 7,7))
  par(mar=c(1,2.25,1.5,1))
  
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$G_{loc,\\alpha, \\epsilon}$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_GOS[1,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_GSS[1,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  points(POW_GT2S[1,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_GT3S[1,], 3:0, pch=16, col=cyan(255), cex=1.5)
  # points(POW_GT4S[1,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_GT5S[1,], 3:0-0.03, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_GSRS[1,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  title(main =TeX("$W=[-0.25, 0.25]^2$"), line=0.5, cex.main=1.5)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_GOB[1,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_GSB[1,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  points(POW_GT4B[1,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_GT6B[1,], 3:0, pch=16, col=cyan(255), cex=1.5)
  # points(POW_GT8B[1,], 3:0-0.03, pch=16, col=teal(255), cex=1.5)
  # points(POW_GT5B[1,], 3:0, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_GSRB[1,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  title(main =TeX("$W=[-0.5, 0.5]^2$"), line=0.5, cex.main=1.5)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$K_{cyl,\\alpha, \\zeta}$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_GOS[2,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_GSS[2,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  points(POW_GT2S[2,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_GT3S[2,], 3:0-0.03, pch=16, col=teal(255), cex=1.5)
  # points(POW_GT4S[2,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_GT5S[2,], 3:0, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_GSRS[2,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_GOB[2,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_GSB[2,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_GT8B[2,], 3:0, pch=16, col=blue(255), cex=1.5)
  points(POW_GT6B[2,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_GT4B[2,], 3:0-0.03, pch=16, col=teal(255), cex=1.5)
  # points(POW_GT5B[2,], 3:0, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_GSRB[2,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$\\Theta-spectrum$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_GOS[3,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_GSS[3,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_GT2S[3,], 3:0, pch=16, col=blue(255), cex=1.5)
  # points(POW_GT3S[3,], 3:0, pch=16, col=cyan(255), cex=1.5)
  # points(POW_GT4S[3,], 3:0, pch=16, col=teal(255), cex=1.5)
  points(POW_GT5S[3,], 3:0-0.03, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_GSRS[3,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_GOB[3,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_GSB[3,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  points(POW_GT4B[3,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_GT5B[3,], 3:0, pch=16, col=cyan(255), cex=1.5)
  # points(POW_GT6B[3,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_GT8B[3,], 3:0-0.03, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_GSRB[3,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot.new()
  plot.new()
  title(xlab="Rejection rate", line=-1, cex.lab=1.25)
  plot.new()
  title(xlab="Rejection rate", line=-1, cex.lab=1.25)
  
  
  par(mar=c(0,2.25,0,1))
  plot.new()
  
  l = 0.2375
  legend(x=-0.03+l, y=1, c("Oracle",  "MC w/ Strauss", "Tiling", "Stoch. rec."), pch=c(NA, 16,16,16), pt.cex=1.5, 
         col=c(NA, magenta(255), teal(255), grey(255)), ncol=4, cex=1.5)
  points_custom(-0.011+l, 0.66, shape=my_shape1, col=black(255), cex=1.5)
  
}




# PLCP

PLCP_B <- read.csv("POW_PB.csv")
PLCP_S <- read.csv("POW_PS.csv")
POW_POS <- PLCP_S[1:3,-1]
POW_PPS <- PLCP_S[4:6,-1]
POW_PTS <- PLCP_S[7:9,-1]
POW_PT2S <- PLCP_S[10:12,-1]
POW_PT3S <- PLCP_S[13:15,-1]
POW_PT4S <- PLCP_S[10:18,-1]
POW_PT5S <- PLCP_S[19:21,-1]
POW_PSS <- PLCP_S[22:24,-1]

POW_POB <- PLCP_B[1:3,-1]
POW_PPB <- PLCP_B[4:6,-1]
POW_PTB <- PLCP_B[7:9,-1]
POW_PT4B <- PLCP_B[10:12,-1]
POW_PT5B <- PLCP_B[13:15,-1]
POW_PT6B <- PLCP_B[16:18,-1]
POW_PT8B <- PLCP_B[19:21,-1]
POW_PSB <- PLCP_B[22:24,-1]

{
  # layout(matrix(c(1,2,3,4,5,6), ncol=1, byrow=FALSE), heights=c(3,3,3,3,3, 2))
  layout(matrix(c(1:13, 13, 13), ncol=3, byrow=TRUE), heights=c(5,5,5,1,1), widths = c(1, 7,7))
  par(mar=c(1,2.25,1.5,1))
  
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$G_{loc,\\alpha, \\epsilon}$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_POS[1,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_PPS[1,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_PTS[1,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_PT2S[1,], 3:0, pch=16, col=blue(255), cex=1.5)
  # points(POW_PT3S[1,], 3:0, pch=16, col=cyan(255), cex=1.5)
  # points(POW_PT4S[1,], 3:0, pch=16, col=teal(255), cex=1.5)
  points(POW_PT5S[1,], 3:0-0.03, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_PSS[1,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  title(main =TeX("$W=[-0.25, 0.25]^2$"), line=0.5, cex.main=1.5)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_POB[1,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_PPB[1,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_PTB[1,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_PT4B[1,], 3:0, pch=16, col=blue(255), cex=1.5)
  # points(POW_PT6B[1,], 3:0, pch=16, col=cyan(255), cex=1.5)
  # points(POW_PT8B[1,], 3:0-0.03, pch=16, col=teal(255), cex=1.5)
  points(POW_PT5B[1,], 3:0, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_PSB[1,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  title(main =TeX("$W=[-0.5, 0.5]^2$"), line=0.5, cex.main=1.5)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$K_{cyl,\\alpha, \\zeta}$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_POS[2,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_PPS[2,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_PTS[2,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_PT2S[2,], 3:0, pch=16, col=blue(255), cex=1.5)
  points(POW_PT3S[2,], 3:0-0.03, pch=16, col=teal(255), cex=1.5)
  # points(POW_PT4S[2,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_PT5S[2,], 3:0, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_PSS[2,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_POB[2,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_PPB[2,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_PTB[2,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_PT8B[2,], 3:0, pch=16, col=blue(255), cex=1.5)
  # points(POW_PT6B[2,], 3:0, pch=16, col=cyan(255), cex=1.5)
  points(POW_PT4B[2,], 3:0-0.03, pch=16, col=teal(255), cex=1.5)
  # points(POW_PT5B[2,], 3:0, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_PSB[2,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$\\Theta-spectrum$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_POS[3,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_PPS[3,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_PTS[3,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  # points(POW_PT2S[3,], 3:0, pch=16, col=blue(255), cex=1.5)
  points(POW_PT3S[3,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_PT4S[3,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_PT5S[3,], 3:0-0.03, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_PSS[3,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points_custom(POW_POB[3,], 3:0, shape=my_shape1, col=black(255), cex=1.5)
  points(POW_PPB[3,], 3:0+0.09, pch=16, col=orange(255), cex=1.5)
  points(POW_PTB[3,], 3:0+0.03, pch=16, col=magenta(255), cex=1.5)
  points(POW_PT4B[3,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_PT5B[3,], 3:0, pch=16, col=cyan(255), cex=1.5)
  # points(POW_PT6B[3,], 3:0, pch=16, col=teal(255), cex=1.5)
  # points(POW_PT8B[3,], 3:0-0.03, pch=16, col=teal(), cex=1.5) # BEST
  points(POW_PSB[3,], 3:0-0.09, pch=16, col=grey(255), cex=1.5)
  
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot.new()
  plot.new()
  title(xlab="Rejection rate", line=-1, cex.lab=1.25)
  plot.new()
  title(xlab="Rejection rate", line=-1, cex.lab=1.25)
  
  
  par(mar=c(0,2.25,0,1))
  plot.new()
  
  l = 0.12
  legend(x=-0.03+l, y=1, c("Oracle",  "MC w/ PLCP", "MC w/ Thomas", "Tiling", "Stoch. rec."), pch=c(NA, 16,16,16,16), pt.cex=1.5, 
         col=c(NA, orange(255),  magenta(255), teal(255), grey(255)), ncol=5, cex=1.5)
  points_custom(-0.011+l, 0.66, shape=my_shape1, col=black(255), cex=1.5)
  
}





###############
#TILING


{
  # layout(matrix(c(1,2,3,4,5,6), ncol=1, byrow=FALSE), heights=c(3,3,3,3,3, 2))
  layout(matrix(c(1:13, 13, 13), ncol=3, byrow=TRUE), heights=c(5,5,5,1,1.5), widths = c(1, 7,7))
  par(mar=c(1,2.25,1.5,1))
  
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$G_{loc,\\alpha, \\epsilon}$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points(POW_LT2S[1,], 3:0+0.17, pch=0, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT3S[1,], 3:0+0.07, pch=1, col=magenta(255), cex=1.5, lwd=2)
  points(POW_LT4S[1,], 3:0-0.03, pch=2, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT5S[1,], 3:0-0.13, pch=3, col=magenta(), cex=1.25, lwd=2) # BEST
  points(POW_GT2S[1,], 3:0+0.15, pch=0, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT3S[1,], 3:0+0.05, pch=1, col=teal(255), cex=1.5, lwd=2)
  points(POW_GT4S[1,], 3:0-0.05, pch=2, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT5S[1,], 3:0-0.15, pch=3, col=teal(), cex=1.25, lwd=2) # BEST
  points(POW_PT2S[1,], 3:0+0.13, pch=0, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT3S[1,], 3:0+0.03, pch=1, col=orange(255), cex=1.5, lwd=2)
  points(POW_PT4S[1,], 3:0-0.07, pch=2, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT5S[1,], 3:0-0.17, pch=3, col=orange(), cex=1.25, lwd=2) # BEST
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  title(main =TeX("$W=[-0.25, 0.25]^2$"), line=0.5, cex.main=1.5)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points(POW_LT4B[1,], 3:0+0.17, pch=0, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT5B[1,], 3:0+0.07, pch=1, col=magenta(255), cex=1.5, lwd=2)
  points(POW_LT6B[1,], 3:0-0.03, pch=2, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT8B[1,], 3:0-0.13, pch=3, col=magenta(), cex=1.25, lwd=2) # BEST
  points(POW_GT4B[1,], 3:0+0.15, pch=0, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT5B[1,], 3:0+0.05, pch=1, col=teal(255), cex=1.5, lwd=2)
  points(POW_GT6B[1,], 3:0-0.05, pch=2, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT8B[1,], 3:0-0.15, pch=3, col=teal(), cex=1.25, lwd=2) # BEST
  points(POW_PT4B[1,], 3:0+0.13, pch=0, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT5B[1,], 3:0+0.03, pch=1, col=orange(255), cex=1.5, lwd=2)
  points(POW_PT6B[1,], 3:0-0.07, pch=2, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT8B[1,], 3:0-0.17, pch=3, col=orange(), cex=1.25, lwd=2) 
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  title(main =TeX("$W=[-0.5, 0.5]^2$"), line=0.5, cex.main=1.5)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$K_{cyl,\\alpha, \\zeta}$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points(POW_LT2S[2,], 3:0+0.17, pch=0, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT3S[2,], 3:0+0.07, pch=1, col=magenta(255), cex=1.5, lwd=2)
  points(POW_LT4S[2,], 3:0-0.03, pch=2, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT5S[2,], 3:0-0.13, pch=3, col=magenta(), cex=1.25, lwd=2) # BEST
  points(POW_GT2S[2,], 3:0+0.15, pch=0, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT3S[2,], 3:0+0.05, pch=1, col=teal(255), cex=1.5, lwd=2)
  points(POW_GT4S[2,], 3:0-0.05, pch=2, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT5S[2,], 3:0-0.15, pch=3, col=teal(), cex=1.25, lwd=2) # BEST
  points(POW_PT2S[2,], 3:0+0.13, pch=0, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT3S[2,], 3:0+0.03, pch=1, col=orange(255), cex=1.5, lwd=2)
  points(POW_PT4S[2,], 3:0-0.07, pch=2, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT5S[2,], 3:0-0.15, pch=3, col=orange(), cex=1.25, lwd=2) # BEST
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points(POW_LT4B[2,], 3:0+0.17, pch=0, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT5B[2,], 3:0+0.07, pch=1, col=magenta(255), cex=1.5, lwd=2)
  points(POW_LT6B[2,], 3:0-0.03, pch=2, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT8B[2,], 3:0-0.13, pch=3, col=magenta(), cex=1.25, lwd=2) # BEST
  points(POW_GT4B[2,], 3:0+0.15, pch=0, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT5B[2,], 3:0+0.05, pch=1, col=teal(255), cex=1.5, lwd=2)
  points(POW_GT6B[2,], 3:0-0.05, pch=2, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT8B[2,], 3:0-0.15, pch=3, col=teal(), cex=1.25, lwd=2) # BEST
  points(POW_PT4B[2,], 3:0+0.13, pch=0, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT5B[2,], 3:0+0.03, pch=1, col=orange(255), cex=1.5, lwd=2)
  points(POW_PT6B[2,], 3:0-0.07, pch=2, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT8B[2,], 3:0-0.17, pch=3, col=orange(), cex=1.25, lwd=2) 
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  par(mar=c(1,2.25,1.5,1))
  plot.new()
  title(ylab=TeX("$\\Theta-spectrum$"), line=0.6, cex.lab=1.4)
  title(ylab="Anisotropy a", line=-0.7, cex.lab=1.25)
  
  
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points(POW_LT2S[3,], 3:0+0.17, pch=0, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT3S[3,], 3:0+0.07, pch=1, col=magenta(255), cex=1.5, lwd=2)
  points(POW_LT4S[3,], 3:0-0.03, pch=2, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT5S[3,], 3:0-0.13, pch=3, col=magenta(), cex=1.25, lwd=2) # BEST
  points(POW_GT2S[3,], 3:0+0.15, pch=0, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT3S[3,], 3:0+0.05, pch=1, col=teal(255), cex=1.5, lwd=2)
  points(POW_GT4S[3,], 3:0-0.05, pch=2, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT5S[3,], 3:0-0.15, pch=3, col=teal(), cex=1.25, lwd=2) # BEST
  points(POW_PT2S[3,], 3:0+0.13, pch=0, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT3S[3,], 3:0+0.03, pch=1, col=orange(255), cex=1.5, lwd=2)
  points(POW_PT4S[3,], 3:0-0.07, pch=2, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT5S[3,], 3:0-0.17, pch=3, col=orange(), cex=1.25, lwd=2) # BEST
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot(c(0,1), c(-0.25,3.25), type="n", yaxt="n", xaxt="n", frame=FALSE, xlab="", ylab="", main="")
  lines(c(0.05,0.05), c(-0.25,3.25), col="grey", lty=1, lwd=2)
  lines(c(0.2,0.2), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.4,0.4), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.6,0.6), c(-0.25,3.25), col="grey", lty=2)
  lines(c(0.8,0.8), c(-0.25,3.25), col="grey", lty=2)
  lines(c(1,1), c(-0.25,3.25), col="grey", lty=2)
  points(POW_LT4B[3,], 3:0+0.17, pch=0, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT5B[3,], 3:0+0.07, pch=1, col=magenta(255), cex=1.5, lwd=2)
  points(POW_LT6B[3,], 3:0-0.03, pch=2, col=magenta(255), cex=1.25, lwd=2)
  points(POW_LT8B[3,], 3:0-0.13, pch=3, col=magenta(), cex=1.25, lwd=2) # BEST
  points(POW_GT4B[3,], 3:0+0.15, pch=0, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT5B[3,], 3:0+0.05, pch=1, col=teal(255), cex=1.5, lwd=2)
  points(POW_GT6B[3,], 3:0-0.05, pch=2, col=teal(255), cex=1.25, lwd=2)
  points(POW_GT8B[3,], 3:0-0.15, pch=3, col=teal(), cex=1.25, lwd=2) # BEST
  points(POW_PT4B[3,], 3:0+0.13, pch=0, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT5B[3,], 3:0+0.03, pch=1, col=orange(255), cex=1.5, lwd=2)
  points(POW_PT6B[3,], 3:0-0.07, pch=2, col=orange(255), cex=1.25, lwd=2)
  points(POW_PT8B[3,], 3:0-0.17, pch=3, col=orange(), cex=1.25, lwd=2) 
  
  # title(ylab="Anisotropy a", line=1.25, cex.lab=1.1)
  axis(2, at=c(3,2,1,0), labels=c("1", "0.8", "0.6", "0.4"),
       pos=-0.03, las=2)
  axis(1, at=c(0,0.05,0.2,0.4,0.6,0.8,1), labels=c("0", TeX("$\\alpha^{(T)}$"), "0.2", "0.4", "0.6", "0.8", "1"))
  
  plot.new()
  plot.new()
  title(xlab="Rejection rate", line=-1, cex.lab=1.25)
  plot.new()
  title(xlab="Rejection rate", line=-1, cex.lab=1.25)
  
  
  par(mar=c(0,2.25,0,1))
  plot.new()
  
  l = 0.18
  legend(x=-0.03+l, y=1, c("",  "lowest  ", "LGCP    ", "second lowest    ", "Gibbs    ", "second highest    ", "PLCP    ", "highest  "), 
         pch=c(NA, 0, 16, 1, 16, 2, 16, 3), pt.cex=c(NA, 1.25, 1.5, 1.5, 1.5, 1.25, 1.5, 1.25),
         pt.lwd = c(NA, 2, NA, 2, NA, 2, NA, 2 ),
         col=c(NA, black(), magenta(255), black(), teal(255), black(), orange(255), black()), ncol=4, cex=1.25)

}




###############
# Summary tables


{
# KEST


allMCB3 <- rbind(POW_LLB[2,], POW_LTB[2,], POW_GSB[2,],POW_PPB[2,], POW_PTB[2,])
allMCB3
apply(allMCB3, 2, mean)
apply(allMCB3, 2, function(x) mean(abs(x-0.05)))
apply(allMCB3, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))




allMCS3 <- rbind(POW_LLS[2,], POW_LTS[2,], POW_GSS[2,], POW_LLS[2,], POW_PTS[2,])
allMCS3
apply(allMCS3, 2, mean)

apply(allMCS3, 2, function(x) mean(abs(x-0.05)))
apply(allMCS3, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))


allMCB4 <- rbind(POW_LTB[2,], POW_GSB[2,], POW_PTB[2,])
allMCB4
apply(allMCB4, 2, mean)

apply(allMCB4, 2, function(x) mean(abs(x-0.05)))
apply(allMCB4, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))



allMCS4 <- rbind(POW_LTS[2,], POW_GSS[2,], POW_PTS[2,])
allMCS4
apply(allMCS4, 2, mean)
apply(allMCS4, 2, function(x) mean(abs(x-0.05)))
apply(allMCS4, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))



aLLSB <- rbind(POW_LSB[2,],  POW_GSRB[2,], POW_PSB[2,])
aLLSB
apply(aLLSB, 2, mean)
apply(aLLSB, 2, function(x) mean(abs(x-0.05)))
apply(aLLSB, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))


aLLSS <- rbind(POW_LSS[2,],  POW_GSRS[2,], POW_PSS[2,])
aLLSS
apply(aLLSS, 2, mean)
apply(aLLSS, 2, function(x) mean(abs(x-0.05)))
apply(aLLSS, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))


#optimal
allTOB <- rbind(POW_LT4B[2,],  POW_GT6B[2,], POW_PT4B[2,])
allTOB
apply(allTOB, 2, mean)
apply(allTOB, 2, function(x) mean(abs(x-0.05)))
apply(allTOB, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))


allTOS <- rbind(POW_LT3S[2,],  POW_GT2S[2,], POW_PT3S[2,])
allTOS
apply(allTOS, 2, mean)
apply(allTOS, 2, function(x) mean(abs(x-0.05)))
apply(allTOS, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))





############################# same for G


allMCB3 <- rbind(POW_LLB[1,], POW_LTB[1,], POW_GSB[1,],POW_PPB[1,], POW_PTB[1,])
allMCB3
apply(allMCB3, 2, mean)
apply(allMCB3, 2, function(x) mean(abs(x-0.05)))
apply(allMCB3, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))



allMCS3 <- rbind(POW_LLS[1,], POW_LTS[1,], POW_GSS[1,], POW_LLS[1,], POW_PTS[1,])
allMCS3
apply(allMCS3, 2, mean)

apply(allMCS3, 2, function(x) mean(abs(x-0.05)))
apply(allMCS3, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))


allMCB4 <- rbind(POW_LTB[1,], POW_GSB[1,], POW_PTB[1,])
allMCB4
apply(allMCB4, 2, mean)

apply(allMCB4, 2, function(x) mean(abs(x-0.05)))
apply(allMCB4, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))



allMCS4 <- rbind(POW_LTS[1,], POW_GSS[1,], POW_PTS[1,])
allMCS4
apply(allMCS4, 2, mean)
apply(allMCS4, 2, function(x) mean(abs(x-0.05)))
apply(allMCS4, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))


aLLSB <- rbind(POW_LSB[1,],  POW_GSRB[1,], POW_PSB[1,])
aLLSB
apply(aLLSB, 2, mean)
apply(aLLSB, 2, function(x) mean(abs(x-0.05)))
apply(aLLSB, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))


aLLSS <- rbind(POW_LSS[1,],  POW_GSRS[1,], POW_PSS[1,])
aLLSS
apply(aLLSS, 2, mean)
apply(aLLSS, 2, function(x) mean(abs(x-0.05)))
apply(aLLSS, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))


#optimal
allTOB <- rbind(POW_LT8B[1,],  POW_GT4B[1,], POW_PT5B[1,])
allTOB
apply(allTOB, 2, mean)
apply(allTOB, 2, function(x) mean(abs(x-0.05)))
apply(allTOB, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))
# 
# 
allTOS <- rbind(POW_LT5S[1,],  POW_GT2S[1,], POW_PT5S[1,])
allTOS
apply(allTOS, 2, mean)
apply(allTOS, 2, function(x) mean(abs(x-0.05)))
apply(allTOS, 2, function(x) sum(abs(x[x>0.05]-0.05))/length(x))









############
# Ambrosia


data <- read.table("doi_10.5063_AA_connolly.206.1-DATA.data", header = TRUE)
uq <- unique(data[,1])
for(i in uq){
  ixs <- which(data[,1]==i)
  if(length(ixs) > 1) data <- data[-ixs[-1],]
}
window <- simWind(xdim=c(-50,50), ydim=c(-50,50))
Ambrosia <-list(pattern=cbind(as.numeric(data[,2])-50, as.numeric(data[,3])-50), winBor=window[[2]],
                lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff)


estLGCP <- infer.LGCP(Ambrosia, model)
sigma_sq_sample <- read.csv("MCMC\\sigma_sq_sample.csv")[,-1]
alpha_sample <- read.csv("MCMC\\alpha_sample.csv")[,-1]
rho_L_sample <- read.csv("MCMC\\rho_L_sample.csv")[,-1]
sigma_sqs <- sigma_sq_sample[6501:12501]
alphas <- alpha_sample[6501:12501]
rho_Ls <- rho_L_sample[6501:12501]
set.seed(57)
index <- sample(1:6000, 1)


{
set.seed(75)
layout(matrix(c(1,2,3,3,4,4,5,6,7,7,8,8,9,9,6,10,11,11,12,12,13), nrow=3, byrow=TRUE),widths = c(0.5, rep(2,6)))

par(mar=rep(0,4))
plot.new()
title(ylab="Monte Carlo (MC)", line=-3, cex.lab=1.5)
par(mar=c(2,1,1,1))
plot.new()
plotSpatPat(rLGCPtransformed("exponential", var= estLGCP[1], scale=estLGCP[2], mu=estLGCP[3], xdim=c(-50,50), ydim=c(-50,50)), cex=0.3)
title(xlab=TeX("LGCP"), line=1, cex.lab=1.5)
plotSpatPat(rDirLine(lineIntensity=rho_Ls[index], alongIntensity=alphas[index], theta=pi/6, xdim=c(-50, 50), ydim=c(-50,50), sigma=sqrt(sigma_sqs[index]), kappa=KofA(1)), cex=0.3)
title(xlab=TeX("PLCP"), line=1, cex.lab=1.5)
plot.new()

par(mar=rep(0,4))
plot.new()
title(ylab="Nonparametric: tiling", line=-2.75, cex.lab=1.5)
par(mar=c(2,1,1,1))
plotSpatPat(singleBoot(Ambrosia,sqrt.n=4), cex=0.3)
title(xlab=TeX("$\\N_{tile}=16$"), line=1, cex.lab=1.5)
plotSpatPat(singleBoot(Ambrosia,sqrt.n=5), cex=0.3)
title(xlab=TeX("$\\N_{tile}=25$"), line=1, cex.lab=1.5)
plotSpatPat(singleBoot(Ambrosia,sqrt.n=6), cex=0.3)
title(xlab=TeX("$\\N_{tile}=36$"), line=1, cex.lab=1.5)
plot.new()
plotSpatPat(singleBoot(Ambrosia,sqrt.n=7), cex=0.3)
title(xlab=TeX("$\\N_{tile}=49$"), line=1, cex.lab=1.5)
plotSpatPat(singleBoot(Ambrosia,sqrt.n=8), cex=0.3)
title(xlab=TeX("$\\N_{tile}=64$"), line=1, cex.lab=1.5)
plot.new()


}
