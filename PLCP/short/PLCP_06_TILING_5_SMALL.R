library(spatstat)
library(circular)
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



# ANGLES MATRIX AUX FUNCTION
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

# ANISOTROPY TESTING PROCEDURE GIVEN ORIGINAL AND REPLICATED VALUES FOR COMPUTATION OF TEST STATISTICS
testAniso <- function(OGv, Vs, DSS){
  if(!(DSS%in%c("G","K","T","P"))) stop("DSS must be either 'G', 'K', 'T', or 'P'.")
  reps <- ncol(Vs)
  mHat <- apply(Vs, 1, mean, na.rm=TRUE)
  vector <- OGv-mHat
  TstatRep <- rep(NA, reps)
  if(!(DSS=="G")){
    CHatDiag <- diag(1/apply(Vs, 1, var, na.rm=TRUE))
    CHatDiag[CHatDiag==Inf] <- 0 
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
  p.val <- sum(as.numeric(Tstat) < TstatRep, na.rm=TRUE)/sum(!is.na(TstatRep))
  return(list(p.val=p.val, Tstat=Tstat, TstatRep=TstatRep))
}

# CYLINDRICAL K-FUNCTION WITH OPTIONAL  CONTRIBUTIONS FOR MPM
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

# LOCAL DIRECTIONAL NEAREST NEIGHBOUR DISTANCE DISTRIBUTION
GLoc<- function(spatPat, alpha, epsilon, r=seq(0,1, length=100), 
                distance=NA, angles=NA, retMatrices=FALSE, contribution=FALSE) {
  if(!is.matrix(distance)){
    distance <- (as.matrix(dist(x=spatPat[[1]], method="euclidean")))
  }
  borders <- cbind(spatPat[[1]][,1]-spatPat$winBor[1,1], spatPat$winBor[2,1]-spatPat[[1]][,1],
                   spatPat[[1]][,2]-spatPat$winBor[1,2], spatPat$winBor[3,2]-spatPat[[1]][,2])
  distBorder <- apply(borders, 1, min)
  
  distQual <- matrix(NA, nrow=nrow(distance), ncol=ncol(distance))
  if(!is.matrix(angles)){
    angles <- angleMatrix(spatPat[[1]])
  }
  for(i in 1:nrow(spatPat[[1]])){
    TF <- ((angles[i,] > alpha-epsilon & angles[i,] < alpha+epsilon ) |
             (angles[i,] > alpha+pi-epsilon & angles[i,] < alpha+pi+epsilon ))
    distQual[i, TF] <- distance[i, TF]
  }
  diag(distQual) <- NA
  minDist <- apply(distQual, 1, function(vec){ifelse(all(is.na(vec)), NA, min(vec, na.rm=TRUE)) })
  
  if(alpha-epsilon>=0 && alpha+epsilon<=pi/2){
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
                    spatPat[[1]][,2]>spatPat[[2]][1,2]+yMar*minDist & spatPat[[1]][,2]<spatPat[[2]][3,2]-yMar*minDist)
  if(length(border)>0){
    for(i in 1:length(r)){
      contributions[border,i] <- (minDist[border]<r[i])/volumes[border]
      contributions[-border,i] <- 0
      contributions[-NAs,i]
      
    }
    numerator <- apply(contributions, 2, sum, na.rm=TRUE)
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
  observedLine <- observedLine[-1,]
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


# PARAMETERS
theta = pi/6
xScalingFactor = 0.6
lineIntensity=16
alongIntensity= 25
sigma=0.015
xdim = ydim = c(-0.5, 0.5)/2

N <- 1000
M <- 1000

Gepsilon <- pi/4
Kaspect <- 0.15


sqrt.n <- 5


r.grid.W <- seq(0, 0.05, length=37)[-1]



registerDoParallel(cores=8)





set.seed(3586)
seeds <- sample(1e6, 1e3)

output <- foreach(l=1:N, .combine=cbind, .export = c("simWind", "angleMatrix","xScalingFactor", "rDirLine", "KofA", "Lest", "sqrt.n",
                                                     "testAniso","cylK", "GLoc",  "stchRecH", "rvonmises",
                                                     "owin", "xdim", "ydim", "lineIntensity", "alongIntensity", "sigma", 
                                                     "theta", "M", "Gepsilon", "Kaspect", 
                                                     "as.ppp", "r.grid.W","boot.test_GK"),
                  .multicombine = TRUE, .maxcombine = N) %dopar% {
                    set.seed(seeds[l])
                    pat <- rDirLine(lineIntensity=lineIntensity, alongIntensity=alongIntensity, theta=theta, xdim=xdim, ydim=ydim,
                                    sigma=sigma, kappa=KofA(xScalingFactor))
                    
                    secDir <- (theta+pi/2)%%pi
                    
                    
                    Gpref <- GLoc(spatPat=pat, alpha=theta, epsilon=Gepsilon, r=r.grid.W, contribution=FALSE, retMatrices=TRUE)
                    Gsec <- GLoc(spatPat=pat, alpha=secDir, epsilon=Gepsilon, r=r.grid.W, contribution=FALSE)
                    Gvec.W <- Gpref$Gloc-Gsec
                    
                    
                    Kpref <- cylK(spatPat=pat, alpha=theta, aspect=Kaspect, r=r.grid.W, contribution=FALSE)
                    Ksec <- cylK(spatPat=pat, alpha=secDir, aspect=Kaspect, r=r.grid.W, contribution=FALSE)
                    Kvec.W <- Kpref-Ksec
                    
                    
                    
                    
                    vTiles <- boot.test_GK(spatPat=pat, sqrt.n=sqrt.n, reps = M, r.grid.W=r.grid.W, 
                                           prefDir=theta, Gepsilon=Gepsilon, Kaspect=Kaspect)
                    
                    
                    test.GW <- testAniso(Gvec.W, vTiles$vGtile.W, "G")
                    test.KW <- testAniso(Kvec.W, vTiles$vKtile.W, "K")
                    
                    print(l)
                    c(test.GW$p.val, test.GW$Tstat, test.GW$TstatRep,
                      test.KW$p.val, test.KW$Tstat, test.KW$TstatRep)
                  }

write.csv(output, "PLCP_06_TILING_5_SMALL_short.csv")