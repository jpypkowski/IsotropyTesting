library(spatstat)
library(circular)
library(doParallel)
library(mvtnorm)

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
  if(is.na(numPoints)) numPoints <- rpois(1, intensity*winSize)
  x <- runif(numPoints, min=xdim[1], max=xdim[2])
  y <- runif(numPoints, min=ydim[1], max=ydim[2])
  return(list(pattern = cbind(x, y), winBor=window[[2]],
              lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff))
}


# THOMAS PROCESS SIMULATION
# FOR BOTH ISOTROPIC AND ANISOTROPIC CASE

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



# INFERENCE FOR THOMAS PROCESS

infer.Thomas <- function(pattern, startpar=c(kappa=10,scale=0.1)){
  estimates <- rep(NA, 3)
  xrange <- c(pattern[[2]][1,1], pattern[[2]][3,1])
  yrange <- c(pattern[[2]][2,2], pattern[[2]][3,2])
  estimates[1:2] <- thomas.estpcf(as.ppp(pattern$pattern, owin(xrange=xrange, yrange=yrange)), startpar=startpar)$par
  estimates[2] <- sqrt(estimates[2])
  estimates[3] <- nrow(pattern$pattern)/estimates[1]/(diff(xrange)*diff(yrange))
  return(estimates)
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

# PERIODOGRAM WITH OPTINAL MPM CONTRIBUTIONS
periodogram <- function(spatPat, p, contribution=FALSE){
  p1 <- p
  p2 <- p
  omega1 <- 2*pi*p1/spatPat$xDiff
  omega2 <- 2*pi*p2/spatPat$yDiff
  grid <- as.matrix(expand.grid(omega1, omega2))
  area <- spatPat$xDiff*spatPat$yDiff
  if(contribution){
    contrib <- matrix(NA, ncol=nrow(grid), nrow=nrow(spatPat$pattern))
    contrib.og <- exp(-1i*grid%*%t(spatPat[[1]]))/sqrt(area)
    contrib.cc <- Re(contrib.og)-1i*Im(contrib.og)
    for(i in 1:nrow(grid)){
      out <- outer(contrib.og[i,], contrib.cc[i,], "*")
      contrib[,i] <- Re(apply(out,1,sum)+apply(out,2,sum))/2
    }
    sdf <- apply(contrib, 2, sum, na.rm=TRUE)
    
  } else {
    DFT <- rowSums(exp(-1i*grid%*%t(spatPat[[1]])))/sqrt(area)
    sdf <- Re(DFT)^2 + Im(DFT)^2
  }
  zero1 <- which(omega1==0)
  zero2 <- which(omega2==0)
  matrixSDF <- matrix(NA, nrow=length(omega1), ncol=length(omega2))
  angles <- matrix(NA, nrow=length(omega1), ncol=length(omega2))
  ang <- atan((grid[,2]/grid[,1]))
  for(i in 1:length(omega1)){
    matrixSDF[i,] <- sdf[((i-1)*length(omega2)+1):((i)*length(omega2))]
    angles[i,] <- ang[((i-1)*length(omega2)+1):((i)*length(omega2))]
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
  if(contribution){
    return(list(sdf = sdf, omega1=omega1, omega2=omega2, p1, p2, angles=angles,
                TContrib = contrib))
  } else {
    return(list(sdf = sdf, omega1=omega1, omega2=omega2, p1, p2, angles=angles))
  }
  
}

# THETA SPECTRUM
thetaSpectrum <- function(prdgrm, thetas, epsilon){ # EPSILON IS h in our specification
  spectrum <- rep(NA, length(thetas))
  for(i in 1:length(thetas)){
    consider <- (abs((prdgrm$angles-thetas[i])%%pi) < epsilon)
    spectrum[i] <- sum(prdgrm$sdf[consider],
                       na.rm=TRUE)/ sum(rowSums(consider, na.rm=TRUE))
  }
  return(spectrum)                 
}

# PARAMETERS
theta = pi/6
xScalingFactor = 1
lineIntensity=16
alongIntensity= 25
sigma=0.015
xdim = ydim = c(-0.5, 0.5)

N <- 1000
M <- 1000

grid_length <- 36
Gepsilon <- pi/4
Kaspect <- 0.15
Tepsilon <- pi*7.5/180
Tp <- -15:15

theta.grid <- seq(0, pi, length=37)[-37]
r.grid.W <- seq(0, 0.25, length=37)[-1]
r.grid.L <- seq(0, 0.05, length=37)[-1]



registerDoParallel(cores=16)




set.seed(6738)
seeds <- sample(1e6, 1e3)
output <- foreach(l=1:N, .combine=cbind, .export = c("simWind", "angleMatrix","xScalingFactor", "rDirLine", "KofA", "Lest", "rmvnorm",
                                                     "testAniso","cylK", "GLoc", "rvonmises",
                                                     "owin", "xdim", "ydim", "lineIntensity", "alongIntensity", "sigma", "thomas.estK",
                                                     "theta", "r.grid.L", "M", "Gepsilon", "Kaspect", 
                                                     "as.ppp", "r.grid.W", "rHomPois", "rThomas", "infer.Thomas"),
                  .multicombine = TRUE, .maxcombine = N) %dopar% {
                    set.seed(seeds[l])
                    pat <- rDirLine(lineIntensity=lineIntensity, alongIntensity=alongIntensity, theta=theta, xdim=xdim, ydim=ydim,
                                    sigma=sigma, kappa=KofA(xScalingFactor))
                   
                    secDir <- (theta+pi/2)%%pi
                    
                    est.param <- infer.Thomas(pat)
                    
                    
                    Gpref <- GLoc(spatPat=pat, alpha=theta, epsilon=Gepsilon, r=r.grid.W, contribution=FALSE, retMatrices=TRUE)
                    Gsec <- GLoc(spatPat=pat, alpha=secDir, epsilon=Gepsilon, r=r.grid.W, contribution=FALSE,
                                 distance=Gpref$distance, angles=Gpref$angles)
                    Gvec.W <- Gpref$Gloc-Gsec
                    
                    Gpref2 <- GLoc(spatPat=pat, alpha=theta, epsilon=Gepsilon, r=r.grid.L, contribution=FALSE,
                                   distance=Gpref$distance, angles=Gpref$angles)
                    Gsec <- GLoc(spatPat=pat, alpha=secDir, epsilon=Gepsilon, r=r.grid.L, contribution=FALSE,
                                 distance=Gpref$distance, angles=Gpref$angles)
                    Gvec.L <- Gpref2-Gsec
                    
                    Kpref <- cylK(spatPat=pat, alpha=theta, aspect=Kaspect, r=r.grid.W, contribution=FALSE)
                    Ksec <- cylK(spatPat=pat, alpha=secDir, aspect=Kaspect, r=r.grid.W, contribution=FALSE)
                    Kvec.W <- Kpref-Ksec
                    
                    Kpref <- cylK(spatPat=pat, alpha=theta, aspect=Kaspect, r=r.grid.L, contribution=FALSE)
                    Ksec <- cylK(spatPat=pat, alpha=secDir, aspect=Kaspect, r=r.grid.L, contribution=FALSE)
                    Kvec.L <- Kpref-Ksec
                    
                    Theta <- periodogram(spatPat=pat, p=Tp, contribution=FALSE)
                    Tvec <- thetaSpectrum(Theta, theta.grid, Tepsilon)
                    
                    vTtile <- matrix(NA, nrow= length(theta.grid), ncol= M)
                    
                    vGtile.W <- vGtile.L <- vKtile.W <- vKtile.L <- matrix(NA, nrow= length(r.grid.W), ncol= M)
                    
                    for(m in 1:M){
                      npoint=0
                      while(npoint<2 || is.null(npoint)){
                        if(est.param[3]>1){
                          recPat <- rThomas(parentIntensity = est.param[1], meanClusterSize = est.param[3],
                                             sigma = est.param[2], xdim = xdim, ydim = ydim, xScalingFactor=1, theta=theta)
                        } else {
                          recPat <- rHomPois(intensity = est.param[1]*est.param[3], xdim = xdim, ydim = ydim)
                        }
                        
                        npoint <- nrow(recPat[[1]])
                      }
                      
                      prefG <- GLoc(spatPat=recPat, alpha=theta, epsilon=Gepsilon, r=r.grid.W, retMatrices = TRUE)
                      secG <- GLoc(spatPat=recPat, alpha=secDir, epsilon=Gepsilon, r=r.grid.W,
                                   distance=prefG$distance, angles=prefG$angles)
                      vGtile.W[,m] <- prefG$Gloc - secG
                      prefG2 <- GLoc(spatPat=recPat, alpha=theta, epsilon=Gepsilon, r=r.grid.L, ,
                                     distance=prefG$distance, angles=prefG$angles)
                      secG <- GLoc(spatPat=recPat, alpha=secDir, epsilon=Gepsilon, r=r.grid.L,
                                   distance=prefG$distance, angles=prefG$angles)
                      vGtile.L[,m] <- prefG2 - secG
                      #K
                      prefK <- cylK(spatPat=recPat, alpha=theta, aspect=Kaspect, r=r.grid.W)
                      secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid.W)
                      vKtile.W[,m] <- prefK-secK
                      prefK <- cylK(spatPat=recPat, alpha=theta, aspect=Kaspect, r=r.grid.L)
                      secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid.L)
                      vKtile.L[,m] <- prefK-secK
                      #Theta spectrum
                      vTtile[,m] <- thetaSpectrum(prdgrm = periodogram(spatPat=recPat, p=Tp),
                                                  thetas = theta.grid, epsilon = Tepsilon)
                    }
                    test.GW <- testAniso(Gvec.W, vGtile.W, "G")
                    test.GL <- testAniso(Gvec.L, vGtile.L, "G")
                    test.KW <- testAniso(Kvec.W, vKtile.W, "K")
                    test.KL <- testAniso(Kvec.L, vKtile.L, "K")
                    test.T <- testAniso(Tvec, vTtile, "T")

                    print(l)
                    c(test.GW$p.val, test.GW$Tstat, test.GW$TstatRep, test.GL$p.val, test.GL$Tstat, test.GL$TstatRep,
                      test.KW$p.val, test.KW$Tstat, test.KW$TstatRep, test.KL$p.val, test.KL$Tstat, test.KL$TstatRep,
                      test.T$p.val, test.T$Tstat, test.T$TstatRep)
                  }

write.csv(output, "PLCP_THOMAS_1_BIG_long_short.csv")