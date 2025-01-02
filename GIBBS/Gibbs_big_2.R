# This file simulates patterns from Gibbs process with different levels of anisotropy



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


BirthDeathMove <- function(interactionF, p, q, M, xdim=c(0,1), ydim=c(0,1), addMargin = 0.05, ..., initn=NULL, initIntensity=NULL, initPat=NULL){
  
  # interactionF - pairwise interaction function
  # p - probability of birth at a birth step
  # q - probability of a Metropolis-Hastings move step
  # M - number of iterations of the algorithm
  # xdim, ydim - observation window limits
  # ... - parameters of the pairwise interaction function
  # initn - number of points in the step-1 pattern
  # initIntensity - intensity for the step-1 Poisson process
  xdim_old <- xdim
  ydim_old <- ydim
  xdim <- xdim+c(-addMargin, addMargin)
  ydim <- ydim+c(-addMargin, addMargin)
  
  if(is.null(initn) & is.null(initIntensity) & is.null(initPat)) stop("Define initn or initIntensity")
  if(!is.null(initPat)){
    curPat <- initPat
  } else  if(is.null(initn)) {
    curPat <- rHomPois(intensity = initIntensity, xdim=xdim, ydim=ydim)
  } else {
    curPat <- rHomPois(numPoints = initn, xdim=xdim, ydim=ydim)
  }
  densities <- rep(NA, M+1)
  densities[1] <- interactionF(spatPat=curPat, log=TRUE, ...)
  
  for(m in 1:M){
    Rm3 <- runif(1)
    n <- nrow(curPat$pattern)
    if (Rm3 <= q && n>0) {
      I <- sample(1:n, 1)
      Rm <- runif(1)
      ksi <- c(runif(1, xdim[1], xdim[2]), runif(1, ydim[1], ydim[2]))
      propPat <- curPat
      propPat$pattern[I,] <- ksi
      r <- exp(interactionF(spatPat=propPat, log=TRUE, ...)-interactionF(spatPat=curPat, log=TRUE, ...))
      if(Rm <= r){
        curPat <- propPat
      } 
    } else {
      Rm1 <- runif(1)
      Rm2 <- runif(1)
      ########### ADAPT
      propPat <- curPat
      if(Rm1 <= p){
        ksi <- c(runif(1, xdim[1], xdim[2]), runif(1, ydim[1], ydim[2]))
        propPat$pattern <- rbind(propPat$pattern, ksi)
        r <- exp(interactionF(spatPat=propPat, log=TRUE, ...)-interactionF(spatPat=curPat, log=TRUE, ...))
        if(Rm2 <= r){
          curPat <- propPat
        } 
      } else {
        if (n>0) {
          eta <- sample(1:n, 1)
          if(n==2){
            propPat$pattern <- matrix(propPat$pattern[-eta,], nrow=1)
          } else if(n==1){
            propPat$pattern <- propPat$pattern[-eta,]
          } else
            propPat$pattern <- propPat$pattern[-eta,]
          r <- exp(interactionF(spatPat=propPat, log=TRUE, ...)-interactionF(spatPat=curPat, log=TRUE,...))
          if(Rm2 <= r){
            curPat <- propPat
          } 
        } 
        
      }
      
    }
    densities[m+1] <- interactionF(spatPat=curPat, log=TRUE, ...)
    # print(c(m,n))
    
  }
  window <- simWind(xdim_old, ydim_old)
  finalPat <- list(pattern=curPat$pattern[(curPat$pattern[,1]>=xdim_old[1] & curPat$pattern[,1]<=xdim_old[2] & curPat$pattern[,2]>=ydim_old[1] & curPat$pattern[,2]<=ydim_old[2]),],
                   winBor=window[[2]], lim=window[[3]],  xDiff=window[[4]], yDiff=window[[5]])
  return(list(spatPat = finalPat, densities = densities))
}

# BDM_BIG

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

theta = pi/6
xScalingFactor = c(1, 0.8, 0.6, 0.4)
beta = 0.5
alpha = -log(beta)
epsilon = 10
sigma = 0.02

xdim=c(-0.5, 0.5)
ydim=c(-0.5, 0.5)

p = 1/3
q = 1/3
sqrt.split = 3
addMargin_small <- addMargin_big  <- 0.05

M_small = 800
M_big= 400 

registerDoParallel(cores=8)



set.seed(7283)
seed.Gibbs <- sample(1e6, 1e3)
SF <- xScalingFactor[1]

Gibbs_big <- foreach(l=1:1000, .combine = rbind, 
                     .export = c("fAnisoLennardJones", "BirthDeathMove", "angleMatrix",
                                 "alpha", "epsilon", "sigma", "seed.Gibbs", "theta", "xScalingFactor",
                                 "xdim", "ydim", "p", "q", "M", "addMargin", "initn"), 
                     .multicombine = TRUE, .maxcombine = 10000) %dopar% {
                       set.seed(seed.Gibbs[l])
                       
                       
                       init <- runif(9, 5, 75)
                       pat <- BDM_BIG(sqrt.split=sqrt.split, interactionF=fAnisoLennardJones, p=p, q=q, M_small=M_small, M_big=M_big, 
                               xdim=c(-0.5, 0.5), ydim=c(-0.5, 0.5),
                               addMargin_small = addMargin_small, addMargin_big = addMargin_big, alpha=-log(beta), 
                               sigma1=sigma*sqrt(2-SF^(1/3)), sigma2=sigma*SF^(1/6),
                               # epsilon1=epsilon/sqrt(SF), epsilon2=epsilon*sqrt(SF),
                               # epsilon1=epsilon*sqrt(SF), epsilon2=epsilon/sqrt(SF),
                               epsilon1=epsilon, epsilon2=epsilon,
                               initn_small=init, theta=theta)[[1]]

                       print(nrow(pat))
                       rbind(c("pattern", "pattern"), pat)
                               

                     }
write.csv(Gibbs_big, "Gibbs1_big_2.csv")


set.seed(1121)
seed.Gibbs <- sample(1e6, 1e3)
SF <- xScalingFactor[2]

Gibbs_big <- foreach(l=1:1000, .combine = rbind, 
                     .export = c("fAnisoLennardJones", "BirthDeathMove", "angleMatrix",
                                 "alpha", "epsilon", "sigma", "seed.Gibbs", "theta", "xScalingFactor",
                                 "xdim", "ydim", "p", "q", "M", "addMargin", "initn"), 
                     .multicombine = TRUE, .maxcombine = 10000) %dopar% {
                       set.seed(seed.Gibbs[l])
                       
                       
                       init <- runif(9, 5, 75)
                       pat <- BDM_BIG(sqrt.split=sqrt.split, interactionF=fAnisoLennardJones, p=p, q=q, M_small=M_small, M_big=M_big, 
                               xdim=c(-0.5, 0.5), ydim=c(-0.5, 0.5),
                               addMargin_small = addMargin_small, addMargin_big = addMargin_big, alpha=-log(beta), 
                               sigma1=sigma*sqrt(2-SF^(1/3)), sigma2=sigma*SF^(1/6),
                               # epsilon1=epsilon/sqrt(SF), epsilon2=epsilon*sqrt(SF),
                               # epsilon1=epsilon*sqrt(SF), epsilon2=epsilon/sqrt(SF),
                               epsilon1=epsilon, epsilon2=epsilon,
                               initn_small=init, theta=theta)[[1]]

                       print(nrow(pat))
                       rbind(c("pattern", "pattern"), pat)
                               

                     }
write.csv(Gibbs_big, "Gibbs08_big_2.csv")


set.seed(1445)
seed.Gibbs <- sample(1e6, 1e3)
SF <- xScalingFactor[3]

Gibbs_big <- foreach(l=1:1000, .combine = rbind, 
                     .export = c("fAnisoLennardJones", "BirthDeathMove", "angleMatrix",
                                 "alpha", "epsilon", "sigma", "seed.Gibbs", "theta", "xScalingFactor",
                                 "xdim", "ydim", "p", "q", "M", "addMargin", "initn"), 
                     .multicombine = TRUE, .maxcombine = 10000) %dopar% {
                       set.seed(seed.Gibbs[l])
                       
                       
                       init <- runif(9, 5, 75)
                       pat <- BDM_BIG(sqrt.split=sqrt.split, interactionF=fAnisoLennardJones, p=p, q=q, M_small=M_small, M_big=M_big, 
                               xdim=c(-0.5, 0.5), ydim=c(-0.5, 0.5),
                               addMargin_small = addMargin_small, addMargin_big = addMargin_big, alpha=-log(beta), 
                               sigma1=sigma*sqrt(2-SF^(1/3)), sigma2=sigma*SF^(1/6),
                               # epsilon1=epsilon/sqrt(SF), epsilon2=epsilon*sqrt(SF),
                               # epsilon1=epsilon*sqrt(SF), epsilon2=epsilon/sqrt(SF),
                               epsilon1=epsilon, epsilon2=epsilon,
                               initn_small=init, theta=theta)[[1]]

                       print(nrow(pat))
                       rbind(c("pattern", "pattern"), pat)
                               

                     }
write.csv(Gibbs_big, "Gibbs06_big_2.csv")


set.seed(854)
seed.Gibbs <- sample(1e6, 1e3)
SF <- xScalingFactor[4]

Gibbs_big <- foreach(l=1:1000, .combine = rbind, 
                     .export = c("fAnisoLennardJones", "BirthDeathMove", "angleMatrix",
                                 "alpha", "epsilon", "sigma", "seed.Gibbs", "theta", "xScalingFactor",
                                 "xdim", "ydim", "p", "q", "M", "addMargin", "initn"), 
                     .multicombine = TRUE, .maxcombine = 10000) %dopar% {
                       set.seed(seed.Gibbs[l])
                       
                       
                       init <- runif(9, 5, 75)
                       pat <- BDM_BIG(sqrt.split=sqrt.split, interactionF=fAnisoLennardJones, p=p, q=q, M_small=M_small, M_big=M_big, 
                               xdim=c(-0.5, 0.5), ydim=c(-0.5, 0.5),
                               addMargin_small = addMargin_small, addMargin_big = addMargin_big, alpha=-log(beta), 
                               sigma1=sigma*sqrt(2-SF^(1/3)), sigma2=sigma*SF^(1/6),
                               # epsilon1=epsilon/sqrt(SF), epsilon2=epsilon*sqrt(SF),
                               # epsilon1=epsilon*sqrt(SF), epsilon2=epsilon/sqrt(SF),
                               epsilon1=epsilon, epsilon2=epsilon,
                               initn_small=init, theta=theta)[[1]]

                       print(nrow(pat))
                       rbind(c("pattern", "pattern"), pat)
                               

                     }
write.csv(Gibbs_big, "Gibbs04_big_2.csv")


set.seed(16578)
seed.Gibbs <- sample(1e6, 1e3)
SF <- xScalingFactor[1]

Gibbs_big <- foreach(l=1:1000, .combine = rbind, 
                     .export = c("fAnisoLennardJones", "BirthDeathMove", "angleMatrix",
                                 "alpha", "epsilon", "sigma", "seed.Gibbs", "theta", "xScalingFactor",
                                 "xdim", "ydim", "p", "q", "M", "addMargin", "initn"), 
                     .multicombine = TRUE, .maxcombine = 10000) %dopar% {
                       set.seed(seed.Gibbs[l])
                       
                       
                       init <- runif(9, 5, 75)
                       pat <- BDM_BIG(sqrt.split=sqrt.split, interactionF=fAnisoLennardJones, p=p, q=q, M_small=M_small, M_big=M_big, 
                               xdim=c(-0.5, 0.5), ydim=c(-0.5, 0.5),
                               addMargin_small = addMargin_small, addMargin_big = addMargin_big, alpha=-log(beta), 
                               sigma1=sigma*sqrt(2-SF^(1/3)), sigma2=sigma*SF^(1/6),
                               # epsilon1=epsilon/sqrt(SF), epsilon2=epsilon*sqrt(SF),
                               # epsilon1=epsilon*sqrt(SF), epsilon2=epsilon/sqrt(SF),
                               epsilon1=epsilon, epsilon2=epsilon,
                               initn_small=init, theta=theta)[[1]]

                       print(nrow(pat))
                       rbind(c("pattern", "pattern"), pat)
                               

                     }
write.csv(Gibbs_big, "Gibbs_oracle_big_2.csv")

