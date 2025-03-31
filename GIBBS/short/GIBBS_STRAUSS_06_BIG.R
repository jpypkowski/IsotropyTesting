library(spatstat)
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


#STRAUSS SIMULATION USING METROPOLIS-HASTINGS

rmhStrauss_rep <- function(beta, gamma, R, xdim, ydim){
  pat <- rmh(list(cif="strauss",par=list(beta=beta,gamma=gamma,r=R),
                  w=c(xdim, ydim)), nrep=1e6, verbose=FALSE)
  window <-simWind(xdim, ydim)
  return(list(pattern=cbind(pat$x, pat$y), winBor=window$winBor, 
              lim=window$lim,  xDiff=window$xDiff, yDiff=window$yDiff))
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



# Inference procedure for Strauss process incl. assessing interaction range

infer.Strauss <- function(pattern){
  estimates <- rep(NA, 3)
  xrange <- c(pattern[[2]][1,1], pattern[[2]][3,1])
  yrange <- c(pattern[[2]][2,2], pattern[[2]][3,2])
  L <- Lest(as.ppp(pattern$pattern, owin(xrange=xrange, yrange=yrange)))
  ix <- which(L$theo - L$trans == max (L$theo - L$trans))
  estimates[3] <- L$r[ix]
  estimates[1:2] <- exp(ppm(as.ppp(pattern$pattern, owin(xrange=xrange, yrange=yrange)) ~ 1, Strauss(estimates[3]))$coef)
  return(estimates)
}





# PARAMETERS
theta = pi/6
sigma=0.015
xdim = ydim = c(-0.5, 0.5)

N <- 1000
M <- 1000

Gepsilon <- pi/4
Kaspect <- 0.15


r.grid.W <- seq(0, 0.05, length=37)[-1]



registerDoParallel(cores=8)



loaded <- read.csv("Gibbs06_big_2.csv")[,-1]
ixs <- which(loaded[,1]=="pattern")

set.seed(9876)
seeds <- sample(1e6, 1e3)
output <- foreach(l=1:N, .combine=cbind, .export = c("simWind", "angleMatrix","Lest", 
                                                     "testAniso","cylK", "GLoc", 
                                                     "owin", "xdim", "ydim", "Strauss", "ppm", "rmhStrauss_rep", "rmh",
                                                     "theta", "M", "Gepsilon", "Kaspect",
                                                     "as.ppp", "r.grid.W", "rHomPois", "infer.Strauss"),
                  .multicombine = TRUE, .maxcombine = N) %dopar% {
                    set.seed(seeds[l])
                    window <- simWind(xdim, ydim)
                    if(l<N){
                      pat <- list(pattern = loaded[(ixs[l]+1):(ixs[l+1]-1),], winBor=window[[2]],
                                  lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff)    
                    } else {
                      pat <- list(pattern = loaded[(ixs[l]+1):nrow(loaded),], winBor=window[[2]],
                                  lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff)  
                    }
                    pat[[1]] <- cbind(as.numeric(pat[[1]][,1]), as.numeric(pat[[1]][,2]))

                  
                    secDir <- (theta+pi/2)%%pi
                    
                    est.param <- infer.Strauss(pat)
                    
                    
                    Gpref <- GLoc(spatPat=pat, alpha=theta, epsilon=Gepsilon, r=r.grid.W, contribution=FALSE, retMatrices=TRUE)
                    Gsec <- GLoc(spatPat=pat, alpha=secDir, epsilon=Gepsilon, r=r.grid.W, contribution=FALSE,
                                 distance=Gpref$distance, angles=Gpref$angles)
                    Gvec.W <- Gpref$Gloc-Gsec
                    
                   
                    
                    Kpref <- cylK(spatPat=pat, alpha=theta, aspect=Kaspect, r=r.grid.W, contribution=FALSE)
                    Ksec <- cylK(spatPat=pat, alpha=secDir, aspect=Kaspect, r=r.grid.W, contribution=FALSE)
                    Kvec.W <- Kpref-Ksec
                    
                   
                    
                    
                    
                    vGtile.W <- vKtile.W <- matrix(NA, nrow= length(r.grid.W), ncol= M)
                    
                    for(m in 1:M){
                      npoint=0
                      while(npoint<2 || is.null(npoint)){
                        
                          if(est.param[2] < 1){
                          recPat <- rmhStrauss_rep(beta = est.param[1], gamma = est.param[2],
                                                 R = est.param[3], xdim = xdim, ydim = ydim)
                        } else {
                          recPat <- rHomPois(intensity=nrow(pat[[1]])/(diff(xdim)*diff(ydim)), xdim=xdim, ydim=ydim)
                        }
                        
                        
                        
                        npoint <- nrow(recPat[[1]])
                      }
                      
                      prefG <- GLoc(spatPat=recPat, alpha=theta, epsilon=Gepsilon, r=r.grid.W, retMatrices = TRUE)
                      secG <- GLoc(spatPat=recPat, alpha=secDir, epsilon=Gepsilon, r=r.grid.W,
                                   distance=prefG$distance, angles=prefG$angles)
                      vGtile.W[,m] <- prefG$Gloc - secG
                      
                      #K
                      prefK <- cylK(spatPat=recPat, alpha=theta, aspect=Kaspect, r=r.grid.W)
                      secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid.W)
                      vKtile.W[,m] <- prefK-secK
                     
                 
                    }
                    test.GW <- testAniso(Gvec.W, vGtile.W, "G")
                    test.KW <- testAniso(Kvec.W, vKtile.W, "K")
                    print(l)
                    c(test.GW$p.val, test.GW$Tstat, test.GW$TstatRep, 
                      test.KW$p.val, test.KW$Tstat, test.KW$TstatRep)
                  }

write.csv(output, "GIBBS_STRAUSS_06_BIG_short.csv")