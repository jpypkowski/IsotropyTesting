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




sigma_sq_sample <- read.csv("sigma_sq_sample.csv")[(7501:12500),-1]
alpha_sample <- read.csv("alpha_sample.csv")[(7501:12500),-1]
rho_L_sample <- read.csv("rho_L_sample.csv")[(7501:12500),-1]


xdim <- ydim <- c(-50, 50)
r.grid.W <- seq(0.125, 25, by=0.125)
M <- 1000
theta <- 0
secDir <- pi/2
Kaspect <- 0.15


set.seed(68)
seeds <- sample(1e6, 1e3)


registerDoParallel(cores=64)
vKtile <- foreach(l=1:M, .combine=cbind, .export = c("simWind", "rDirLine", "cylK", "rvonmises",
                                                     "xdim", "ydim", "sigma_sq_sample", 
                                                     "theta", "M", "Kaspect", "secDir", "r.grid.W",
                                                     "seeds", "alpha_sample", "rho_L_sample"), 
                  .multicombine = TRUE, .maxcombine = M) %dopar% {
                    set.seed(seeds[l])
                    npoint <- 0
                    while(npoint<2){
                      index <- sample(1:5000, 1)
                      recPat <- rDirLine(lineIntensity = rho_L_sample[index], alongIntensity = alpha_sample[index],
                                         sigma = sqrt(sigma_sq_sample[index]), xdim = xdim, ydim = ydim, kappa=KofA(1), theta=pi/6)
                      npoint <- nrow(recPat[[1]])
                    }                    
                    prefK <- cylK(spatPat=recPat, alpha=theta, aspect=Kaspect, r=r.grid.W)
                    secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid.W)
                    print(l)
                    prefK-secK
                  }

write.csv(vKtile, "vKtile_PLCP.csv")


