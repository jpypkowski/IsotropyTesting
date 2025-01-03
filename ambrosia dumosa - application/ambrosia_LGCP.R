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

rLGCPtransformed <- function(model, mu, xdim, ydim, ...){
  pat <- rLGCP(model=model, mu=mu, win=owin(xrange=xdim, yrange=ydim), ...)
  window <-simWind(xdim, ydim)
  return(list(pattern = cbind(pat$x, pat$y), winBor=window[[2]],
              lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff))
}

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


data <- read.table("doi_10.5063_AA_connolly.206.1-DATA.data", header = TRUE)
uq <- unique(data[,1])
for(i in uq){
  ixs <- which(data[,1]==i)
  if(length(ixs) > 1) data <- data[-ixs[-1],]
}
window <- simWind(xdim=c(-50,50), ydim=c(-50,50))
Ambrosia <-list(pattern=cbind(as.numeric(data[,2])-50, as.numeric(data[,3])-50), winBor=window[[2]],
                lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff)

model="exponential"
estLGCP <- infer.LGCP(Ambrosia, model)
rd <- -estLGCP[2]*log(0.1)
rmax <- 1.7
r.grid.W <- seq(0.05, 1.7, by=0.05)
M <- 1000
theta <- 0
secDir <- pi/2
Kaspect <- 0.15


set.seed(7841)
seeds <- sample(1e6, 1e3)


registerDoParallel(cores=64)
vKtile <- foreach(l=1:M, .combine=cbind, .export = c("simWind", "rLGCPtransformed", "cylK",
                                           "rLGCP", "owin", "xdim", "ydim", "model", "estLGCP", 
                                           "theta", "M", "Kaspect", "secDir", "r.grid.W",
                                           "as.ppp", "seeds"), 
        .multicombine = TRUE, .maxcombine = M) %dopar% {
            set.seed(seeds[l])
            recPat <- rLGCPtransformed(model, mu=estLGCP[3], xdim=c(-50, 50), ydim=c(-50, 50), var=estLGCP[1], scale=estLGCP[2])
            prefK <- cylK(spatPat=recPat, alpha=theta, aspect=Kaspect, r=r.grid.W)
            secK <- cylK(spatPat=recPat, alpha=secDir, aspect=Kaspect, r=r.grid.W)
            print(l)
            prefK-secK
        }

write.csv(vKtile, "vKtile_LGCP.csv")