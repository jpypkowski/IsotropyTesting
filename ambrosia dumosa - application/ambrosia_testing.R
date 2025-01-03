# Ambrosia Dumosa

library(spatstat)

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


# load and clean data
data <- read.table("doi_10.5063_AA_connolly.206.1-DATA.data", header = TRUE)
uq <- unique(data[,1])
 
for(i in uq){  # remove repeated observations
  ixs <- which(data[,1]==i)
  if(length(ixs) > 1) data <- data[-ixs[-1],]
}
window <- simWind(xdim=c(-50,50), ydim=c(-50,50)) #re-centred to align with the PLCP inference code
Ambrosia <-list(pattern=cbind(as.numeric(data[,2])-50, as.numeric(data[,3])-50), winBor=window[[2]],
                lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff)





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
  p.val <- (1+sum(as.numeric(Tstat) < TstatRep, na.rm=TRUE))/(1+sum(!is.na(TstatRep)))
  return(list(p.val=p.val, Tstat=Tstat, TstatRep=TstatRep))
}

########
# LGCP #
########


estLGCP <- infer.LGCP(Ambrosia, covariance="exponential") # estimating parameters and assessing dependency range
rd <- -estLGCP[2]*log(0.1)
rmax <- 1.7
r.grid.W <- seq(0.05, 1.7, by=0.05)
Kaspect = 0.15
theta <- 0
secDir <- pi/2


prefK <- cylK(spatPat=Ambrosia, alpha=theta, aspect=Kaspect, r=r.grid.W)
secK <- cylK(spatPat=Ambrosia, alpha=secDir, aspect=Kaspect, r=r.grid.W)
vecK_LGCP <- prefK-secK

vKtileW_LGCP <- read.csv("vKtile_LGCP.csv")[,-1]
test_LGCP <- testAniso(OGv=vecK_LGCP, Vs=vKtileW_LGCP, DSS="K")$p.val 


########
# PLCP #
########


sigma_sq_sample <- read.csv("sigma_sq_sample.csv")[,-1]
sigma_sqs <- sigma_sq_sample[6501:12500]
quantile(sqrt(sigma_sqs), c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975,1)) # assessing the dependency range


vKtileW_PLCP <- read.csv("vKtile_PLCP.csv")[seq(2,200, by=2),-1]
dim(vKtileW_PLCP)

r.grid.W <- seq(0.25, 25, by=0.25)
prefK <- cylK(spatPat=Ambrosia, alpha=theta, aspect=Kaspect, r=r.grid.W)
secK <- cylK(spatPat=Ambrosia, alpha=secDir, aspect=Kaspect, r=r.grid.W)
vecK_PLCP_tile <- prefK-secK

test_PLCP_25 <- testAniso(OGv=vecK_PLCP_tile, Vs=vKtileW_PLCP, DSS="K")$p.val
test_PLCP_19 <- testAniso(OGv=vecK_PLCP_tile[1:76], Vs=vKtileW_PLCP[1:76,], DSS="K")$p.val
test_PLCP_13 <- testAniso(OGv=vecK_PLCP_tile[1:52], Vs=vKtileW_PLCP[1:52,], DSS="K")$p.val


##########
# TILING #
##########


vKtile_4 <- read.csv("vKtile_k4_025.csv")[,-1]
vKtile_5 <- read.csv("vKtile_k5_025.csv")[,-1]
vKtile_6 <- read.csv("vKtile_k6_025.csv")[,-1]
vKtile_7 <- read.csv("vKtile_k7_025.csv")[,-1]
vKtile_8 <- read.csv("vKtile_k8_025.csv")[,-1]


test_k4 <- testAniso(OGv=vecK_LGCP_tile, Vs=vKtile_4, DSS="K")$p.val
test_k5 <- testAniso(OGv=vecK_LGCP_tile, Vs=vKtile_5, DSS="K")$p.val
test_k6 <- testAniso(OGv=vecK_LGCP_tile, Vs=vKtile_6, DSS="K")$p.val
test_k7 <- testAniso(OGv=vecK_LGCP_tile, Vs=vKtile_7, DSS="K")$p.val
test_k8 <- testAniso(OGv=vecK_LGCP_tile, Vs=vKtile_8, DSS="K")$p.val


