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

  # INFERENCE FOR THOMAS PROCESS

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



  # PARAMETERS
  theta = pi/6
  xScalingFactor = 0.8
  var=3
  scale= 0.03
  mu=log(400)-var/2
  model="exponential"
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



  registerDoParallel(cores=16)




  set.seed(789)
  seeds <- sample(1e6, 1e3)
  output <- foreach(l=1:N, .combine=cbind, .export = c("simWind", "angleMatrix","xScalingFactor", "rLGCPtransformed", "rLGCP", "Lest", 
                                                      "testAniso","cylK", "GLoc", "periodogram", "thetaSpectrum", "model",
                                                      "owin", "xdim", "ydim", "var", "scale", "mu", "lgcp.estpcf", "rAnisotrop",
                                                      "theta", "grid_length", "M", "Gepsilon", "Kaspect", "Tepsilon", "Tp", "theta.grid",
                                                      "as.ppp", "r.grid.W", "rHomPois", "rThomas", "infer.LGCP"),
                    .multicombine = TRUE, .maxcombine = N) %dopar% {
                      set.seed(seeds[l])
                      pat <- rAnisotrop(theta=theta, xScalingFactor, process=rLGCPtransformed, 
                                        xdim=xdim, ydim=ydim, model=model, mu=mu, var=var, scale=scale)
                      L <- Lest(as.ppp(pat[[1]], W=owin(xdim, ydim)))
                      max.R <- min(diff(xdim)/4, L$r[1+which(((cumsum(abs(L$theo-L$trans)))/L$r)[-1]==max(((cumsum(abs(L$theo-L$trans)))/L$r)[-1]))])
                      r.grid.L <- seq(0, max.R, length=grid_length+1)[-1]
                      secDir <- (theta+pi/2)%%pi
                      
                      est.param <- infer.LGCP(pat, model)
                      
                      
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
                      
                      vGtile.W <- vGtile.L <- vKtile.W <- vKtile.L <- matrix(NA, nrow= length(r.grid.W), ncol= M)
                      vTtile <- matrix(NA, nrow= length(theta.grid), ncol= M)
                      
                      for(m in 1:M){
                        npoint=0
                        if(l==314) set.seed(3274)
                        while(npoint<2 || is.null(npoint)){
                          recPat <- rLGCPtransformed(var = est.param[1], scale = est.param[2], model=model,
                                                    mu = est.param[3], xdim = xdim, ydim = ydim)
                          
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


write.csv(output, "LGCP_LGCP_08_BIG.csv")