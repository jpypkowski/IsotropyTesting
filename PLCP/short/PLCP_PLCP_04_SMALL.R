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
    # print(sum(is.na(dens_prop)))
    # print(c(log_p_ratio, log_int_ratio, prod_sum_ratio, log(threshold)))
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
    
    
    
  } # END LOOP
  return(list(sigma_sq_sample=sigma_sq_sample, b1_sample=b1_sample, a2_sample=a2_sample,
              alpha_sample=alpha_sample, rho_L_sample=rho_L_sample, acccept=accept))
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
xScalingFactor = 0.4
lineIntensity=16
alongIntensity= 25
sigma=0.015
xdim = ydim = c(-0.5, 0.5)/2

N <- M <- 1000

Gepsilon <- pi/4
Kaspect <- 0.15


iterB <- 10000

r.grid.W <- seq(0, 0.05, length=37)[-1]



registerDoParallel(cores=8)
                       




set.seed(16)
seeds <- sample(1e6, 1e3)

output <- foreach(l=1:N, .combine=cbind, .export = c("simWind", "angleMatrix", "xScalingFactor", "rDirLine", "KofA", "Lest",
                                                     "testAniso","cylK", "GLoc", "rvonmises",
                                                     "owin", "xdim", "ydim", "lineIntensity", "alongIntensity", "sigma",
                                                     "theta", "M", "Gepsilon", "Kaspect", 
                                                     "as.ppp", "r.grid.W", "iterB", "int_sum", "sum_log_matrix",
                                                     "Rbirth", "Rdeath", "infer.PLCP"),
                  .multicombine = TRUE, .maxcombine = N) %dopar% {
                    set.seed(seeds[l])
                    pat <- rDirLine(lineIntensity=lineIntensity, alongIntensity=alongIntensity, theta=theta, xdim=xdim, ydim=ydim,
                                    sigma=sigma, kappa=KofA(xScalingFactor))
                    
                    secDir <- (theta+pi/2)%%pi
                    lines <- matrix(data=c(0.6, cos(pi*0.3), sin(pi*0.3), pi*0.3,
                          0.5, cos(11*pi/6), sin(11*pi/6), 11*pi/6,
                          -0.5, cos(2*pi/11), sin(2*pi/11), 2*pi/11,
                          1.5, cos(pi/6), sin(pi/6), pi/6,
                          -0.5, cos(pi/6), sin(pi/6), pi/6,
                          0.5, cos(pi*0.8), sin(pi*0.8), pi*0.8,
                          -2.4, cos(pi*0.9), sin(pi*0.9), pi*0.9,
                          -5, cos(pi/90), sin(pi/90), pi/90
                        ),ncol=4, byrow=TRUE)
                    lines[,1] <- lines[,1]/2
                    params_inf <- infer.PLCP(spatPat=pat, a1=1+1e-6, b1=1e-12, a2=1+1e-6, b2=1e-12, M=iterB, M_integral=33^2,
                                             prop_sigma=1e-5, max_sigma_sq=0.3, W_margin=0.05,
                                             sigma_sq_init=0.0225^2, lines=lines, birth_prob=1/3, death_prob=1/3)
                    Gpref <- GLoc(spatPat=pat, alpha=theta, epsilon=Gepsilon, r=r.grid.W, contribution=FALSE, retMatrices=TRUE)
                    Gsec <- GLoc(spatPat=pat, alpha=secDir, epsilon=Gepsilon, r=r.grid.W, contribution=FALSE)
                    Gvec.W <- Gpref$Gloc-Gsec
                    
                    Kpref <- cylK(spatPat=pat, alpha=theta, aspect=Kaspect, r=r.grid.W, contribution=FALSE)
                    Ksec <- cylK(spatPat=pat, alpha=secDir, aspect=Kaspect, r=r.grid.W, contribution=FALSE)
                    Kvec.W <- Kpref-Ksec
                    
                    vGtile.W <- vKtile.W <-  matrix(NA, nrow= length(r.grid.W), ncol= M)
                    for(m in 1:M){
                      npoint=0
                      while(npoint<2){
                        index <- sample((iterB/2+1):iterB, 1)
                        recPat <- rDirLine(lineIntensity = params_inf$rho_L[index], alongIntensity = params_inf$alpha[index],
                                           sigma = sqrt(params_inf$sigma_sq_sample[index]), xdim = xdim, ydim = ydim, kappa=KofA(1), theta=pi/6)
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

write.csv(output, "PLCP_PLCP_04_SMALL_short.csv")
