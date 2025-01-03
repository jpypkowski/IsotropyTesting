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
              alpha_sample=alpha_sample, rho_L_sample=rho_L_sample, accept=accept, lines=lines))
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



lines <- matrix(data=c(0.5, cos(11*pi/6), sin(11*pi/6), 11*pi/6,
                       -0.5, cos(2*pi/11), sin(2*pi/11), 2*pi/11,
                       0.5, cos(pi*0.8), sin(pi*0.8), pi*0.8,
                       0.0, cos(pi*0.7), sin(pi*0.7), pi*0.7,
                       -5, cos(pi/90), sin(pi/90), pi/90,
                       0.0, cos(pi*0.3), sin(pi*0.3), pi*0.3,
                       -0.2, cos(pi*0.85), sin(pi*0.85), pi*0.85,
                       0.1, cos(pi*0.45), sin(pi*0.45), pi*0.45),ncol=4, byrow=TRUE)

a=55

lines <- cbind((lines[,1])*100, lines[,-1])
set.seed(563)
while(nrow(lines)<40){
  u_angle <- runif(1, 0, 2*pi)
  
  u_prop <- c(cos(u_angle), sin(u_angle))
  
  if((u_angle > 0 && u_angle <= pi/2) || (u_angle >= pi && u_angle <= 3*pi/2)){
    y_prop <- runif(1, -a/tan(u_angle) - a, a/tan(u_angle)+a)
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
    if(y_prop > 0){
      psi <- u_angle %% pi
      x <- y_prop
    } else {
      psi <- (pi-u_angle) %% pi 
      x <- -y_prop
    }
    if((x-a)*tan(psi)<=a) {
      continue <- TRUE
    } else{
      continue <- FALSE
    }
  } 
  if(continue == TRUE){
    lines <- rbind(lines, c(y_prop, u_prop, u_angle))
  }
}


set.seed(356)
infer_PLCP_ambrosia <- infer.PLCP(Ambrosia, a1=1+1e-6, b1=1e-12, a2=1+1e-6, b2=1e-12, M=12500, M_integral=33^2, 
           prop_sigma=3e-1, max_sigma_sq=30, W_margin=5,
           sigma_sq_init=10, lines=lines, birth_prob=1/3, death_prob=1/3)


write.csv(infer_PLCP_ambrosia$sigma_sq_sample, "sigma_sq_sample.csv")
write.csv(infer_PLCP_ambrosia$alpha_sample, "alpha_sample.csv")
write.csv(infer_PLCP_ambrosia$rho_L_sample, "rho_L_sample.csv")
