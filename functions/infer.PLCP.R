# BAYESIAN INFERENCE FOR PLCP + ALL AUXILIARY FUNCTIONS
# AS PER HYBRID MARKOV CHAIN MONTE CARLO ALGORITHM (Møller, Safavimanesh, & Rasmussen 2016)
# arXiv:1503.07423v5
# THE MAIN FUNCTION IS LOCATED ON THE BOTTOM OF THIS FILE precede by all auxiliary functions


# component of $R_{\sigma^2}$: $\int_W f\{p_{u_j^\perp}(x_i-y_i\|\sigma^2}dx as a vector with entries indexed by j
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
    vector <- sum(dnorm(proj_dist, mean = 0, sd = sqrt(sigma_sq))) # in this case 'vector' is a scalar!
  }
  
  return(vector)
}

# entries for computation of the product at the end of $R_{\sigma^2}$ (see 'prod_sum_ratio' in the main function)
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
    func_vals <- dnorm(proj_dist, mean = 0, sd = sqrt(sigma_sq)) # in this case 'func_vals' is a vector, not a matrix!
  }
  
  return(func_vals)
}

# acceptance probability of a proposed line birth (Eq. 30)
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
  print(sum_log)
  if(log){
    return(list(logRb = log_1st_bit + log_2nd_bit + sum_log, int = int_prop_l, dens = dens_prep_l))
  } else {
    return(list(Rb = exp(log_1st_bit + log_2nd_bit  + sum_log), int = int_prop_l, dens = dens_prep_l))
  }
}

# acceptance probability of a proposed line death (Eq. 31)
Rdeath <- function(rm, l_rm, lambda_J, rho_L, alpha, dens_mat, int_vals, point_area, k, log=FALSE){
  print(rm)
  print(dens_mat)
  log_1st_bit <- log(rho_L * lambda_J * abs(l_rm[3]) / (k + 1))
  print(log_1st_bit)
  log_2nd_bit <- -alpha*int_vals[rm]*point_area
  print(log_2nd_bit)
  if(sum(dens_mat[,-rm]) > 0){
    if(ncol(dens_mat)==2){
      sum_log <- sum(log(1 + dens_mat[,rm]/sum(dens_mat[,-rm])), na.rm=TRUE)
      print(sum_log)
    } else {
      sum_log <-sum(log(1 + dens_mat[,rm]/apply(dens_mat[,-rm], MARGIN = 1, sum)), na.rm=TRUE)
    }
    if(log){
      print(1)
      return(- log_1st_bit - log_2nd_bit - sum_log)
    } else {
      print(2)
      return(1/exp(log_1st_bit + log_2nd_bit  + sum_log))
    }
    
  } else{
    if(log){
      print(3)
      return(-Inf)
    } else {
      print(4)
      return(0)
    }
  }
}


# THE MAIN FUNCITON
# spatPat        - point pattern specified using simWind
# a1, b1, a2, b2 - parameters of priors for parameter $\nu$ (our notation) 
#                  or $\alpha$ (notation by Møller et al., 2016; see above )
# M              - number of iterations of the algorithm
# M_integral     - number of points used to approximate integrals over the observation window
# prop_sigma     - 1/2 of maximum spread of the uniform distribution used to propose $\sigma^2$
# max_sigma_sq   - a parameter used to limit $\sigma^2$ for faster convergence
# W_margin       - a distance by which the observation window is extended by for the simulation of latent lines
# sigma_sq_init  - initial value of $\sigma^2$
# lines          - set of initial lines (see files titled as "PLCP_PLCP_..." for an example)
# birth_prob     - probability of a birth proposal at birth-death-move step
# death_prob     - probability of a death proposal at birth-death-move step

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
    
    b1_sample[iter+1] <- b1_sample[1] + sum(ints_cur) * point_area
    a2_sample[iter+1] <- a2_sample[1] + nrow(lines)
    
    
    
    alpha_sample[iter] <- rgamma(1, shape=a1_post, rate=b1_sample[iter+1])
    rho_L_sample[iter] <- rgamma(1, shape=a2_sample[iter+1], rate=b2_post)
    
    
    # Metropolis-Hastings updating
    # Updating sigma_sq
    
    # proposal ~ Uniform distribution
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
    
    
    logR <- log_p_ratio + log_int_ratio + prod_sum_ratio # logarithm of $R_{\sigma^2}$
    
    threshold <- runif(1, 0, 1)
    if(log(threshold) < logR) { # accepting or rejecting the proposal for $\sigma^2$
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
        if((x-a)*tan(psi)<=a) { # check if the proposal is within requirements
          continue <- TRUE
        } else {
          continue <- FALSE
        }
      }
      
      

      if(continue){
        l_prop <- c(y_prop, u_prop, u_angle)
        logRb <- Rbirth(l_prop, points=spatPat$pattern, lines, lambda_J, rho_L = rho_L_sample[iter], 
                        alpha = alpha_sample[iter], sigma_sq = sigma_sq_sample[iter+1], MC_points, dens_mat = dens_cur, point_area=point_area,
                        log=TRUE)
        
        threshold <- runif(1,0,1) # accepting or rejecitng the proposed line birth
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
      if(log(threshold) < logRd) { # accepting or rejecting a line death
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
                          log=TRUE) # numerator as probability of accepting a new line
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
        # denominator is the reciprocal od the probability of death of the line proposed to move
        log_denom <- Rdeath(rm, l_rm, lambda_J, rho_L = rho_L_sample[iter], alpha = alpha_sample[iter], 
                            dens_mat = dens_cur, int_vals = ints_cur, point_area=point_area, k = nrow(lines), log=TRUE)
        
        
        logRm <- log_num[[1]]+log_denom #multiplying acceptance probabilites of birth and death
        print(c(log_num[[1]], log_denom, logRm ))
        if(log(threshold) < logRm) {
          lines <- rbind(l_less, l_prop)
          ints_cur <- c(ints_cur[-rm], log_num$int)#### ADD THE NEW ONE AS IN BIRTH
          dens_cur <- cbind(dens_cur[,-rm], log_num$dens)
        }
      }
      
    }
    
    # INCLUDE THE BELOW IF YOU WANT TO TRACK THE SAMPLED VALUES WHILE THE ALGORITHM RUNS
    
    # if(iter %% 500 == 0){
    #   par(mfrow=c(2,2))
    #   plot(sigma_sq_sample[1:(iter+1)], type="l")
    #   plot(alpha_sample[1:(iter+1)], type="l")
    #   plot(rho_L_sample[1:(iter+1)], type="l")
    # }
    
    
  } # END LOOP
  return(list(sigma_sq_sample=sigma_sq_sample, b1_sample=b1_sample, a2_sample=a2_sample,
              alpha_sample=alpha_sample, rho_L_sample=rho_L_sample, accept=accept))
}