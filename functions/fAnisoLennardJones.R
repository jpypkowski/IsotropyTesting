# Energy function with anisotropic 6-12 Lennard-Jones energy function
# spatPat  - a point pattern specified using simWind
# theta    - angle $\theta$
# alpha    - parameter alpha of the isotropic component
# epsilon1 - interaction strength $\rho_1$ (NOT $\epsilon$) 
# epsilon2 - interaction strength $\rho_2$ (NOT $\epsilon$) 
# sigma1   - scale parameter $\sigmna_1$ 
# sigma2   - scale parameter $\sigmna_2$ 
# log      - set to TRUE to compute a logarithm of energy rather than energy

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
