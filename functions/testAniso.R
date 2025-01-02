# ANISOTROPY TESTING PROCEDURE GIVEN ORIGINAL AND REPLICATED VALUES FOR COMPUTATION OF TEST STATISTICS
# OGv - the vector v of the directional summary statistic used to construct a test statistic
#       computed from the observed pattern
# Vs  - vectors of the DSS obtained from pattern replicates (in a matrix form)
# DSS - a name of the DSS used; use "G", "K", or "T" (for Theta spectrum)


testAniso <- function(OGv, Vs, DSS){
  if(!(DSS%in%c("G","K","T","P"))) stop("DSS must be either 'G', 'K', 'T', or 'P'.")
  reps <- ncol(Vs)
  mHat <- apply(Vs, 1, mean, na.rm=TRUE)
  vector <- OGv-mHat
  TstatRep <- rep(NA, reps)
  if(!(DSS=="G")){
    CHatDiag <- diag(1/apply(Vs, 1, var, na.rm=TRUE))
    CHatDiag[CHatDiag==Inf] <- 0 # if the variance is 0, a contribution at the corresponding distance 
                                 # will not be include to not overwhelm the test statistic
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