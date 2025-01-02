# THETA SPECTRUM
# prdgrm  - periodogram computed from the pattern
# thetas  - angles for which the function is computed
# epsilon - bandwidth $h$


thetaSpectrum <- function(prdgrm, thetas, epsilon){ 
  spectrum <- rep(NA, length(thetas))
  for(i in 1:length(thetas)){
    consider <- (abs((prdgrm$angles-thetas[i])%%pi) < epsilon)
    spectrum[i] <- sum(prdgrm$sdf[consider],
                       na.rm=TRUE)/ sum(rowSums(consider, na.rm=TRUE))
  }
  return(spectrum)                 
}
