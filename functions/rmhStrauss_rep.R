# STRAUSS SIMULATION USING METROPOLIS-HASTINGS
# beta, gamma, R - parameters as obtained from infer. Strauss (in the same order)
# xdim, ydim     - observation window bounds

rmhStrauss_rep <- function(beta, gamma, R, xdim, ydim){
  pat <- rmh(list(cif="strauss",par=list(beta=beta,gamma=gamma,r=R),
                  w=c(xdim, ydim)), nrep=1e6, verbose=FALSE)
  window <-simWind(xdim, ydim)
  return(list(pattern=cbind(pat$x, pat$y), winBor=window$winBor, 
              lim=window$lim,  xDiff=window$xDiff, yDiff=window$yDiff))
}