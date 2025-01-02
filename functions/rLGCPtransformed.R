# Adjusting spatstat's simulation of LGCP to our data format
# model      - variance model, e.g. "exponential"
# mu         - mean of the underlying Gaussian field
# xdim, ydim - observation window bounds in x and y dimensions
# ...        - further arguments to be passed to spatsat's rLGCP;
#              in particular 'scale' and 'var' for exponential variance

rLGCPtransformed <- function(model, mu, xdim, ydim, ...){
  pat <- rLGCP(model=model, mu=mu, win=owin(xrange=xdim, yrange=ydim), ...)
  window <-simWind(xdim, ydim)
  return(list(pattern = cbind(pat$x, pat$y), winBor=window[[2]],
              lim=window[[3]], xDiff=window$xDiff, yDiff=window$yDiff))
}