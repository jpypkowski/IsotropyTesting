# PLOTTING A SPATIAL PATTERN ON A RECTANGULAR OBSERVATION WINDOW
plotSpatPat <- function(spatPat, axes=FALSE, 
                        xlab="", ylab="", ...){
  plot(spatPat$pattern, frame.plot = FALSE, axes=axes,
       xlim=spatPat$lim[,1], ylim=spatPat$lim[,2], 
       xlab=xlab, ylab=ylab, ...)
  polygon(x=spatPat$winBor[,1], y=spatPat$winBor[,2], 
          lty=1, lwd=1, border="darkgrey")
}