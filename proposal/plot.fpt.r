#' Plot FPT
#'
#' A modified version of \code{\link[adehabitatLT]{plot.fipate}} from \code{adehabitatLT} that plots the fpt for multiple radii. 
#'
#' @param x a \code{fipati} (first passage time) class object
#' @param radii a vector of radii for x

plot.fpt <- function (x, radii, main = "", ...) 
{
  if (!inherits(x, "fipati")) 
    stop("x should be of class 'fipati'")
  
  for(i in 1:length(radii))
  {
    u <- x[[1]]
    plot(attr(u, "date"), u[, i], 
         ylab = "FPT", main=paste(main, radii[i], sep=""), ...)
    lines(attr(u, "date"), u[, i])
  }
}
