\name{plot.plsreg1}
\alias{plot.plsreg1}
\title{Plot PLSR1 basic results}
\description{
  Plot method for objects of class \code{"plsreg1"}
}
\usage{
  \method{plot}{plsreg1}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{"plsreg1"} returned by \code{\link{plsreg1}}.}
  \item{\dots}{Further arguments are ignored.}
}
\details{
  The function \code{plot.plsreg1} displays four plots: \cr
  1) the circle of correlations between the variables and the first two components (scores) \cr
  2) the plot of components (t1,t2) and the plot of components (t1,u1) \cr
  3) the plot of T2 Hotelling Confidence ellipse \cr
  4) the dot chart comparing the response variable (Y) with the predicted values (Y-predicted) \cr
}
\author{Gaston Sanchez}
\seealso{\code{\link{plsreg1}},  \code{\link{plsreg2}}}
\examples{
  \dontrun{
  ## example of plsreg1 
  data(vehicles)
  pls1 <- plsreg1(vehicles[,1:12], vehicles[,13])
  plot(pls1)
  }
}
