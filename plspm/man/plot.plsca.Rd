\name{plot.plsca}
\alias{plot.plsca}
\title{Plot PLSCA basic results}
\description{
  Plot method for objects of class \code{"plsca"}
}
\usage{
  \method{plot}{plsca}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{"plsca"} returned by \code{\link{plsca}}.}
  \item{\dots}{Further arguments are ignored.}
}
\details{
  The function \code{plot.plsca} displays four plots: \cr
  1) the circles of correlations of the X-scores (T-components) and the Y-scores (U-components) \cr
  2) the plot of T-components (t1,t2) and the plot of U-components (u1,u1) \cr
  3) barplots of the explained variances and communalities \cr
  4) the dot charts with (t1,u1) and (t2,u2) \cr

}
\author{Gaston Sanchez}
\seealso{\code{\link{plsca}}}
\examples{
  \dontrun{
  ## example of PLSCA with the vehicles dataset
  data(vehicles)
  can <- plsca(vehicles[,1:12], vehicles[,13:16])
  plot(can)
  }
}
