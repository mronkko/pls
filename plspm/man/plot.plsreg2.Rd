\name{plot.plsreg2}
\alias{plot.plsreg2}
\title{Plot PLSR2 basic results}
\description{
  Plot method for objects of class \code{"plsreg2"}
}
\usage{
  \method{plot}{plsreg2}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{"plsreg2"} returned by \code{\link{plsreg2}}.}
  \item{\dots}{Further arguments are ignored.}
}
\details{
  The function \code{plot.plsreg2} displays four plots: \cr
  1) the circle of correlations between the variables and the first two components (scores) \cr
  2) the plot of weights and loadings (w,c) \cr
  3) the plot of components (t1,t2) and the plot of components (u1,u2) \cr
  4) the plot of components (t1,u1) and the plot of components (t2,u2) \cr
}
\author{Gaston Sanchez}
\seealso{\code{\link{plsreg2}},  \code{\link{plsreg1}}}
\examples{
  \dontrun{
  ## example of plsreg2
  data(vehicles)
  pls2 <- plsreg2(vehicles[,1:12], vehicles[,13:16])
  plot(pls2)
  }
}

