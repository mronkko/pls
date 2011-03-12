\name{plot.nipals}
\alias{plot.nipals}
\title{Plot NIPALS basic results}
\description{
  Plot method for objects of class \code{"nipals"}
}
\usage{
  \method{plot}{nipals}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{"nipals"} returned by \code{\link{nipals}}.}
  \item{\dots}{Further arguments are ignored.}
}
\details{
  The function \code{plot.nipals} displays three plots: \cr
  1) screeplots of the pseudo eigenvalues \cr
  2) the circle of correlations between variables and the first two components (scores) \cr
  3) the plot of the two first scores, and the plot of the two first loadings 
}
\author{
  Gaston Sanchez
}
\seealso{\code{\link{nipals}}}
\examples{
  \dontrun{
  ## example of NIPALS algorithm
  data(wines)
  nip1 <- nipals(wines[,-1], nc=5)
  plot(nip1)
  }
}
