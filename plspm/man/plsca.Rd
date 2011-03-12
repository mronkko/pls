\name{plsca}
\alias{plsca}
\title{PLS-CA: Partial Least Squares Canonical Analysis}
\description{
  Performs partial least squares canonical analysis for two blocks of data. Compared to PLSR2, the blocks of variables in PLS-CA play a symmetric role (i.e. there is neither predictors nor predictands)
}
\usage{
  plsca(X, Y, nc = NULL, scaled = TRUE)
}
\arguments{
  \item{X}{A numeric matrix or data frame (X-block).}
  \item{Y}{A numeric matrix or data frame (Y-block).}
  \item{nc}{The number of extracted PLS components (\code{NULL} by default)}
  \item{scaled}{A logical value indicating whether scaling data should be performed (\code{TRUE} by default).}
}
\details{
  Arguments \code{X} and \code{Y} must contain more than one variable.    

  No missing data are allowed.

  When \code{nc=NULL} the number of components is determined by taking the minimum between the number of columns from \code{X} and \code{Y}.

  When \code{scaled=TRUE} the data is scaled to standardized values (mean=0, variance=1). Otherwise the data will only be centered (mean=0).
}
\value{
  An object of class \code{"plsca"}, basically a list with the following elements:
  \item{x.scores}{scores of the X-block (also referred to as T components).}
  \item{x.wgs}{weights of the X-block.}
  \item{x.loads}{loadings of the X-block.}
  \item{y.scores}{scores of the Y-block (also referred to as U components).}
  \item{y.wgs}{weights of the Y-block.}
  \item{y.loads}{loadings of the Y-block.}
  \item{cor.xt}{correlations between X and T.}
  \item{cor.yu}{correlations between Y and U.}
  \item{cor.tu}{correlations between T and U.}
  \item{cor.xu}{correlations between X and U.}
  \item{cor.yt}{correlations between Y and T.}
  \item{R2X}{explained variance of X by T.}
  \item{R2Y}{explained variance of Y by U.}
  \item{com.xu}{communality of X with U.}
  \item{com.yt}{communality of Y with T.}
}
\references{ 
  Tenenhaus, M. (1998) \emph{La Regression PLS. Theorie et Pratique}. Editions TECHNIP, Paris.
}
\author{Gaston Sanchez}
\seealso{\code{\link{print.plsca}},\code{\link{plot.plsca}}}
\examples{
  \dontrun{
  ## example of PLSCA with the vehicles dataset
  data(vehicles)
  can <- plsca(vehicles[,1:12], vehicles[,13:16])
  can
  plot(can)
  }
}
