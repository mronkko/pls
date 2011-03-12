\name{nipals}
\alias{nipals}
\title{NIPALS: Non-linear Iterative Partial Least Squares}
\description{
  Principal Component Analysis with NIPALS algorithm
}
\usage{
  nipals(x, nc = 2, scaled = TRUE)
}
\arguments{
  \item{x}{A numeric matrix or data frame.}
  \item{nc}{Number of components kept in the results (by default 2)}
  \item{scaled}{A logical value indicating whether scaling data is performed (\code{TRUE} by default).}
}
\details{
  The function \code{nipals} performs Principal Component Analysis of a data matrix that may contain missing data.
}
\value{
  An object of class \code{"nipals"}, basically a list with the following elements:
  \item{values}{The pseudo eigenvalues.}
  \item{scores}{The extracted scores.}
  \item{loadings}{The loadings.}
  \item{cor.sco}{Correlations between the variables and the scores.}
  \item{disto}{Squared distance of the observations to the origin.}
  \item{contrib}{Contributions of the observations (rows).}
  \item{cos}{Squared cosinus.}
  \item{dmod}{Distance to the Model.}
  When the analyzed data contain missing values, the help interpretation tools (e.g. \code{cor.sco, disto, contrib, cos, dmod}) may not be meaningful, that is to say, some of the results may not be coherent. 
}
\references{
  Tenenhaus, M. (1998) \emph{La Regression PLS. Theorie et Pratique.} Editions TECHNIP, Paris.

  Tenenhaus, M. (2007) \emph{Statistique. Methodes pour decrire, expliquer et prevoir}. Dunod, Paris.
}
\author{Gaston Sanchez}
\seealso{\code{\link{plot.nipals}}}
\examples{
  \dontrun{
  ## example of NIPALS algorithm
  data(wines)
  nip1 <- nipals(wines[,-1], nc=5)
  plot(nip1)

  ## USArrests data vary 
  nip2 <- nipals(USArrests)
  plot(nip2)
  }
}
