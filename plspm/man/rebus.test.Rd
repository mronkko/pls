\name{rebus.test}
\alias{rebus.test}
\title{Permutation Test for REBUS Multi-Group Comparison}
\description{
  Performs permutation tests for comparing pairs of groups from a REBUS object.
}
\usage{
  rebus.test(pls, reb, Y = NULL)
}
\arguments{
  \item{pls}{ Object of class \code{"plspm"} returned by \code{\link{plspm}}. }
  \item{reb}{ Object of class \code{"rebus"} returned by either \code{\link{rebus.pls}} or \code{\link{it.reb}}. }
  \item{Y}{ Optional dataset (matrix or data frame) used when argument \code{dataset=NULL} inside \code{pls}. }
}
\details{
  A permutation test on path coefficients, loadings, and GoF index is applied to the classes obtained from REBUS, by comparing two classes at a time. That is to say, a permutation test is applied on pair of classes. The number of permutations in each test is 100. \cr
  The number of classes handled by \code{rebus.test} is limited to 6. 

  When the object \code{pls} does not contain a data matrix (i.e. \code{pls$data=NULL}), the user must provide the data matrix or data frame in \code{Y}.
}
\value{
  An object of class \code{"rebus.test"}, basically a list containing the results of each pair of compared classes.
  In turn, each element of the list is also a list with the results for the path coefficients, loadings, and GoF index.
}
\references{
  Chin, W.W. (2003) A permutation procedure for multi-group comparison of PLS models. In: Vilares, M., Tenenhaus, M., Coelho, P., Esposito Vinzi, V., Morineau, A. (Eds.) \emph{PLS and Related Methods - Proceedings of the International Symposium PLS03.} Decisia, pp. 33-43. 
}
\author{Laura Trinchera, Gaston Sanchez}
\seealso{ \code{\link{rebus.pls}}, \code{\link{local.models}} }
\examples{
  \dontrun{
  ## example of rebus analysis with simulated data
  data(sim.data)
  ## Calcualte plspm
  sim.inner <- matrix(c(0,0,0,0,0,0,1,1,0),3,3,byrow=TRUE)
  dimnames(sim.inner) <- list(c("Price","Quality","Satisfaction"),
                            c("Price","Quality","Satisfaction"))
  sim.outer <- list(c(1,2,3,4,5),c(6,7,8,9,10),c(11,12,13)) 
  sim.mod <- c("A","A","A")  ## reflective indicators
  sim.global <- plspm(sim.data, inner=sim.inner, 
                      outer=sim.outer, modes=sim.mod)
  sim.global
  ## Cluster analysis on residuals of global model
  sim.res.clus <- res.clus(sim.global)
  ## Iterative steps of REBUS algorithm on 2 classes
  rebus.sim <- it.reb(sim.global, sim.res.clus, nk=2, 
                     stop.crit=0.005, iter.max=100)
  sim.permu <- rebus.test(sim.global, rebus.sim)
  sim.permu
  sim.permu$test_1_2
  ## or equivalently
  sim.permu[[1]]
  }
}