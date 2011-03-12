\name{it.reb}
\alias{it.reb}
\title{
  Response Based Unit Segmentation (REBUS) Algorithm - Iterative Steps
}
\description{
  REBUS-PLS is an iterative algorithm for performing response based clustering in a PLS-PM framework. 
  This function allows to perform the iterative steps of the REBUS-PLS Algorithm (Steps 5-9). 
  It provides summarized results for final local models and the final partition of the units. 
  \cr
  Before running this function, it is necessary to run the \code{\link{res.clus}} function to choose the number of classes 
  to take into account.
}
\usage{
  it.reb(pls, hclus.res, nk, Y = NULL, stop.crit = 0.005, iter.max = 100)
}
\arguments{
  \item{pls}{ Object of class \code{"plspm"} returned by \code{\link{plspm}}. }
  \item{hclus.res}{ Object of class \code{"hclust"} returned by \code{\link{res.clus}} }
  \item{nk}{ Number of classes to take into account. This value should be defined according to the dendrogram 
        obtained by performing \code{\link{res.clus}}.}
  \item{Y}{ Optional dataset (matrix or data frame) used when argument \code{dataset=NULL} inside \code{pls}. }
  \item{stop.crit}{ Number indicating the stop criterion for the iterative algorithm. The author suggests using the threshold
       of less than 0.05\% of units changing class from one iteration to the other as stopping rule.}
  \item{iter.max}{ An integer indicating the maximum number of iterations. Default value = 100.}
}
\details{
  This function allows us to perform Steps 5-9 and the iteration procedure of the REBUS-PLS Algorithm.   
  Moreover, this function provides the final class membership for each unit and summary results 
  for the final local models. \cr
  A threshold of 6 units per class is fixed. If there is a class with less than 6 units, the algorithm stops.
  For more details refer to \code{\link{rebus.pls}}.
  
  \emph{For expert users: if you want to test REBUS-PLS on several number of classes you need to run this
  function several times by changing} \code{nk}. \cr

  When the object \code{pls} does not contain a data matrix (i.e. \code{pls$data=NULL}), the user must provide the data matrix or data frame in \code{Y}.
}
\value{
  An object of class \code{"rebus"}, basically a list with the following elements:
  \item{loadings}{Matrix of standardized loadings (i.e. correlations with LVs.) for each local model.}
  \item{path.coefs}{Matrix of path coefficients for each local model.}
  \item{quality}{Matrix containing the average communalities, the average redundancies, the R2 values,
       and the GoF index for each local model.}
  \item{segments}{Vector defining the class membership of each unit.} 
  \item{origdata.clas}{ The numeric matrix with original data and with a new column defining class 
        membership of each unit.}  
}
\references{ 
  Esposito Vinzi, V., Trinchera, L., Squillacciotti, S., and Tenenhaus, M. (2008)
  REBUS-PLS: A Response-Based Procedure for detecting Unit Segments in PLS Path Modeling.
  \emph{Applied Stochastic Models in Business and Industry (ASMBI)}, \bold{24}, pp. 439-458. 
  
  Trinchera, L. (2007) Unobserved Heterogeneity in Structural Equation Models: a new approach to 
  latent class detection in PLS Path Modeling. \emph{Ph.D. Thesis}, University of Naples "Federico II", 
  Naples, Italy.

  \url{http://www.fedoa.unina.it/view/people/Trinchera,_Laura.html}
}
\author{Laura Trinchera, Gaston Sanchez}
\seealso{\code{\link{plspm}},  \code{\link{rebus.pls}},  \code{\link{res.clus}}}
\examples{
  \dontrun{
  ## example of rebus analysis with simulated data
  data(sim.data)
  ## First compute GLOBAL model
  sim.inner <- matrix(c(0,0,0,0,0,0,1,1,0),3,3,byrow=TRUE)
  dimnames(sim.inner) <- list(c("Price","Quality","Satisfaction"),
                            c("Price","Quality","Satisfaction"))
  sim.outer <- list(c(1,2,3,4,5),c(6,7,8,9,10),c(11,12,13)) 
  sim.mod <- c("A","A","A")  ## reflective indicators
  sim.global <- plspm(sim.data, inner=sim.inner, 
                      outer=sim.outer, modes=sim.mod)
  sim.global
  ## Then compute cluster analysis on residuals of global model
  sim.res.clus <- res.clus(sim.global)
  ## To conclude run iteration algorithm
  rebus.sim <- it.reb(sim.global, sim.res.clus, nk=2, 
                     stop.crit=0.005, iter.max=100)
  ## You can also compute complete outputs 
  ## for local models by running:
  local.rebus <- local.models(sim.global, rebus.sim)
  }
}
