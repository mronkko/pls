\name{res.clus}
\alias{res.clus}
\title{Cluster Analysis on communality and structural residuals of a PLS-PM object}
\description{
  This function computes communality and structural residuals from the global model 
  and performs a Hierarchical Cluster Analysis on these residuals according to the REBUS algorithm.
}
\usage{
  res.clus(pls, Y = NULL)
}
\arguments{
  \item{pls}{ Object of class \code{"plspm"} returned by \code{\link{plspm}} .}
  \item{Y}{ Optional dataset (matrix or data frame) used when argument \code{dataset=NULL} inside \code{pls}. }
}
\details{
  The function \code{res.clus} comprises the second and third steps of the REBUS-PLS Algorithm. 
  In particular it computes communality and structural residuals (as defined in Trinchera (2007) and
  Esposito Vinzi \emph{et al.} (2008)) of each unit from the global model (step two of REBUS-PLS Algorithm).
  Then it performs a Hierarchical Cluster Analysis on these residuals (step three of REBUS-PLS Algorithm).
  As a result, this function directly provides a dendrogram obtained from a Hierarchical Cluster Analysis. 

  The number of classes (K) to be taken into account during the successive steps of the REBUS-PLS 
  Algorithm (performed by \code{\link{it.reb}}), and the initial composition of the classes are obtained according
  to the results of the Hierarchical Cluster Analysis. Users must choose K according to this dendrogram and use it as an argument in the 
  \code{\link{it.reb}} function.

  When the object \code{pls} does not contain a data matrix (i.e. \code{pls$data=NULL}), the user must provide the data matrix or data frame in \code{Y}.
}
\value{
  An Object of class \code{"hclust"} containing the results of the Hierarchical Cluster Analysis 
  on the communality and structural residuals.
}
\references{ 
  Esposito Vinzi, V., Trinchera, L., Squillacciotti, S., and Tenenhaus, M. (2008)
  REBUS-PLS: A Response-Based Procedure for detecting Unit Segments in PLS Path Modeling.
  \emph{Applied Stochastic Models in Business and Industry (ASMBI)}, \bold{24}, pp. 439-458. 
  
  Trinchera, L. (2007) Unobserved Heterogeneity in Structural Equation Models: a new approach to 
  latent class detection in PLS Path Modeling. \emph{Ph.D. Thesis}, University of Naples "Federico II", 
  Naples, Italy.
}
\author{Laura Trinchera, Gaston Sanchez}
\seealso{\code{\link{resclus.plot}}, \code{\link{it.reb}}, \code{\link{plspm}}}
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
  ## Then compute cluster analysis on the residuals of global model
  sim.res.clus <- res.clus(sim.global)
  }
}