\name{rebus.pls}
\alias{rebus.pls}
\title{ Response Based Unit Segmentation (REBUS) Algorithm - complete - }
\description{
  REBUS-PLS is an iterative algorithm for performing response based clustering in a 
  PLS-PM framework. 
  \cr
  This function allows to perform all the steps of the REBUS-PLS Algorithm. 
}
\usage{
  rebus.pls(pls, Y = NULL, stop.crit = 0.005, iter.max = 100)
}
\arguments{
  \item{pls}{ Object of class \code{"plspm"} returned by \code{\link{plspm}}. }
  \item{Y}{ Optional dataset (matrix or data frame) used when argument \code{dataset=NULL} inside \code{pls}. }
  \item{stop.crit}{ Number indicating the stop criterion for the iterative algorithm. The author suggests using the threshold
       of less than 0.05\% of units changing class from one iteration to the other as stopping rule.}
  \item{iter.max}{ An integer indicating the maximum number of iterations. Default value = 100.}
}
\details{
  REBUS-PLS is an iterative algorithm that starting from the global model allows us to detect local 
  models performing better than the global model. The algorithm is composed by 9 steps:
  \cr
  Step 1: The first step of the REBUS-PLS algorithm involves computing the global model on all the observed 
  units, by performing a simple PLS Path Modeling analysis (\code{\link{plspm}}). 
  \cr
  Steps 2 and 3: In the second step, the communality and the structural residuals of each unit from the global 
  model are obtained. The number of classes (K) to be taken into account during the successive
  iterations and the initial composition of the classes are obtained by performing a Hierarchical 
  Cluster Analysis on the computed residuals (both from the measurement and the structural models).
  \cr
  Step 4: Choose the number of classes (K) to take into account according to the dendrogram obtained
  at Step 3.
  \cr
  Step 5: Once the number of classes to consider is chosen, the initial composition of the classes are 
  obtained according to results of the Hierarchical Cluster Analysis.
  \cr
  Step 6: A PLS Path Modeling analysis is performed on each formed class and K provisional local 
  models are estimated. 
  \cr
  The group-specific parameters computed at the previous step are used to 
  compute the communality and the structural residuals of each unit from each local model. 
  \cr
  Step 7: The Closeness Measure (CM) of each unit from each local model is computed. 
  \cr
  Step 8: Each unit is, therefore, assigned to the closest local model, i.e. to the model from which 
  each unit shows the smallest CM value. 
  \cr
  Step 9: Once the composition of the classes is updated, K new local models are estimated.
  
  Step 5-9 have to be iterate until there is convergence on class composition. 
  The author suggests using the threshold of less than 0.05\% of units changing class from one 
  iteration to the other as stopping rule. \cr
  In addition, a threshold of 6 units per class is fixed. If there is a class with less than 6 units, the algorithm stops.

  This function allows us to perform in a unique function all the steps of the REBUS-PLS Algorithm. 
  \cr
  In particular, firstly it computes the communality and structural residuals from the global model 
  and performs a Hierarchical Cluster Analysis on this residuals (by applying \code{\link{res.clus}} 
  function). Then it shows the obtained dendrogram and asks the user to choose the number of clusters 
  to be take into account.This value need to be an integer value higher than 1.
  Once the number of classes is chosen, this function performs the iteration steps of 
  the REBUS-PLS Algorithm and directly provides the final class membership for each unit, 
  summary results for the final local models, and the Group Quality Index (GQI) (Trinchera, 2007)
  ((by applying \code{\link{it.reb}} function). 
  
  Once the convergence is achieved you should:
  \cr
  (i) Estimate the final local models
  \cr
  (ii) Validate the Group Quality Index
  
  Both these last two steps can be performed by running \code{\link{local.models}}. 
  
  N.B. Once the global model is computed it is possible to obtain the same results as the ones obtained
  by running \code{rebus.pls} function, by running the function \code{\link{res.clus}} and then the 
  function \code{\link{it.reb}}. 
  \cr
  If you want to test REBUS-PLS on several number of classes you need to run \code{\link{it.reb}}
  function several times by changing nk. This allows you to avoid computing the residuals each time from the
  global model.
}
\value{
  An object of class \code{"rebus"}, basically a list with the following elements:
  \item{loadings}{Matrix of standardized loadings (i.e. correlations with LVs.) for each local model.}
  \item{path.coefs}{Matrix of path coefficients for each local model.}
  \item{quality}{Matrix containing the average communalities, the average redundancies, the R2 values,
       and the GoF values for each local model.}
  \item{segments}{Vector defining for each unit the class membership.} 
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
\seealso{\code{\link{plspm}}, \code{\link{res.clus}}, \code{\link{resclus.plot}}, \code{\link{it.reb}}, \code{\link{rebus.test}}, \code{\link{local.models}}}
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
  ## run rebus.pls function on the sim.data and choose the number of classes
  ## to be taken into account according to the displayed dendrogram.
  rebus.sim <- rebus.pls(sim.global, stop.crit = 0.005, iter.max = 100)
  ## You can also compute complete outputs 
  ## for local models by running:
  local.rebus <- local.models(sim.global, rebus.sim)
  ##
  ## ONCE COMPUTED THE GLOBAL MODEL IT IS POSSIBLE TO OBTAIN THE SAME RESULTS AS THE ONES
  ## OBTAINED BY RUNNING rebus.pls FUNCTION BY RUNNING THE FUNCTION res.clus
  ## AND THEN THE FUNCTION it.reb.
  ##
  ## Example:
  ## Once the global model is calculated you need to perform cluster analysis on residuals of the global model:
  sim.res.clus <- res.clus(sim.global)
  ## and then run the iteration algorithm:
  rebus.sim.twofunc <- it.reb(sim.global, sim.res.clus, nk=2, 
                     stop.crit=0.005, iter.max=100)
  ##
  ## rebus.sim = rebus.sim.twofunc
  }
}
