\name{plspm}
\alias{plspm}
\alias{plspm.fit}
\title{PLS-PM: Partial Least Squares Path Modeling}
\description{
  Estimate path models with latent variables by partial least squares approach
}
\usage{
  plspm(x, inner, outer, modes = NULL, scheme = "centroid", scaled = TRUE, boot.val = FALSE,
        br = NULL, plsr = FALSE, tol = 0.00001, iter = 100, dataset = TRUE)

  plspm.fit(x, inner, outer, modes = NULL, scheme = "centroid", scaled = TRUE, tol = 0.00001, iter = 100)
}
\arguments{
  \item{x}{A numeric matrix or data frame containing the manifest variables.}
  \item{inner}{A square (lower triangular) boolean matrix representing the inner model (i.e. the path relationships betwenn latent variables).}
  \item{outer}{List of vectors with column indices from x representing the outer model (i.e. which manifest variables correspond to the latent variables).}
  \item{modes}{A character vector indicating the type of measurement for each latent variable.
               \code{"A"} for reflective measurement or \code{"B"} for formative measurement (\code{NULL} by default).}
  \item{scheme}{A string of characters indicating the type of inner weighting scheme. Possible values are \code{"centroid"}, \code{"factor"}, or \code{"path"}.}
  \item{scaled}{A logical value indicating whether scaling data is performed (\code{TRUE} by default).}
  \item{boot.val}{A logical value indicating whether bootstrap validation is performed (\code{FALSE} by default).}
  \item{br}{An integer indicating the number bootstrap resamples. Used only when \code{boot.val=TRUE}.}
  \item{plsr}{A logical value indicating whether pls regression is applied (\code{FALSE} by default).}
  \item{tol}{Decimal value indicating the tolerance crietrion for the iterations (\code{tol=0.00001} by default).}
  \item{iter}{An integer indicating the maximum number of iterations (\code{iter=100} by default).}
  \item{dataset}{A logical value indicating whether the data matrix should be retrieved (\code{TRUE} by default).}
}
\details{
  The function \code{plspm} estimates a path model by partial least squares approach while providing all results. \cr

  The function \code{plspm.fit} performs the basic PLS algorithm and provides simple results (e.g. outer weights, LVs scores, path coefficients, R2, and loadings). \cr
  
  The argument \code{inner.mat} is a matrix of zeros and ones that indicates the structural relationships between latent variables. This must be a lower triangular matrix. \code{inner.mat} will contain a 1 when column \code{j} affects row \code{i}, 0 otherwise. \cr

  The argument \code{sets} is a list of vectors of indices indicating the sets of manifest variables associated to the latent variables. The length of \code{sets} must be equal to the number of rows of \code{inner.mat}. \cr

  The argument \code{modes} is a character vector indicating the type of measurement for each latent variable. A value of \code{"A"} is used when a latent variable has reflective manifest variables, and a value of \code{"B"} is used when the latent variable has formative manifest variables. The length of \code{modes} must be equal to the number of rows of \code{inner.mat} (i.e. equal to the length of \code{sets}). \cr

  The argument \code{scaled} is \code{TRUE} by default. This means that data in \code{x} is scaled to standardized values (i.e. mean=0 and variance=1, calculating the variance dividing by \code{N} instead of \code{N-1}). Unless the data has the same scale for all variables, it is strongly recommended to keep this argument unchanged. \cr

  When bootstrap validation is performed, the default number of re-samples is 100. However, \code{br} can be specified in a range from 100 to 1000. \cr

  The argument \code{plsr} gives the option to calculate the path coefficients by means of pls regression. \cr

  The argument \code{tol} can be specified in a range from 0 to 0.001. \cr

  The minimum value of \code{iter} is 100. \cr
}
\value{
  An object of class \code{"plspm"}. When the function \code{plspm.fit} is called, it returns a list with basic results:
  \item{outer.mod}{Results of the outer (measurement) model. Includes: outer weights, standardized loadings, communalities, and redundancies.}
  \item{inner.mod}{Results of the inner (structural) model. Includes: path coefficients and R-squared for each endogenous latent variable.}
  \item{latents}{Matrix of standardized latent variables (variance=1 calculated divided by \code{N}) obtained from centered data (mean=0).}
  \item{scores}{Matrix of latent variables used to estimate the inner model. If \code{scaled=FALSE} then \code{scores} are latent variables calculated with the original data (non-stardardized). If \code{scaled=TRUE} then \code{scores} and \code{latents} have the same values.}
  \item{out.weights}{Vector of outer weights.}
  \item{loadings}{Vector of standardized loadings (i.e. correlations with LVs.)}
  \item{path.coefs}{Matrix of path coefficients (this matrix has a similar form as \code{inner.mat}).}
  \item{r.sqr}{Vector of R-squared coefficients.}
  \cr
  If the function \code{plspm} is called, the previous list of results also contains the following elements:
  \item{outer.cor}{Correlations between the latent variables and the manifest variables (also called crossloadings).}
  \item{inner.sum}{Summarized results by latent variable of the inner model. Includes: type of LV, type of measurement, number of indicators, R-squared, average communality, average redundancy, and average variance extracted}
  \item{effects}{Path effects of the structural relationships. Includes: direct, indirect, and total effects.} 
  \item{unidim}{Results for checking the unidimensionality of blocks (These results are only meaningful for reflective blocks).}
  \item{gof}{Table with indexes of Goodness-of-Fit. Includes: absolute GoF, relative GoF, outer model GoF, and inner model GoF.}
  \item{data}{Data matrix containing the manifest variables used in the model. Only when \code{dataset=TRUE}.}
  \item{boot}{List of bootstrapping results; only available when argument \code{boot.val=TRUE}.}
}
\references{
  Tenenhaus, M., Esposito Vinzi, V., Chatelin Y.M., and Lauro, C. (2005) PLS path modeling. \emph{Computational Statistics & Data Analysis}, \bold{48}, pp. 159-205. 

  Tenenhaus, M., and Pages, J. (2001) Multiple factor analysis combined with PLS path modelling. Application to the analysis of relationships between physicochemical variables, sensory profiles and hedonic judgements. \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{58}, pp. 261-273.

  Tenenhaus, M., and Hanafi, M. A bridge between PLS path modeling and multi-block data analysis. \emph{Handbook on Partial Least Squares (PLS): Concepts, methods, and applications.} Springer: In press.

  Lohmoller, J.-B. (1989) \emph{Latent variables path modelin with partial least squares.} Heidelberg: Physica-Verlag.

  Wold, H. (1985) Partial Least Squares. In: Kotz, S., Johnson, N.L. (Eds.), \emph{Encyclopedia of Statistical Sciences}, Vol. 6. Wiley, New York, pp. 581-591.

  Wold, H. (1982) Soft modeling: the basic design and some extensions. In: K.G. Joreskog & H. Wold (Eds.), \emph{Systems under indirect observations: Causality, structure, prediction}, Part 2, pp. 1-54. Amsterdam: Holland.
}
\author{
  Gaston Sanchez 

  \url{http://www.scribd.com/plspm}
}
\seealso{\code{\link{print.plspm}}, \code{\link{summary.plspm}}, \code{\link{plot.plspm}}.}
\examples{
  \dontrun{
  ## example of PLS-PM in ecological analysis
  ## model with three LVs and formative indicators
  data(arizona)
  ari.inner <- matrix(c(0,0,0,0,0,0,1,1,0),3,3,byrow=TRUE)
  dimnames(ari.inner) <- list(c("ENV","SOIL","DIV"),c("ENV","SOIL","DIV"))
  ari.outer <- list(c(1,2),c(3,4,5),c(6,7,8)) 
  ari.mod <- c("B","B","B")  ## formative indicators
  res1 <- plspm(arizona, inner=ari.inner, outer=ari.outer, modes=ari.mod,
  scheme="factor", scaled=TRUE, plsr=TRUE)
  res1
  summary(res1)
  
  ## typical example of PLS-PM in customer satisfaction analysis
  ## model with six LVs and reflective indicators
  data(satisfaction)
  IMAG <- c(0,0,0,0,0,0)
  EXPE <- c(1,0,0,0,0,0)
  QUAL <- c(0,1,0,0,0,0)
  VAL  <- c(0,1,1,0,0,0)
  SAT  <- c(1,1,1,1,0,0) 
  LOY  <- c(1,0,0,0,1,0)
  sat.inner <- rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
  sat.outer <- list(1:5,6:10,11:15,16:19,20:23,24:27)
  sat.mod <- rep("A",6)   ## reflective indicators
  res2 <- plspm(satisfaction, sat.inner, sat.outer, sat.mod, scaled=FALSE, boot.val=TRUE)
  summary(res2)
  plot(res2)
  
  ## example of PLS-PM in sensory analysis
  ## estimate a path model for the orange juice data
  data(orange)
  senso.inner <- matrix(c(0,0,0,1,0,0,1,1,0),3,3,byrow=TRUE)
  dimnames(senso.inner) <- list(c("PHYCHEM","SENSORY","HEDONIC"),
                              c("PHYCHEM","SENSORY","HEDONIC"))
  senso.outer <- list(1:9,10:16,17:112)
  senso.mod <- rep("A",3)
  res3 <- plspm.fit(orange, senso.inner, senso.outer, senso.mod, 
                scheme="centroid", scaled=TRUE)
  summary(res3)

  ## example of PLS-PM in multi-block data analysis
  ## estimate a path model for the wine data set 
  ## requires package FactoMineR 
  library(FactoMineR)
  data(wine)
  SMELL <- c(0,0,0,0)
  VIEW <- c(1,0,0,0)
  SHAKE <- c(1,1,0,0)
  TASTE <- c(1,1,1,0)
  wine.inner <- rbind(SMELL,VIEW,SHAKE,TASTE)
  wine.outer <- list(3:7,8:10,11:20,21:29)
  wine.mods <- rep("A",4)
  # using function plspm.fit (basic pls algorithm)
  res4 <- plspm.fit(wine, wine.inner, wine.outer, wine.mods, scheme="centroid")                
  plot(res4, what="all", arr.pos=.4, box.prop=.4, cex.txt=.8)
  }
}
