\name{plsreg1}
\alias{plsreg1}
\title{PLS-R1: Partial Least Squares Regression 1}
\description{
  Calculates partial least squares regression for the univariate case (i.e. one response variable)
}
\usage{
  plsreg1(x, y, nc = 2, cv = FALSE)
}
\arguments{
  \item{x}{A numeric matrix or data frame with the predictor variables (which may contain missing data).}
  \item{y}{A numeric vector for the reponse or predictand variable.}
  \item{nc}{The number of extracted PLS components (2 by default).}
  \item{cv}{A logical value indicating whether cross-validation should be performed (\code{FALSE} by default).}
}
\details{
  The minimum number of PLS components \code{nc} to be extracted is 2.

  The argument \code{x} may contain missing data. Conversely, the argument \code{y} must not contain missing values.

  The data is scaled to standardized values (mean=0, variance=1). 

  The argument \code{cv} gives the option to perform leave-one-out cross validation to choose the best number of PLS components.
}
\value{
  An object of class \code{"plsreg1"}, basically a list with the following elements:
  \item{scores}{PLS components.}
  \item{x.loads}{loadings of the predictor variables.}
  \item{y.loads}{loadings of the predictand variable.}
  \item{u.scores}{u scores of the predictand variable.}
  \item{raw.wgs}{weights to calculate the PLS scores with the deflated matrices of predictor variables.}
  \item{mod.wgs}{modified weights to calculate the PLS scores with the matrix of predictor variables.}
  \item{std.coef}{Vector of standardized regression coefficients.}
  \item{coeffs}{Vector of regression coefficients (used with the original data scale).}
  \item{R2}{Vector of PLS R-squared.}
  \item{y.pred}{Vector of predicted values.}
  \item{resid}{Vector of residuals.}
  \item{cor.sco}{Correlations between the variables and the PLS components.}
  \item{T2}{Table of Hotelling T2 values (used to detect atypical observations).}
  \item{Q2}{Table with the cross validation results. Includes: PRESS, RSS, Q2, and cummulated Q2. Only available when \code{cv=TRUE}}
}
\references{
  Geladi, P., and Kowalski, B. (1986) Partial Least Squares Regression: A Tutorial. \emph{Analytica Chimica Acta}, \bold{185}, pp. 1-17.

  Tenenhaus, M. (1998) \emph{La Regression PLS. Theorie et Pratique}. Editions TECHNIP, Paris.
  
  Tenenhaus, M., Gauchi, J.-P., and Menardo, C. (1995) Regression PLS et applications. \emph{Revue de statistique appliquee}, \bold{43}, pp. 7-63. 

  Valencia, J.L., Diaz-Llanos, F.J. (2004) \emph{Metodos de Prediccion en Situaciones Limite}. Editorial La Muralla, S.A. Madrid.
}
\author{Gaston Sanchez}
\seealso{\code{\link{print.plsreg1}},  \code{\link{plot.plsreg1}},  \code{\link{plsreg2}}.}
\examples{
  \dontrun{
  ## example of PLSR1 with the vehicles dataset
  ## predictand variable: price of vehicles
  data(vehicles)
  pls1 <- plsreg1(vehicles[,1:12], vehicles[,13], cv=TRUE)
  pls1
  plot(pls1)
  }
}
