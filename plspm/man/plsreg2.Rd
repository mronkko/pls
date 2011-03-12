\name{plsreg2}
\alias{plsreg2}
\title{PLS-R2: Partial Least Squares Regression 2}
\description{
  Calculates partial least squares regression for the multivariate case (i.e. more than one response variable)
}
\usage{
  plsreg2(X, Y, nc = 2)
}
\arguments{
  \item{X}{A numeric matrix or data frame containing the predictor variables.}
  \item{Y}{A numeric matrix or data frame containing the predictand variables.}
  \item{nc}{The number of extracted PLS components (2 by default)}
}
\details{
  The minimum number of PLS components \code{nc} is 2.

  The data is scaled to standardized values (mean=0, variance=1). 

  No missing data are allowed.

  Argument \code{Y} must contain more than one variable. If \code{Y} is a vector, you may use the function \code{\link{plsreg1}}.
}
\value{
  An object of class \code{"plsreg2"}, basically a list with the following elements:
  \item{x.scores}{components of the predictor variables.}
  \item{x.loads}{loadings of the predictor variables.}
  \item{y.scores}{components of the predictand variables.}
  \item{y.loads}{loadings of the predictand variables.}
  \item{raw.wgs}{weights to calculate the PLS scores with the deflated matrices of predictor variables.}
  \item{mod.wgs}{modified weights to calculate the PLS scores with the matrix of predictor variables.}
  \item{cor.tx}{modified weights to calculate the PLS scores with the matrix of predictor variables.}
  \item{cor.ty}{modified weights to calculate the PLS scores with the matrix of predictor variables.}
  \item{std.coef}{Vector of standardized regression coefficients.}
  \item{coeffs}{Vector of regression coefficients (used with the original data scale).}
  \item{y.pred}{Vector of predicted values.}
  \item{resid}{Vector of residuals.}
  \item{expvar}{table with R-squared coefficients.}
  \item{Q2}{table of Q2 indexes (i.e. leave-one-out cross validation).}
  \item{Q2cum}{table of cummulated Q2 indexes.}
  \item{VIP}{Variable Importance for Projection.}
}
\references{ 
  Geladi, P., and Kowlaski, B. (1986) Partial Least Squares Regression: A Tutorial. \emph{Analytica Chimica Acta}, \bold{185}, pp. 1-17.

  Hoskuldsson, A. (1988) PLS Regression Methods. \emph{Journal of Chemometrics}, \bold{2}, pp. 211-228.

  Tenenhaus, M. (1998) \emph{La Regression PLS. Theorie et Pratique}. Editions TECHNIP, Paris.
  
  Valencia, J.L., Diaz-Llanos, F.J. (2004) \emph{Metodos de Prediccion en Situaciones Limite}. Editorial La Muralla, S.A. Madrid.
}
\author{Gaston Sanchez}
\seealso{\code{\link{print.plsreg2}},  \code{\link{plot.plsreg2}},  \code{\link{plsreg1}}.}
\examples{
  \dontrun{
  ## example of PLSR2 with the vehicles dataset
  data(vehicles)
  pls2 <- plsreg2(vehicles[,1:12], vehicles[,13:16])
  pls2
  plot(pls2)
  }
}
