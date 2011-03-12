\name{plot.plspm}
\alias{plot.plspm}
\title{Plot plspm path diagrams}
\description{
  Plot method for objects of class \code{"plspm"}
}
\usage{
  \method{plot}{plspm}(x, what = "inner", how = "joint", arr.pos = 0.5, 
     box.prop = 0.5, box.cex = 1, cex.txt = 1, \dots)
}
\arguments{
  \item{x}{An object of class \code{"plspm"} returned by \code{\link{plspm}}.}
  \item{what}{A character string indicating what models should be ploted. Options: \code{"inner"} (defualt value), \code{"weights"}, \code{"loadings"}, and \code{"all"}.}
  \item{how}{A character string indicating how the diagrams should be ploted. Options: \code{"joint"} (defualt value), and \code{"split"}.}
  \item{arr.pos}{Relative position of arrowheads on arrows.}
  \item{box.prop}{Length/width ratio of boxes.}
  \item{box.cex}{Relative size of text in boxes.}
  \item{cex.txt}{Relative size of text on arrows.}
  \item{\dots}{Further arguments are ignored.}
}
\details{
  The function \code{plot.plspm} allows to display path diagrams of a plspm model. \cr
  If \code{what="inner"} only a diagram of the inner model is displayed. \cr
  If \code{what="weights"} diagrams of the outer weights for each block of variables are displayed. \cr
  If \code{what="loadings"} diagrams of the loadings for each block of variables are displayed. \cr
  If \code{what="all"} diagrams of the inner model, outer weights, and loadings are displayed. \cr
  When \code{how="split"} all diagrams are displayed separately by each block of variables: \cr
}
\note{
  Function \code{plot.plspm} is based on the function \code{\link{plotmat}} of package \code{diagram}. \cr
  \url{http://cran.r-project.org/web/packages/diagram/vignettes/diagram.pdf}  
}
\seealso{\code{\link{plspm}}}
\examples{
  \dontrun{
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
  res2 <- plspm(satisfaction, sat.inner, sat.outer, sat.mod, scheme="centroid", 
                scaled=FALSE)
  ## plot diagram of the inner model
  plot(res2)
  ## plot diagrams of both the inner model and outer model (loadings and weights)
  plot(res2, what="all")
  }
}
