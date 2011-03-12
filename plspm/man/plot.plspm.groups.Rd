\name{plot.plspm.groups}
\alias{plot.plspm.groups}
\title{Plot results from a group comparison test in PLS-PM}
\description{
  Plot method for objects of class \code{"plspm.groups"}
}
\usage{
  \method{plot}{plspm.groups}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{"plspm.groups"} returned by \code{\link{plspm.groups}}.}
  \item{\dots}{Further arguments are ignored.}
}
\details{
  The function \code{plot.plspm.groups} displays a barplot with the path coefficients of the compared groups.
}
\seealso{\code{\link{plspm.groups}},  \code{\link{plspm}}}
\examples{
  \dontrun{
  ## example with customer satisfaction analysis
  ## group comparison based on the segmentation variable "gender"
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
  pls <- plspm(satisfaction, sat.inner, sat.outer, sat.mod, scheme="factor", 
                scaled=FALSE)
  ## permutation test with 100 permutations
  res.group <- plspm.groups(pls, satisfaction$gender, 
               method="permutation", reps=100)
  res.group
  plot(res.group)
  }
}
