\name{summary.plspm}
\alias{summary.plspm}
\title{Summarizing partial least quares path modeling analysis}
\description{
  Summary method for objects of class \code{"plspm"}.
}
\usage{
  \method{summary}{plspm}(object, \dots)
}
\arguments{
  \item{object}{Object of class \code{"plspm"}.}
  \item{\dots}{Further arguments are ignored.}
}
\details{
  The function \code{summary.plspm} gives the detailed numerical results of the \code{\link{plspm}} function. 

  First, the model specification is shown followed by the definition of the blocks of variables.

  In second place the unidimensionality results of the blocks of manifest variables are displayed. 

  Then the outer (measurement) model and crossloadings are presented. 

  In fourth place, the inner (structural) model results as well as the correlations among latent variables are shown. 

  Next, a summary of the inner model and the goodness-of-fit indexes are presented.

  Finally, a table with the path effects of the structural relationships is given. 

  In addition, if bootstrap validation was performed, the bootstrapping results are listed.
}
