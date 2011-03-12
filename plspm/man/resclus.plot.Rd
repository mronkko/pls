\name{resclus.plot}
\alias{resclus.plot}
\title{
  Colored dendrogram for communality and structural residuals according to the REBUS algorithm
}
\description{
  Plot function for objects of class \code{"hclust"} returned by \code{res.clus}
}
\usage{
  resclus.plot(x, k=2, col.up="black", col.down=rainbow(k), lty.up=2, lty.down=1, lwd.up=1, lwd.down=1,
               type="rectangle", knot.pos="mean", show.labels=FALSE, only.tree=FALSE, members, ...)
}
\arguments{
  \item{x}{An object of class \code{"hclust"} obtained from \code{res.clus}. }
  \item{k}{The number of classes.}
  \item{col.up}{The color to be used for the lines above the cut of the dendrogram. Default "black".}
  \item{col.down}{The colors to be used for the classes below the cut of the dendrogram.}
  \item{lty.up}{The line type above the cut of the dendrogram. 
                Line types are specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash).}
  \item{lty.down}{The line type for the classes below the cut of the dendrogram.}
  \item{lwd.up}{The line width above the cut of the dendrogram. Line width must be a positive number, defaulting to 1. }
  \item{lwd.down}{The line width for the classes (below the cut of the dendrogram).}
  \item{type}{The type of dendrogram. Can be either "rectangle" or "triangle".}
  \item{knot.pos}{Possition of the knots. Can be "mean", "bary", "left", "right" or "random".}
  \item{show.labels}{Logical value indicating whether the labels of the objects should be printed. Default \code{FALSE}.}
  \item{only.tree}{Logical value indicating whether only the dendrogram should be printed. Default \code{FALSE}.}
  \item{members}{\code{NULL} or a vector with length size of a dissimilarity structure as produced by \code{\link{dist}}.}
  \item{\dots}{Further arguments are ignored.}
}
\details{
  The function \code{resclus.plot} displays a dendrogram that helps to visualize the partitions in a dendrogram from \code{res.clus} results.
}
\references{
  Modified version of the functions created by Romain Fracois. Source code available at:
  \url{http://addictedtor.free.fr/graphiques/RGraphGallery.php?graph=79}
}
\seealso{\code{\link{res.clus}}, \code{\link{it.reb}}}
\examples{
  \dontrun{
  ## example of rebus analysis with simulated data
  data(sim.data)
  ## First compute GLOBAL model
  sim.mat <- matrix(c(0,0,0,0,0,0,1,1,0),3,3,byrow=TRUE)
  dimnames(sim.mat) <- list(c("Price","Quality","Satisfaction"),
                            c("Price","Quality","Satisfaction"))
  sim.sets <- list(c(1,2,3,4,5),c(6,7,8,9,10),c(11,12,13)) 
  sim.mod <- c("A","A","A")  ## reflective indicators
  sim.global <- plspm(sim.data, sim.mat, sim.sets, sim.mod)
  sim.global
  ## Then compute cluster analysis on the residuals of global model
  sim.res.clus <- res.clus(sim.global)
  ## Plot dendrogram with 4 classes specifying colors
  dev.new()
  resclus.plot(sim.res.clus, k=4, col.up="gray", col.down=c("orange","royalblue","green3","red"))
  ## Plot triangular dendrogram with 4 classes, labels shown, and widther lines
  dev.new()
  resclus.plot(sim.res.clus, k=4, lwd.down=2, col.up="gray", show.labels=TRUE, type="triangle")
  }
}
