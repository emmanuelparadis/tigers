\name{polygonOverlap}
\alias{polygonOverlap}
\alias{decomposePolygon}
\title{Decomposition and Overlap of Polygons}
\description{
  \code{decomposePolygon} decomposes a polygon into convex subpolygons.

  \code{polygonOverlap} finds the intersection of two polygons.
}
\usage{
decomposePolygon(x, y = NULL, method = 1, quiet = FALSE)
polygonOverlap(A, B)
}
\arguments{
  \item{x, y}{the coordinates of the points given in the usual way in
    \R.}
  \item{method}{the method used for triangulation (see
    \code{\link{triangulate}}).}
  \item{quiet}{if the polygon is convex, a warning message is issued
    unless this option is switched to \code{TRUE}.}
  \item{A, B}{two two-column matrices giving the coordinates of two
    polygons.}
}
\details{
  Both functions require the polygons to be in counterclockwise order
  (which is checked and arranged internally if needed).

  The method in \code{decomposePolygon} is from Hertel and Mehlhorn
  (1983).

  The method in \code{polygonOverlap} is based on first decomposing the
  two polygons into convex subpolygons, then computing their
  intersections with \code{\link{convexPolygonOverlap}}. The results is
  a list of polygons. A different algorithm is sketched in Chamberlain
  and Duquette (2007).
}
\note{
  These two functions are still in development.
}
\value{
  \code{decomposePolygon} returns a two-column matrix with integers
  where each row gives the indices of two vertices of the input polygon
  defining a diagonal; the set of these diagonals define convex
  subpolygons.

  \code{polygonOverlap} returns a list of polygons each defined by a
  two-column numeric matrix giving the coordinates of the vertices.
}
\author{Emmanuel Paradis}
\references{
  Chamberlain, R. G. and Duquette, W. H. (2007) Some algorithms for
  polygons on a sphere. JPL Open Repository. <doi:2014/41271>

  Hertel, S. and Mehlhorn, K. (1983) Fast triangulation of simple
  polygons. In: \emph{Foundations of Computation Theory.} Ed. Karpinski,
  M. Springer, Berlin, pp. 207--218. <doi:10.1007/3-540-12689-9_105>
}
\seealso{\code{\link{convexPolygonOverlap}}, \code{\link{is.clockwise}}}
\examples{
## same polygon than in ?triangulate
XY <- rbind(c(0, 0), c(1, 0), c(.25, .25), c(.5, .5),
            c(1.2, .8), c(1, .78), c(0, 1))
decomposePolygon(XY) # similar to the output of triangulate()
## "lift up" one vertex:
XYb <- XY
XYb[6, 2] <- 1.2
decomposePolygon(XYb) # one diagonal less

## A is concave, B is convex:
A <- rbind(c(0, 1.5), c(2, 1), c(0.5, 1.5), c(2, 2))
B <- rbind(c(1, 0), c(3, 0), c(3, 3), c(1, 3))
AB <- polygonOverlap(A, B)
plot(rbind(A, B), , "n", asp = 1)
polygon(A)
polygon(B)
lapply(AB, polygon, col = "gold")
}
\keyword{manip}