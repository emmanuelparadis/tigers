\name{is.insidePolygon}
\alias{is.insidePolygon}
\title{Test If a Point Is Inside a Polygon}
\description{
  This function tests if a point is inside a polygon.
}
\usage{
is.insidePolygon(XY, points)
}
\arguments{
  \item{XY}{A two-column matrix giving the coordinates of a polygon.}
  \item{points}{a vector with two values giving the coordinates of a
    point, or a matrix with two columns.}
}
\details{
  The algorithm is based on ``ray-tracing'': a segment is traced between
  \code{points} and an arbitrary point far from the polygon. If this
  segment intersects an odd number of edges of the polygon, then
  \code{points} is inside the polygon.

  The polygon must be open and can be in either clockwise or
  counterclockwise order. If the polygon is closed, it is modified
  internally without warning (the original polygon is not modified).
}
\value{a logical vector.}
\author{Emmanuel Paradis}
\seealso{\code{\link{is.open}}}
\examples{
XY <- rbind(c(0, 0), c(0, 1), c(1, 1), c(1, 0))
stopifnot(is.insidePolygon(XY, c(0.5, 0.5)))
stopifnot(!is.insidePolygon(XY, c(1.5, 1.5)))
}
\keyword{manip}