\name{triangulate}
\alias{triangulate}
\title{Triangulate a Polygon}
\description{
  Performs the decomposition of a polygon into triangles.
}
\usage{
triangulate(x, y = NULL, method = 1)
}
\arguments{
  \item{x, y}{the coordinates of the points given in the usual way in
    \R.}
  \item{method}{an integer between 1 and 4 specifying the triangulation
    method.}
}
\details{
  The following methods are available:

  \itemize{
    \item 1: the triangles are created in successive order from the
    first appropriate angle (i.e., an ear) encountered in the polygon.
    \item 2: the triangles are created to favour thin triangles.
    \item 3: the triangles are created to favour fat triangles.
    \item 4: the triangles are created to favour regular-looking
      triangles based on their determinant.
  }

  These methods have different requirements: method 1 needs the polygon
  to be closed, whereas the other methods need it to be open; method 2
  needs the polygon to be in counterclockwise order, and method 3 needs
  it to be in clockwise order (the other methods are not sensitive to
  this order). These requirements are checked before performing the
  triangulation and the polygon is changed internally (without warning
  since the original polygon is not modified) if necessary.
}
\note{
  The internal codes need to be checked and tested again.
}
\references{
  Toussaint, G. (1991) Efficient triangulation of simple
  polygons. \emph{Visual Computer}, \bold{7}, 280--295.
  <doi:10.1007/BF01905693>
}
\value{
  a three-column matrix giving the indices of the vertices in each
  triangle (i.e., each row a is a triangle).
}
\author{Emmanuel Paradis}
\examples{
XY <- rbind(c(0, 0),
            c(1, 0),
            c(.25, .25),
            c(.5, .5),
            c(1.2, .8),
            c(1, .78),
            c(0, 1))
(tri <- triangulate(XY))
plot(XY, type = "n", asp = 1)
for (i in 1:nrow(tri))
    polygon(XY[tri[i, ], ], border = "white", col = "green", lwd = 2)
polygon(XY, lwd = 4, border = "lightgrey")
text(XY, labels = 1:nrow(XY), cex = 1.2)
}
\keyword{manip}