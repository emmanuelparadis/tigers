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
  \item{method}{an integer (1 or 2) specifying the triangulation
    method.}
}
\details{
  The following methods are available:

  \itemize{
    \item 1: the triangles are created in successive order from the
    first appropriate angle (i.e., an ear) encountered in the polygon.
    \item 2: the triangles are created to favour regularly-shaped
    triangles.
  }

  Method 1 needs the polygon to be closed, whereas method 2 needs it to
  be open and in counterclockwise order (method 1 is not sensitive to
  this order). These requirements are checked before performing the
  triangulation and the polygon is changed internally (without warning
  since the original polygon is not modified) if necessary.
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
