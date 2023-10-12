\name{chullPolygon}
\alias{chullPolygon}
\title{Convex Hull of Polygon}
\description{
  Finds the convex hull of a polygon.

  Note that the function \code{chull} (see link below) finds the convex
  hull of a set of points and is about twice slower than the present one
  when applied to a polygon.
}
\usage{
chullPolygon(x, y = NULL)
}
\arguments{
  \item{x, y}{the coordinates of the points given in the usual way in \R.}
}
\details{
  This internal implementation requires the polygon to be open and in
  clockwise order (a crash will happen otherwise). Clockwise order is
  checked and possibly handled before calling the C code.
}
\references{
  Graham, R. L. and Yao, F. F. (1983) Finding the convex hull of a
  simple polygon. \emph{Journal of Algorithms}, \bold{4}, 324--331.
  <doi:10.1016/0196-6774(83)90013-5>
}
\value{
  a vector of integers which give the indices of the vertices of the
  input polygon defining the convex hull.
}
\seealso{
  \code{\link[grDevices]{chull}}
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
(i <- chullPolygon(XY))
plot(XY, type = "n", asp = 1)
polygon(XY, lwd = 5, border = "lightgrey")
text(XY, labels = 1:nrow(XY), cex = 2/1.5)
polygon(XY[i, ], border = "blue", lty = 2, lwd = 3)
}
\keyword{manip}