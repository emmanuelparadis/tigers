\name{area}
\alias{area}
\title{Area of Polygon}
\description{
  This function computes the area of a polygon with Euclidean coordinates
  (e.g., UTM).
}
\usage{
area(x, y = NULL)
}
\arguments{
  \item{x, y}{the coordinates of the points given in the usual way in
    \R.}
}
\details{
  The unit of the area are the squared unit of input coordinates by
  default.
}
\value{
  a single numeric value giving the area of the polygon.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{geod}}}
\examples{
XY <- rbind(c(0, 0),
            c(1, 0),
            c(.25, .25),
            c(.5, .5),
            c(1.2, .8),
            c(1, .78),
            c(0, 1))
area(XY)
}
\keyword{manip}