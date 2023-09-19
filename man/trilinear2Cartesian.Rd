\name{trilinear2Cartesian}
\alias{trilinear2Cartesian}
\alias{Cartesian2trilinear}
\title{Trilinear Coordinates}
\description{
  \code{trilinear2Cartesian} calculates the coordinates of a point
  inside a triangle given three values interpreted as proportions.

  \code{Cartesian2trilinear} does the reverse operation.
}
\usage{
trilinear2Cartesian(p, X)
Cartesian2trilinear(xy, X)
}
\arguments{
  \item{p}{a vector with three numeric values (see details).}
  \item{X}{a numeric matrix with 3 rows and 2 columns giving the
    coordinates of the triangle.}
  \item{xy}{a vector with two numeric values (Cartesian coordinates).}
}
\details{
  The values in \code{p} do not need to sum to one since they are scaled
  internally.

  The triangle defined by \code{X} can be of any type. The coordinates
  returned by \code{trilinear2Cartesian} is always inside the triangle.

  \code{Cartesian2trilinear} does not check if \code{xy} is inside the
  triangle.
}
\value{
  \code{trilinear2Cartesian} returns a numeric matrix with a single row
  and two columns giving the coordinates of the point.

  \code{Cartesian2trilinear} returns a  a numeric matrix with a single row
  and three columns.
}
\author{Emmanuel Paradis}
\references{
  \url{https://en.wikipedia.org/wiki/Trilinear_coordinates}
}
\examples{
## rectangular triangle (counterclockwise):
X <- rbind(c(0, 0), c(0, 1), c(1, 0))
plot(X, , "n", asp = 1)
polygon(X)

h <- sqrt(2) # hypothenuse length

points(trilinear2Cartesian(c(1, 1, 1), X)) # incenter
points(trilinear2Cartesian(c(1, h, h), X), pch = 2) # centroid
points(trilinear2Cartesian(c(h, 1, 1), X), pch = 3) # symmedian point
## the 3 midpoints:
points(trilinear2Cartesian(c(0, h, h), X), pch = 7)
points(trilinear2Cartesian(c(1, 0, h), X), pch = 7)
points(trilinear2Cartesian(c(1, h, 0), X), pch = 7)

legend("topright", ,
       c("incenter", "centroid", "symmedian point", "midpoints"),
       pch = c(1:3, 7))

f <- c(0.1, 0.3, 0.6)
o <- trilinear2Cartesian(f, X)
p <- Cartesian2trilinear(o, X)
p - f # < 1e-15
stopifnot(all.equal(as.vector(p), f))
}
\keyword{manip}
