\name{great_circle_line}
\alias{great_circle_line}
\alias{gcl}
\title{Great Circle Line}
\description{
  This function calculates the coordinates of the line on the surface of
  a sphere between two points. All coordinates are in decimal degrees.

  \code{gcl} is simply an alias.
}
\usage{
great_circle_line(x0, y0, x1, y1, linear = FALSE, npoints = 100)
gcl(x0, y0, x1, y1, linear = FALSE, npoints = 100)
}
\arguments{
  \item{x0, y0, x1, y1}{the coordinates of the two points similar to
    \code{\link[graphics]{segments}}.}
  \item{linear}{a logical value.}
  \item{npoints}{an integer giving the number of points where the
    coordinates are calculated (should be at least two).}
}
\details{
  The interval between \code{x0} and \code{x1} is split into regular
  segments, then the latitudes are computed, by default, using a
  great circle formula (Chamberlain and Duquette, 2007).

  If \code{linear = TRUE}, the coordinates are treated as linear (i.e.,
  Euclidean).
}
\value{
  a numeric matrix with two columns and colnames 'x' and 'y'.
}
\author{Emmanuel Paradis}
\references{
  Chamberlain, R. G. and Duquette, W. H. (2007) Some algorithms for
  polygons on a sphere. JPL Open Repository. <doi:2014/41271>
}
\seealso{\code{\link{geod}}}
\examples{
X1 <- 3; Y1 <- 49 # Paris
X2 <- 101; Y2 <- 13 # Bangkok
## if (require(maps))  map() else
plot(c(-180, 180), c(-90, 90), "n")
text(X1, Y1, "Paris")
text(X2, Y2, "Bangkok")
lines(gcl(X1, Y1, X2, Y2), col = "blue", lwd = 2)
lines(gcl(X1, Y1, X2, Y2, linear = TRUE), col = "red", lwd = 2)

## assess the error implied by using linear interpolation for the
## diagonal of a 1 degree by 1 degree square near the equator:
xya <- gcl(0, 0, 1, 1)
xyb <- gcl(0, 0, 1, 1, TRUE)
## the error in degrees:
error <- xya[, "y"] - xyb[, "y"]
plot(xya[, "x"], error * 3600, "o",
     xlab = "Longitude (degrees)", ylab = "Error (arc-seconds)")

## max (vertical) distance between these 2 curves:
geod(c(0.5, 0.5), c(0.5, 0.5 + max(error))) # ~6.5 m
## NOTE: the actual shortest (orthogonal) distance
## between these two curves is ~4.6 m
## (assuming the vertical distance helps to define a rectangular
## triangle, we have: 0.5 * sqrt(6.5^2 * 2)) ~ 4.6)

## NOTE2: dividing the coordinates by 10 results in dividing
## these deviations by 1000
}
\keyword{manip}
