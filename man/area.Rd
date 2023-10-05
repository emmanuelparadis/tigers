\name{area}
\alias{area}
\title{Area of Polygon}
\description{
  This function computes the area of a polygon where the coordinates are
  either Euclidean (e.g., UTM) or angular (e.g., longitude and latitude).
}
\usage{
area(x, y = NULL, lonlat = FALSE, spheroid = TRUE)
}
\arguments{
  \item{x, y}{the coordinates of the points given in the usual way in
    \R.}
  \item{lonlat}{if \code{TRUE}, the coordinates are treated as
    longitudes and latitudes (in decimal degrees). By default, Euclidean
    coordinates are assumed.}
  \item{spheroid}{by default, the algorithm assuming a spheroid-Earth
    model is used. If \code{FALSE}, a spherical-Earth model is assumed
    (see details; ignored if \code{lonlat = FALSE}).}
}
\details{
  If \code{lonlat = TRUE}, the coordinates must be in decimal degrees
  with longitude as the first column, and latitude as the second one.

  In the case of angular coordinates (if \code{lonlat = TRUE}), the
  default method assumes the Earth is a spheroid with parameters taken
  from IERS Conventions (2010): the equatorial radius is
  6,378,136.6 m, and the inverse of the Earth's flattening is
  298.25642 (Paradis 2023).

  The alternative algorithm (if \code{spheroid = FALSE}) is from
  Chamberlain and Duquette (2007). These authors reported that this
  formula does not work with very large regions of the Earth,
  particularly including the pole(s). However, it seems that if the
  result is multiplied by two, a correct answer is found (see examples
  below), although it is not clear whether this is a robust result. The
  radius of the Earth in this formula is fixed internally to 6,371,008 m
  which is the mean radius (i.e., (2*a + b)/3, with a and b are the
  equatorial and polar radius, respectively). An R version of this
  algorithm is given below.

  Some tests (not reported here) showed that the results from
  Chamberlain and Duquette's formula agree with the formula given in
  Paradis (2020) (which works only with spherical trapezes) to 12 digits
  after the decimal place.

  The unit of the area are the squared unit of input coordinates by
  default, or squared metre if angular coordinates are used.
}
\references{
  Chamberlain, R. G. and Duquette, W. H. (2007) Some algorithms for
  polygons on a sphere. JPL Open Repository. <doi:2014/41271>

  IERS Conventions (2010) IERS Technical Note 36. Verlag des Bundesamts
  f\"ur 232 Kartographie und Geod\"asie, Frankfurt am Main.
  \url{https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html}.

  Paradis, E. (2020) Modelling transition in land cover highlights
  forest losses and gains in Southeast Asia. \emph{Biodiversity and
    Conservation}, \bold{29}, 2539--2551.

  Paradis, E. (2023) Fast and efficient computation of geodesic polygon
  areas. Submitted to \emph{Journal of Geodesy}.
}
\value{
  a single numeric value.
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
area(XY, lonlat = TRUE)

## the R version of Chamberlain and Duquette's algorithm:
foo <- function(lon, lat, R = 6378.1) {
    x <- lon * pi/180 # degrees -> radians
    y <- lat * pi/180 #
    n <- length(x) # number of vertices in the polygon
    prv <- c(n, 1:(n - 1))
    nxt <- c(2:n, 1)
    abs(0.5 * R^2 * sum((x[prv] - x[nxt]) * sin(y)))
}

## a region 1 degree by 1 degree near the equator:
foo(c(0, 1, 1, 0), c(0, 0, 1, 1))

## same than: area(c(0, 1, 1, 0), c(0, 0, 1, 1), lonlat = TRUE)

## what is the approximation if we consider this region a square
## of side 111.19 km (see ?geod):
111.19^2

## this "triangle" is a quarter of the northern hemisphere:
## (see above Details for why the result is multiplied by 2)
2 * foo(c(0, 90, 0), c(0, 0, 90))

## same than: area(c(0, 90, 0), c(0, 0, 90), lonlat = TRUE)

## the area of a sphere of radius r is: 4 pi r^2, thus a quarter
## of an hemisphere is:
(4 * pi * 6378.1^2) / 8
}
\keyword{manip}