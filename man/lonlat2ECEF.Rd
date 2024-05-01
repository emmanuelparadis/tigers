\name{lonlat2ECEF}
\alias{lonlat2ECEF}
\alias{ECEF2lonlat}
\title{Conversions of Coordinates}
\description{
  Conversion between lon-lat and ECEF (Earth-centered, Earth-fixed)
  coordinates.
}
\usage{
lonlat2ECEF(lon, lat = NULL, alt = 0, as.matrix = TRUE)
ECEF2lonlat(x, y = NULL, z = NULL)
}
\arguments{
  \item{lon, lat}{coordinates given as two vectors or as a matrix (or
    a data frame).}
  \item{alt}{a vector of altitudes.}
  \item{as.matrix}{a logical value specifying whether to return the
    results in a matrix. The default is to return them in a list.}
  \item{x, y, z}{three numeric vectors of the same length, or only
    \code{x} as a three-column matrix.}
}
\details{
  This function uses an ellopsoid model of Earth shape with the
  following parameters: equatorial radius = 6,378,137 m, polar radius
  = 6,356,752 m.
}
\value{
  \code{lonlat2ECEF} returns a matrix with three columns or a list with
  three vectors; \code{ECEF2lonlat} returns a matrix with three columns.
}
\author{Emmanuel Paradis}
\references{
  Blewitt, G. (2024) An improved equation of latitude and a global
  system of graticule distance coordinates. \emph{Journal of Geodesy},
  \bold{98}, 6.
}
\examples{
## From Blewitt's Table 1:
DF <- expand.grid(alt = c(1, 4, 10, 40) *1e3,
                  lat = seq(15, 75, 15),
                  lon = 0)
DF <- as.matrix(DF[3:1])

x <- lonlat2ECEF(DF[, -3], alt = DF[, 3], as.matrix = TRUE)

DFbis <- ECEF2lonlat(x)
xbis <- lonlat2ECEF(DFbis[, -3], alt = DFbis[, 3], as.matrix = TRUE)

summary(x - xbis)
summary(DF - DFbis)
}
\keyword{manip}
