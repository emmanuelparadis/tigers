\name{axisMap}
\alias{axisMap}
\title{Draw axes on Maps}
\description{
  Draw checkerboard-style axes on a map with appropriate scales and
  annotations.
}
\usage{
axisMap(latitude = FALSE, width = 0.05, len = 1,
        cols = c("black", "white"), \dots)
}
\arguments{
  \item{latitude}{by default axes are drawn for the longitudes (above
    and below the map). Use \code{latitude = TRUE} to draw axes for the
    latitudes (left and right sides of the map).}
  \item{width}{the width of the bands.}
  \item{len}{the length (increment) of the bands.}
  \item{cols}{the alternate colours of the bands.}
  \item{\dots}{further arguments passed to \code{axis} (see examples).}
}
\details{
  The axes are drawn as bands of width given by the argument of the same
  name and alternate colors according to \code{len}.
}
\author{Emmanuel Paradis}
\examples{
n <- 100
plot(runif(n, -180, 180), runif(n, -90, 90), pch = 3)
axisMap(len = 10)
## 'las' is passed with '...'
axisMap(TRUE, len = 10, las = 0)
}
\keyword{aplot}
