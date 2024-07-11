\name{buffer}
\alias{buffer}
\title{Buffer Around Polygons}
\description{
  Define a buffer zone around a polygon so that the outer points of the
  buffer are at a given distance to the nearest point of the input
  polygon.
}
\usage{
buffer(x, r, smoothing = 1, unit = "degree", quiet = FALSE)
}
\arguments{
  \item{x}{a two-column numeric matrix giving the coordinates of the
    vertices of the polygon.}
  \item{r}{a single numeric value giving the distance of the buffer
    zone.}
  \item{smoothing}{a numeric value coding for the smooting of the arcs
    of the buffer. Larger values result in less smoothing; set
    \code{smoothing = Inf} for no smoothing at all (see examples).}
  \item{unit}{a character string, one of \code{"degree"}, \code{"km"},
    or \code{"m"}. This argument defines a coefficient multiplied by the
    previous one to give better results.}
  \item{quiet}{a logical value specifying whether to print the progress
    of the computations.}
}
\details{
  This version assumes that the coordinates in \code{x} are Euclidean
  (see examples for conversion from longitude-latitude to UTM).

  The code may not work very well on polygons with complicated shapes
  and when the buffer size (argument \code{r}) is too large. This is
  still in progress.
}
\value{
  a two-column numeric matrix.
}
\author{Emmanuel Paradis}
\examples{
data(Anduki)
r0 <- 0.005
z <- buffer(Anduki, r0)
plot(z, , "n", asp = 1, xlab = "Longitude", ylab = "Latitude")
R <- seq(r0, r0/5, length.out = 5)
COLS <- R * 1000 + 1
for (i in seq_along(R))
    polygon(buffer(Anduki, R[i], quiet = TRUE), col = COLS[i])
polygon(Anduki, col = "white")
legend("topleft", , c(paste("r =", R), "Anduki Forest Reserve"),
       pt.bg = c(COLS, "white"), pt.cex = 2.5, pch = 22, bty = "n")
title("Buffers around the Anduki Forest Reserve")

## a more realistic application with the same data:
x <- lonlat2UTM(Anduki) # convert to UTM coordinates
z <- buffer(x, 100)
plot(z, , "n", asp = 1, xlab = "Easting", ylab = "Northing")
title("Buffer zone of 100 m (UTM coordinates)")
polygon(z, col = "lightgrey", lwd = 1/3)
polygon(x, col = "white", lwd = 1/3)


## a concave hexagon
x <- rbind(c(0, 0), c(0, 1), c(0.5, 0.75),
           c(1, 1), c(1, 0), c(.5, .5))
r0 <- 0.1
z <- buffer(x, r0)
plot(z, , "n", asp = 1)
R <- seq(r0, r0/10, length.out = 10)
COLS <- R * 100 + 1
for (i in seq_along(R))
    polygon(buffer(x, R[i], quiet = TRUE), col = COLS[i])
polygon(x, col = "white")

\dontrun{
## more fancy:
plot(z, , "n", asp = 1)
R <- seq(r0, r0/10, length.out = 1000)
COLS <- rainbow(1000)
for (i in seq_along(R))
    polygon(buffer(x, R[i], quiet = TRUE), col = COLS[i], border = NA)
polygon(x, col = "white", border = NA)
}

## Effect of smoothing
layout(matrix(1:4, 2, 2, TRUE))
for (sm in c(1, 30, 50, Inf)) {
    z <- buffer(x, 0.5, sm, quiet = TRUE)
    plot(z, , "n", asp = 1, ann = FALSE, axes = FALSE, )
    polygon(z, col = "lightblue")
    title(paste("smoothing =", sm))
    polygon(x, col = "white")
}
layout(1)
}
\keyword{manip}
