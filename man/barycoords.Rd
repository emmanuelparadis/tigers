\name{barycoords}
\alias{barycoords}
\title{Computes Barycentric Coordinates}
\description{
  The barycentric coordinates of a point inside a polygon are weighted
  coordinates of the vertices of this polygon. The algorithm implemented
  in this function works for any concave or convex polygon (Hormann and
  Floater, 2006).
}
\usage{
barycoords(XY, point)
}
\arguments{
  \item{XY}{A two-column matrix giving the coordinates of a polygon.}
  \item{point}{a vector with two values giving the coordinates of a point.}
}
\details{
  If the polygon is a triangle, the \code{\link{trilinear2Cartesian}}
  can be used instead.

  The polygon must be open (see \code{\link{is.open}}), and can be
  either in clockwise or in counterclockwise order (see
  \code{\link{is.clockwise}}).

  For the moment, the function is not vectorized with respect to
  \code{point}, so it must be called for each point separately (see
  examples). This is likely to change in the future.
}
\value{
  a numeric vector with length equal to the number of vertices in the
  polygon.
}
\references{
  Hormann, K. and Floater, M. S. (2006) Mean value coordinates for
  arbitrary planar polygons. \emph{ACM Transactions on Graphics}
  \bold{25}, 1424--1441. <doi:10.1145/1183287.1183295>
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{trilinear2Cartesian}}
}
\examples{
## a square:
xy <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))

## a small function to get the coordinates directly:
f <- function(Pxy) barycoords(xy, Pxy)
## the CMYK scale:
F <- col2rgb(c("cyan", "magenta", "yellow", "black"))

n <- 1e5L
## random points in the square
Pxys <- matrix(runif(2 * n), n, 2)
system.time(res <- t(apply(Pxys, 1, f))) # < 1 sec
colnames(res) <- as.character(1:4)

## all rows should (approximately) sum to one:
all.equal(rowSums(res), rep(1, n), tol = 1e-15)

## transform the barycentric coordinates into colours:
COLS <- t(F \%*\% t(res)) / 255
rgbCOLS <- apply(COLS, 1, function(x) do.call(rgb, as.list(x)))
## add transparency:
rgbCOLS <- paste0(rgbCOLS, "33")
## plot the results:
plot(0:1, 0:1, "n", asp = 1, ann = FALSE, axes = FALSE)
points(Pxys, pch = ".", col = rgbCOLS, cex = 20)
## the visual effect is nicer with n <- 1e6L above and cex = 7
## in the last command

## the example below follows the same logic than the previous one

## an 8-vertex polygon:
xy <- cbind(c(0, 0.5, 1, 3, 1, 0.5, 0, -2),
            c(0, -2, 0, 0.5, 1, 3, 1, 0.5))

## random points in the square and in the 4 triangles:
Pxys <- rbind(matrix(runif(2 * n), n, 2),
              rpit(n, xy[1:3, ]),
	      rpit(n, xy[3:5, ]),
	      rpit(n, xy[5:7, ]),
	      rpit(n, xy[c(7:8, 1), ]))

system.time(res <- t(apply(Pxys, 1, f))) # < 5 sec

colnames(res) <- as.character(1:8)
F <- col2rgb(c("black", "red", "orange", "green",
               "yellow", "blue", "purple", "white"))
## F <- col2rgb(rainbow(8)) # alternative
COLS <- t(F \%*\% t(res)) / 255.001
rgbCOLS <- apply(COLS, 1, function(x) do.call(rgb, as.list(x)))
rgbCOLS <- paste0(rgbCOLS, "33") # add transparency
plot(xy, , "n", asp = 1, ann = FALSE, axes = FALSE)
points(Pxys, pch = ".", col = rgbCOLS, cex = 5)
}
\keyword{manip}
