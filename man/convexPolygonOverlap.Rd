\name{convexPolygonOverlap}
\alias{convexPolygonOverlap}
\title{Overlap of Two Convex Polygons}
\description{
  Find the intersection of two convex polygons.
}
\usage{
convexPolygonOverlap(A, B)
}
\arguments{
  \item{A, B}{two two-column matrices giving the coordinates of two
    polygons.}
}
\details{
  The intersection of two overlapping convex polygons is a single convex
  polygon.

  The two input polygons must be in clockwise order.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{is.clockwise}}, \code{\link{polygonOverlap}}}
\examples{
X <- matrix(rnorm(3800), ncol = 2)
A <- X[chull(X), ]
Y <- matrix(rnorm(3800), ncol = 2)
B <- Y[chull(Y), ]

plot(rbind(A, B), type = "n", asp = 1)
polygon(A)
COLS <- c("blue", "red")
text(A, labels = 1:nrow(A), font = 2, cex = 1.5, col = COLS[1])
polygon(B)
text(B, labels = 1:nrow(B), font = 2, cex = 1.5, col = COLS[2])
legend("topleft", , c("A", "B"), text.font = 2, text.col = COLS)
O <- .Call("convexPolygonOverlap_Call", A, B)
polygon(O, border = NA, col = rgb(1, 1, 0, 0.5))
}
\keyword{manip}