\name{redundantVertices}
\alias{redundantVertices}
\alias{is.clockwise}
\alias{is.convex}
\alias{is.open}
\alias{revPolygon}
\title{Redundant Vertices in a Polygon}
\description{
  Tests and optionally correct for redundant vertices in a polygon.

  The other functions test some features of a polygon.

  \code{revPolygon} reverses the order of the vertices (i.e., swiching
  between clockwise and counterclockwise orders).
}
\usage{
redundantVertices(x, tol = 1e-8, check.only = FALSE)
is.clockwise(x)
is.convex(x)
is.open(x)
revPolygon(x, copy = TRUE)
}
\arguments{
  \item{x}{a two-column matrix.}
  \item{tol}{the tolerance to consider two vertices identical.}
  \item{check.only}{a logical value.}
  \item{copy}{by default, a new polygon is created; if \code{FALSE}, the
    vertex order is reversed within the same polygon.}
}
\details{
  If \code{check.only} is \code{TRUE}, the function prints the
  diagnostics and nothing is returned. Otherwise, the possibly corrected
  matrix is returned.
}
\references{
  The method for \code{is.clockwise} is from:

  \url{https://en.wikipedia.org/wiki/Curve_orientation}
}
\author{Emmanuel Paradis}
\seealso{\code{\link{haveOverlap}}}
\keyword{manip}
