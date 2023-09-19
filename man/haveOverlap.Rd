\name{haveOverlap}
\alias{haveOverlap}
\alias{samePolygons}
\title{Compare Two Polygons}
\description{
  These functions compare two polygons.
}
\usage{
haveOverlap(A, B)
samePolygons(A, B, digits = 10)
}
\arguments{
  \item{A, B}{Two two-column matrices giving the coordinates of two
    polygons.}
  \item{digits}{the number of digits considered when comparing the
    coordinates.}
}
\value{a single logical value}
\author{Emmanuel Paradis}
\seealso{\code{\link{redundantVertices}}}
\keyword{manip}
