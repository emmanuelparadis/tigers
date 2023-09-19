\name{distance_to_line}
\alias{distance_to_line}
\alias{dtl}
\alias{distance_to_arc}
\alias{dta}
\title{Distance to Line}
\description{
  These functions calculate the shortest distances from a set of points
  to a line (in Euclidean coordinates) or an arc (in angular coordinates).

  \code{dtl} and \code{dtl} are aliases to \code{distance_to_line} and
  \code{distance_to_arc}, respectively.
}
\usage{
distance_to_line(x, y = NULL, x0, y0, x1, y1,
                 alpha = NULL, beta = NULL)
dtl(x, y = NULL, x0, y0, x1, y1, alpha = NULL, beta = NULL)
distance_to_arc(x, y = NULL, x0, y0, x1, y1, prec = 0.001)
dta(x, y = NULL, x0, y0, x1, y1, prec = 0.001)
}
\arguments{
  \item{x, y}{the coordinates of the points given in the usual way in \R.}
  \item{x0, y0, x1, y1}{the coordinates of two points defining the line
    similar to \code{\link[graphics]{segments}}. These are ignored if
    \code{alpha} and \code{beta} are given.}
  \item{alpha, beta}{alternatively to the previous arguments, the
    parameters of the line (\code{beta} is the slope).}
  \item{prec}{the precision of the estimated distances (see details).}
}
\details{
  \code{distance_to_line} uses Euclidean geometry (see references). The
  coordinates can be in any units.

  \code{distance_to_arc} uses distances along arcs on the (Earth)
  sphere. The coordinates must be in decimal degrees. The calculations
  are done by iterations using intervals of decreasing lengths along the
  arc. The iterations are stopped when the required precision is reached
  (see argument \code{prec}).
}
\value{
  a numeric vector giving the distances; \code{distance_to_line} returns
  them in the same unit than the input data; \code{distance_to_arc}
  returns them in kilometres (km).
}
\author{Emmanuel Paradis}
\references{
  \url{https://en.wikipedia.org/wiki/Distance\_from\_a\_point\_to\_a\_line}
}
\seealso{
  \code{\link{great_circle_line}}, \code{\link{geoTrans}},
  \code{\link{geod}}
}
\examples{
## distance from the topleft corner of the unity square to the diagonal:
(d <- dtl(matrix(c(1, 0), , 2), NULL, 0, 0, 1, 1))
all.equal(d, sqrt(2)/2)

## see also ?great_circle_line
x <- y <- 0:10/10
dta(x, y, 0, 0, 1, 1)
}
\keyword{manip}
