\name{RMA}
\alias{RMA}
\title{Reduced Major Axis}
\description{
  Computes the coefficients of the reduced major axis (RMA) of a set of
  points.
}
\usage{
  RMA(x, y = NULL)
}
\arguments{
  \item{x, y}{the coordinates of the points given in the usual way in \R.}
}
\details{
  The RMA is found by solving a polynomial equation of degree two, so
  there are actually two solutions which are both returned. It is
  usually straightforward to find the appropriate solution.
}
\references{
  \url{https://mathworld.wolfram.com/LeastSquaresFittingPerpendicularOffsets.html}
}
\value{
  a matrix with two rows and two columns named alpha and beta for the
  intercepts and slopes, respectively.
}
\author{Emmanuel Paradis}
\examples{
x <- 1:1000
y <- x + rnorm(1000, 5)
RMA(x, y) # same than RMA(cbind(x, y))
}
\keyword{robust}