\name{fast2waytable}
\alias{fast2waytable}
\title{Fast Two-Way Contingency Tables}
\description{
  Contingency tables of two vectors in a much faster way than the
  default \code{\link[base]{table}}.
}
\usage{
fast2waytable(x, y, levels = NULL)
}
\arguments{
  \item{x, y}{two vectors of the same length with integers (or values
    that can be coerced to).}
  \item{levels}{the unique values found in both vectors.}
}
\details{
  If \code{levels} is known (and given), the running times are
  considerably shorter (up to several ten times; see examples).

  \code{NA}'s are not handled in this version.
}
\value{
  a square matrix with the counts.
}
\author{Emmanuel Paradis}
\seealso{\code{\link[base]{table}}, \code{\link[base]{tabulate}}}
\examples{
z <- 11:20 * 10L
n <- 13e6
x <- sample(z, n, TRUE)
y <- sample(z, n, TRUE)
system.time(res1 <- fast2waytable(x, y)) # ~ 0.4 sec
system.time(res2 <- fast2waytable(x, y, z)) # ~ 0.02 sec
system.time(res3 <- table(x, y)) # ~ 1.8 sec
all(res1 == res2)
all(res1 == res3)
}
\keyword{manip}
