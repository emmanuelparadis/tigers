\name{random_point_in_triangle}
\alias{random_point_in_triangle}
\alias{rpit}
\title{Random Points in Triangle}
\description{
  Generates random points inside a triangle using Osada et al.'s (2002,
  Sect. 4.2) method.
}
\usage{
  random_point_in_triangle(n, X, rfun1 = runif, rfun2 = runif)
  rpit(n, X, rfun1 = runif, rfun2 = runif)
}
\arguments{
  \item{n}{an integer giving the number of points to generate.}
  \item{X}{a numeric matrix with 3 rows and 2 columns giving the
    coordinates of the triangle.}
  \item{rfun1}{a function generating random values in [0,1]. By default,
    the values are generated under a uniform distribution.}
  \item{rfun2}{same as the previous argument (see details).}
}
\details{
  By default, the points are uniformly distributed in the triangle.
  The \code{\link[stats]{Beta}} function offers an interesting
  alternative to  generate points concentrated in a specific part of the
  triangle (see examples).
}
\value{
  A numeric matrix with \code{n} rows and two columns giving the
  coordinates of the points.
}
\references{
  Osada, R., Funkhouser, T., Chazelle, B., and Dobkin, D. (2002) Shape
  distributions. \emph{ACM Transactions on Graphics}, \bold{21},
  807--832. <doi:10.1145/571647.571648>
}
\author{Emmanuel Paradis}
\examples{
## a random triangle in [0,1]^2:
P <- matrix(runif(6), 3, 2)

## n points uniformly distributed in the triangle P:
n <- 10000
x <- rpit(n, P)

layout(matrix(1:2, 1))

plot(P, type = "n", asp = 1)
polygon(P, col = "yellow", border = NA)
points(x, pch = ".", col = "blue")

## using Beta distributions:
foo <- function(n) rbeta(n, 1, 10)
bar <- function(n) rbeta(n, 1, 1)
y <- rpit(n, P, foo, bar)

plot(P, type = "n", asp = 1)
polygon(P, col = "yellow", border = NA)
points(y, pch = ".", col = "blue")

layout(1)
}
\keyword{datagen}
