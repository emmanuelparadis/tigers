\name{HausdorffDistance}
\alias{HausdorffDistance}
\title{Hausdorff Distance}
\description{
  Computes the Hausdorff distance between two polygons. The distances
  can be directed (i.e., asymmetric) or not.
}
\usage{
  HausdorffDistance(A, B, directed = FALSE)
}
\arguments{
  \item{A, B}{two two-column matrices giving the coordinates of two
    polygons.}
  \item{directed}{a logical value. By default, the undirected distance
    is returned.}
}
\details{
  If \code{directed = TRUE}, the order of the two polygons is important.
}
\value{a single numeric value.}
\author{Emmanuel Paradis}
\examples{

A <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
B <- A
B[, 1] <- B[, 1] + 2
B[c(1, 4), 1] <- 1.15

plot(rbind(A, B), type = "n", asp = 1)
COLS <- c("blue", "red")
polygon(A, border = COLS[1], lwd = 3)
polygon(B, border = COLS[2], lwd = 3)
text(mean(A[, 1]), mean(A[, 2]), "A", font = 2, col = COLS[1])
text(mean(B[, 1]), mean(B[, 2]), "B", font = 2, col = COLS[2])

(H <- HausdorffDistance(A, B))
(HAB <- HausdorffDistance(A, B, TRUE))
(HBA <- HausdorffDistance(B, A, TRUE))
arrows(0, 0.75, 1.15, 0.75, length = 0.1, code = 3)
text(0.5, 0.85, paste("H(A->B)", "=", HAB))
arrows(1, 0.15, 3, 0.15, length = 0.1, code = 3)
text(2, 0.25, paste("H(B->A)", "=", HBA))
text(1.5, -0.5, paste("H = max(H(A->B), H(B->A))", "=", H))
}
\keyword{manip}