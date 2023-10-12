## random_point_in_triangle.R (2023-09-19)

##   Simulate Random Points in a Triangle

## Copyright 2023 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../DESCRIPTION for licensing issues.

## https://www.cs.princeton.edu/~funk/tog02.pdf
## <doi:10.1145/571647.571648>

rpit <- function(n, X, rfun1 = runif, rfun2 = runif)
{
    if (nrow(X) != 3 || ncol(X) != 2)
        stop("X should be a numeric matrix with 3 rows and 2 columns")
    storage.mode(X) <- "double"

    sr1 <- sqrt(rfun1(n))
    r2 <- rfun2(n)

    A <- 1 - sr1
    B <- sr1 * (1 - r2)
    C <- r2 * sr1

    xy <- matrix(NA_real_, n, 2L)
    for (j in 1:2)
        xy[, j] <- A * X[1L, j] + B * X[2L, j] + C * X[3L, j]
    xy
}

random_point_in_triangle <- rpit
