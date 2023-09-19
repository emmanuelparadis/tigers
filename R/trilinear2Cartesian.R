## trilinear2Cartesian.R (2023-09-19)

##   Trilinear Coordinates

## Copyright 2023 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../DESCRIPTION for licensing issues.

trilinear2Cartesian <- function(p, X)
{
    if (nrow(X) != 3 || ncol(X) != 2)
        stop("X should be a matrix with 3 rows and 2 columns")
    A <- X[1L, ]; B <- X[2L, ]; C <- X[3L, ]
    ## sidelines:
    a <- sqrt((B[1] - C[1])^2 + (B[2] - C[2])^2)
    b <- sqrt((A[1] - C[1])^2 + (A[2] - C[2])^2)
    c <- sqrt((A[1] - B[1])^2 + (A[2] - B[2])^2)

    if (length(p) != 3) stop("p should have 3 values")
    ## scaling p:
    p <- p / sum(p)

    ax <- a * p[1L]; by <- b * p[2L]; cz <- c * p[3L]
    denom <- ax + by + cz
    rbind(A * ax / denom + B * by /denom + C * cz / denom)
}

Cartesian2trilinear <- function(xy, X)
{
    if (nrow(X) != 3 || ncol(X) != 2)
        stop("X should be a matrix with 3 rows and 2 columns")
    A <- X[1L, ]; B <- X[2L, ]; C <- X[3L, ]

    if (length(xy) != 2) stop("xy should have 2 values")

    ## orthogonal (unsigned) distance between the point xy and the
    ## line defined by the points P and Q (the order of the latter is
    ## unimportant):
    getDist <- function(xy, P, Q) {
        beta <- (Q[2] - P[2]) / (Q[1] - P[1])
        if (beta == 0) return(abs(xy[2] - P[2]))
        if (is.infinite(beta)) return(abs(xy[1] - P[1]))
        alpha <- P[2] - beta * P[1]
        abs((alpha + beta * xy[1] - xy[2]) / sqrt(1 + beta^2))
    }
    res <- c(getDist(xy, B, C), getDist(xy, A, C), getDist(xy, A, B))
    rbind(res / sum(res))
}
