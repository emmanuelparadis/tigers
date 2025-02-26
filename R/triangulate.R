## triangulate.R (2025-02-26)

##   Polygon Triangulation

## Copyright 2023-2025 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../DESCRIPTION for licensing issues.

is.open <- function(x) any(x[1L, ] != x[nrow(x), ])

.open <- function(XY) XY[-nrow(XY), , drop = FALSE]

.close <- function(XY) rbind(XY, XY[1L, , drop = FALSE])

triangulate <- function(x, y = NULL, method = 1)
{
    method <- as.integer(method)
    if (method < 1 || method > 2)
        stop("method should be 1 or 2")
    if (!is.null(y)) x <- cbind(x, y)
    o <- is.open(x)
    cw <- is.clockwise(x)
    wasClosed <- wasReversed <- FALSE # two indicators
    if (method == 1) {
        if (o) {
            x <- .close(x)
            wasClosed <- TRUE
            N <- nrow(x) # one row was added
        }
    } else {
        if (!o) x <- .open(x)
    }
    if (method == 2 && cw) {
        x <- revPolygon(x)
        wasReversed <- TRUE
        n <- nrow(x) - as.integer(wasClosed)
    }
    res <- .Call(triangulate_Call, x, method) + 1L
    if (wasClosed) res[res == N] <- 1L
    if (wasReversed) res[] <- (n:1)[as.vector(res)]
    res
}
