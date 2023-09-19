## triangulate.R (2023-09-19)

##   Polygon Triangulation

## Copyright 2023 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../DESCRIPTION for licensing issues.

is.open <- function(x) any(x[1L, ] != x[nrow(x), ])

.open <- function(XY) XY[-nrow(XY), , drop = FALSE]

.close <- function(XY) rbind(XY, XY[1L, , drop = FALSE])

triangulate <- function(x, y = NULL, method = 1)
{
    method <- as.integer(method)
    if (method > 4)
        stop("method should be an integer between 1 and 4")
    if (!is.null(y)) x <- cbind(x, y)
    o <- is.open(x)
    cw <- is.clockwise(x)
    if (method == 1) {
        if (o) x <- .close(x)
    } else {
        if (!o) x <- .open(x)
    }
    if ((method == 2 && cw) || (method == 3 && !cw))
        x <- revPolygon(x)
    .Call(triangulate_Call, x, method) + 1L
}
