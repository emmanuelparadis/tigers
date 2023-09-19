## convexPolygonOverlap.R (2023-09-19)

##   Polygon Decomposition and Overlap

## Copyright 2023 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../DESCRIPTION for licensing issues.

convexPolygonOverlap <- function(A, B)
{
    .Call(convexPolygonOverlap_Call, A, B)
}

decomposePolygon <- function(x, y = NULL, method = 1, quiet = FALSE)
{
    if (!is.null(y)) x <- cbind(x, y)
    if (is.convex(x)) {
        if (!quiet) warning("polygon is convex: returning NULL")
        return(NULL)
    }
    method <- as.integer(method)
    if (is.clockwise(x)) x <- revPolygon(x)
    tri <- triangulate(x, method = method)

    n <- nrow(x)
    nm1 <- n - 1L

    ## extract the diagonals implied by the triangulation
    diags <- matrix(NA_integer_, 0, 2)
    for (k in 1:nrow(tri)) {
        z <- tri[k, ]
        for (i in 1:2) {
            for (j in (i + 1):3) {
                tmp <- abs(z[i] - z[j])
                if (tmp != 1 && tmp != nm1)
                    diags <- rbind(diags, sort(z[c(i, j)]))
            }
        }
    }
    diags <- unique(diags)

    prv <- c(n, 1:nm1)
    iii <- 1:n
    nxt <- c(2:n, 1)

    ## returns 1, 0, or - 1, depending on whether c is to the right
    ## of, collinear with, or to the left of L[a, b]
    angle <- function(a, b, c)
        sign((x[c, 1L] - x[a, 1L]) * (x[b, 2L] - x[c, 2L]) -
             (x[c, 1L] - x[b, 1L]) * (x[a, 2L] - x[c, 2L]))

    beta <- mapply(angle, a = prv, b = iii, c = nxt)
    del <- beta[diags[, 1L]] + beta[diags[, 2L]] == -2L
    if (any(del)) diags <- diags[!del, , drop = FALSE]
    diags
}

polygonOverlap <- function(A, B)
{
    ## 1. get the polygons from the decomposition in convex polygons
    foo <- function(XY) {
        decomp <- decomposePolygon(XY, quiet = TRUE)
        if (is.null(decomp)) return(list(XY))
        vrtx <- seq_len(nrow(XY))
        N <- nrow(decomp)
        res <- list(vrtx)
        for (i in 1:N) {
            m <- lapply(res, function(x) match(decomp[i, ], x))
            w <- which(sapply(m, function(x) !anyNA(x)))
            path <- res[[w]]
            tmp <- m[[w]]
            a <- tmp[1L]
            b <- tmp[2L]
            n <- length(path)
            res[[w]] <- list(path[a:b], path[c(b:n, 1:a)])
            res <- unlist(res, recursive = FALSE)
        }
        lapply(res, function(i) XY[i, ])
    }
    D1 <- foo(A)
    D2 <- foo(B)
    n1 <- length(D1)
    n2 <- length(D2)

    ## 2. get the intersection between both sets of convex polygons
    for (i in 1:n1)
        if (!is.clockwise(D1[[i]]))
            revPolygon(D1[[i]], FALSE)
    for (i in 1:n2)
        if (!is.clockwise(D2[[i]]))
            revPolygon(D2[[i]], FALSE)
    POLY <- vector("list", n1 * n2)
    k <- 0L
    for (i in 1:n1) {
        for (j in 1:n2) {
            k <- k + 1L
            POLY[[k]] <- convexPolygonOverlap((D1[[i]]), (D2[[j]]))
        }
    }
    POLY
}
