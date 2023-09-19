## FUNs.R (2023-09-19)

##   Various Functions

## Copyright 2023 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../DESCRIPTION for licensing issues.

## P, Q: two-column matrices with coordinates of the vertices
## returns TRUE if the bounding boxes of two polygons overlap
## (returns FALSE if the boxes are contiguous)
.overlappingBbox <- function(P, Q)
{
    tmp <- P[, 1L]; xminP <- min(tmp); xmaxP <- max(tmp)
    tmp <- Q[, 1L]; xminQ <- min(tmp); xmaxQ <- max(tmp)
    if (xminP >= xmaxQ) return(FALSE)
    if (xminQ >= xmaxP) return(FALSE)
    tmp <- P[, 2L]; yminP <- min(tmp); ymaxP <- max(tmp)
    tmp <- Q[, 2L]; yminQ <- min(tmp); ymaxQ <- max(tmp)
    if (yminP >= ymaxQ) return(FALSE)
    if (yminQ >= ymaxP) return(FALSE)
    TRUE
}

## P, Q: two-column matrices with coordinates of the vertices
haveOverlap <- function(A, B)
{
    if (nrow(A) < 3 || nrow(B) < 3) {
        warning("at least one polygon has less than 3 vertices: returning NULL")
        return(NULL)
    }
    if (!.overlappingBbox(A, B)) return(FALSE)
    if (samePolygons(A, B)) return(TRUE)
    if (.Call(haveOverlapTwoPolygons, A, B)) return(TRUE)
    .Call(haveOverlapTwoPolygons, B, A)
}

## x: a two-column matrix
is.clockwise <- function(x)
### https://en.wikipedia.org/wiki/Curve_orientation
{
    ## if (is.null(y)) {
    y <- x[, 2L]
    x <- x[, 1L]
    ## }
    n <- length(x)
    ## make sure the polygon is open:
    if (x[1L] == x[n] && y[1L] == y[n]) {
        x <- x[-n]
        y <- y[-n]
        n <- n - 1L
    }
    b <- which.min(x)
    a <- b - 1L
    if (a < 1) a <- n
    c <- b + 1L
    if (c > n) c <- 1L
    det <- (x[b] - x[a])*(y[c] - y[a]) - (x[c] - x[a])*(y[b] - y[a])
    det < 0
}

is.convex <- function(x)
{
    n <- nrow(x)
    if (n < 4) return(TRUE)

    test <- function() sign((X[c] - X[a]) * (Y[b] - Y[c]) - (X[c] - X[b]) * (Y[a] - Y[c]))

    X <- x[, 1L]; Y <- x[, 2L]
    a <- 1L; b <- 2L; c <- 3L
    init <- test()

    ## in case aligned vertices
    while (!init && c < n) {
        a <- b
        b <- c
        c <- c + 1L
        init <- test()
    }

    while (c < n) {
        a <- b
        b <- c
        c <- c + 1L
        if (test() != init) return(FALSE)
    }
    ## the last 2 triplets:
    a <- n - 1L; b <- n; c <- 1L
    if (test() != init) return(FALSE)
    a <- n; b <- 1L; c <- 2L
    if (test() != init) return(FALSE)
    TRUE
}

revPolygon <- function(x, copy = TRUE)
{
    copy <- as.integer(copy)
    if (storage.mode(x) != "double")
        stop("values in 'x' should be stored as double")
    dx <- dim(x)
    if (is.null(dx)) {
        res <- .Call(rev_Call, x, copy)
    } else {
        if (length(dx) != 2 || dx[2] != 2)
            stop("'x' should be either a vector or a matrix with 2 columns")
        res <- .Call(rev_2cols_Call, x, copy)
    }
    if (copy) res else invisible(res)
}

## A, B: two-column matrices with coordinates of the vertices
##       (the polygons may be "closed" or "open")
## digits: precision of the coordinates
samePolygons <- function(A, B, digits = 10)
{
    ## make sure the polygons are open before comparing them:
    n <- nrow(A)
    if (A[1L, 1L] == A[n, 1L] && A[1L, 2L] == A[n, 2L]) A <- A[-n, ]
    n <- nrow(B)
    if (B[1L, 1L] == B[n, 1L] && B[1L, 2L] == B[n, 2L]) {
        B <- B[-n, ]
        n <- n - 1L
    }
    if (nrow(A) != n) return(FALSE)
    X <- paste(round(A[, 1L], digits = digits), round(A[, 2L], digits = digits))
    Y <- paste(round(B[, 1L], digits = digits), round(B[, 2L], digits = digits))
    if (!is.clockwise(A)) X <- rev(X)
    if (!is.clockwise(B)) Y <- rev(Y)
    m <- match(X, Y)
    if (anyNA(m)) return(FALSE)
    dm <- diff(m)
    sdm <- sum(dm == 1)
    if (sdm == n - 1) return(TRUE)
    if (sdm == n - 2 && sum(dm == -(n - 1)) == 1) return(TRUE)
    FALSE
}

redundantVertices <- function(x, tol = 1e-8, check.only = FALSE)
{
    res <- .Call(redundant_vertices, x, tol, check.only)
    if (check.only) invisible(res) else res
}

chullPolygon <- function(x, y = NULL)
{
    if (!is.null(y)) x <- cbind(x, y)
    cwo <- is.clockwise(x)
    if (!cwo) {
        n <- nrow(x)
        x <- x[n:1, ]
    }
    ## NOTE: the 2nd argument is unused
    res <- .Call(convex_hull_C, x, 1L) + 1L
    if (!cwo) n - res + 1L else res
}

is.insidePolygon <- function(XY, points)
{
    if (!is.open(XY)) XY <- .open(XY)
    .Call(InsidePolygon_Call, XY, points)
}

barycoords <- function(XY, point)
{
    n <- nrow(XY)
    if (length(point) != 2)
        stop("'point' must have 2 values")
    BC <- numeric(n)
    res <- .C(test_Hormann_Floater, XY, point, BC, n)
    res[[3]]
}

area <- function(x, y = NULL, lonlat = FALSE, spheroid = TRUE)
{
    if (!is.null(y)) x <- cbind(x, y)
    if (lonlat) x <- x * pi/180
    .Call(area_Call, x, as.integer(c(lonlat, spheroid)))
}
