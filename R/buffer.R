## buffer.R (2024-07-11)

##   Buffer Around Polygons

## Copyright 2024 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../COPYING for licensing issues.

.polar2rect  <- function (r, angle)
    list(x = r * cos(angle), y = r * sin(angle))

buffer <- function(x, r, smoothing = 1, unit = "degree", quiet = FALSE)
{
    unit <- match.arg(unit, c("degree", "km", "m"))
    inc <- switch(unit, "m" = 1/r, "km" = 1/(r*1000), "degree" = 0.01)
    inc <- inc * smoothing

    if (!is.clockwise(x)) x <- revPolygon(x)
    if (is.open(x)) x <- .close(x)
    n <- nrow(x) # assumed closed
    x0 <- x[-n, 1L]
    y0 <- x[-n, 2L]
    x1 <- x[-1L, 1L]
    y1 <- x[-1L, 2L]

    theta <- atan2(y1 - y0, x1 - x0) + pi/2
    A <- B <- .polar2rect(r, theta)
    B$x <- B$x + x0
    B$y <- B$y + y0
    A$x <- A$x + x1
    A$y <- A$y + y1

    from_to_loop <- function(i, j, n) {
        if (j == n + 1L) return(1:i)
        if (j < i) return(j:i)
        c(j:n, 1:i)
    }

    updateDF <- function(a, b, c, d) {
        i <- length(DF[[1L]]) + 1L
        DF[[1L]][i] <<- a
        DF[[2L]][i] <<- b
        DF[[3L]][i] <<- c
        DF[[4L]][i] <<- d
    }

    smooth <- function(x0, y0, x1, y1, xr, yr) {
        th0 <- atan2(y0 - yr, x0 - xr)
        th1 <- atan2(y1 - yr, x1 - xr)
        if (th0 == th1) return(invisible())
        if (th0 < 0) th0 <- th0 + 2 * pi
        if (th1 < 0) th1 <- th1 + 2 * pi
        delta.th <- if (th0 > th1) th0 - th1 else th0 + (2 * pi - th1)
        if (abs(delta.th - pi*2) < 1e-4) return(invisible())
        n <- delta.th %/% inc
        rest <- delta.th %% inc
        if (rest > inc / 10) n <- n - 1
        if (n < 1) return(invisible())
        if (!rest) inc <- delta.th / n
        angles <- th0 - inc * 1:(n - 1)
        angles <- angles %% (2 * pi)
        s <- which(angles > pi)
        if (length(s)) angles[s] <- angles[s] - 2 * pi
        res <- .polar2rect(r, angles)
        res$x <- res$x + xr
        res$y <- res$y + yr
        as.vector(rbind(res$x, res$y))
    }

    n <- n - 1L # closed -> open
    ## start from the westmost vertex
    start <- which.min(B$x)
    ARC <- FALSE
    if (A$x[i <- which.min(A$x)] < B$x[start]) {
        start <- i
        ARC <- !ARC
    }
    stop <- start
    if (stop == -1) stop <- n

    ## indicator which is switched to TRUE after the origin of the
    ## polygon 'x' has been reached
    READY <- FALSE
    testBREAK <- function() j >= stop && READY
    if (start == 1)
        testBREAK <- function() j == stop && length(FINAL) > 4

    DF0 <- list(index = integer(), x = numeric(),
                y = numeric(), arcORedge = logical())

    i <- start
    FINAL <- if (ARC) c(A$x[i], A$y[i]) else c(B$x[i], B$y[i])
    x0 <- FINAL[1]
    y0 <- FINAL[2]
    ## transpose FINAL at the end
    if (!quiet) {
        cat(n, "vertices\n")
        cat(start, ".........", n, " 1..........", stop, "\n|", sep = "")
    }
    block.size <- ceiling(n / 37)
    Nblocks <- 0L
    Nb2 <- 0L

    repeat {
        if (!quiet) {
            if (i > start && (i - start) %/% block.size > Nblocks) {
                cat("=")
                Nblocks <- Nblocks + 1L
            }
            if (i < start && i %/% block.size > Nb2) {
                cat("=")
                Nb2 <- Nb2 + 1L
            }
        }
        j <- i + as.integer(ARC)
        if (j == n + 1L) j <- 1L
###
        if (ARC) {
            x1 <- B$x[j]
            y1 <- B$y[j]
        } else {
            x1 <- A$x[j]
            y1 <- A$y[j]
        }

if (abs(x0 - x1) >= 1e-8 || abs(y0 - y1) >= 1e-8) {
###
        SI <- c(x0, y0, x1, y1)
        DF <- DF0
        if (ARC) {
            arc <- c(x[i + 1L, ], r, SI)
            ## arc-arc comparisons
            K <- from_to_loop(stop, j + 1L, n)
            for (k in K) {
                k2 <- k + 1L
                if (k2 == n + 1L) k2 <- 1L
                arc2 <- c(x[k2, ], r, A$x[k], A$y[k], B$x[k2], B$y[k2])
                o <- intersArcArc(arc, arc2)
                if (o[[4]]) {
                    updateDF(k, o[[3]][1], o[[3]][2], TRUE)
                    if (o[[4]] > 1)
                        updateDF(k, o[[3]][3], o[[3]][4], TRUE)
                }
            }
            ## arc-edge comparisons
            K <- from_to_loop(stop, j + 1L, n)
            for (k in K) {
                seg <- c(B$x[k], B$y[k], A$x[k], A$y[k])
                o <- .Call("inters_seg_arc", seg, arc)
                if (o[[2]]) {
                    updateDF(k, o[[1]][1], o[[1]][2], FALSE)
                if (o[[2]] > 1)
                    updateDF(k, o[[1]][3], o[[1]][4], FALSE)
                }
            }
        } else {
            ## edge-edge comparisons
            K <- from_to_loop(stop, j + 1L, n)
            for (k in K) {
                SII <- c(B$x[k], B$y[k], A$x[k], A$y[k])
                o <- .Call("inters_seg_seg", SI, SII)
                if (o[[2]]) updateDF(k, o[[1]][1], o[[1]][2], FALSE)
            }
            ## edge-arc comparisons
            K <- from_to_loop(stop, j, n)
            for (k in K) {
                k2 <- k + 1L
                if (k2 == n + 1L) k2 <- 1L
                arc <- c(x[k2, ], r, A$x[k], A$y[k], B$x[k2], B$y[k2])
                o <- .Call("inters_seg_arc", SI, arc)
                if (o[[2]]) {
                    updateDF(k, o[[1]][1], o[[1]][2], TRUE)
                    if (o[[2]] > 1)
                        updateDF(k, o[[1]][3], o[[1]][4], TRUE)
                }
            }
        }
###
        if (nINTS <- length(DF$index)) {
            ## find the nearest intersection
            if (nINTS > 1) {
                ##repeat {
                if (ARC) {
                    ##w <- which.min(atan2(x[i, 2L] - x[i + 1L, 2L], x[i, 1L] - x[i + 1L, 1L]) - atan2(DF$y - x[i + 1L, 2L], DF$x - x[i + 1L, 1L]))
                    ## fixed 2024-07-11:
                    w <- which.min(atan2(x[i, 2L] - y0, x[i, 1L] - x0) - atan2(DF$y - y0, DF$x - x0))
                } else {
                    w <- which.min((x0 - DF$x)^2 + (y0 - DF$y)^2)
                }
                ##    if (DF$index[w] < i) DF <- lapply(DF, "[", -w) else break
                ##}
                DF <- lapply(DF, "[", w) # just get the closest point
            }
            x1 <- DF$x
            y1 <- DF$y
        }
###
        if (ARC) {
            o <- smooth(x0, y0, x1, y1, x[i + 1L, 1L], x[i + 1L, 2L])
            FINAL <- c(FINAL, o)
        }
###
        FINAL <- c(FINAL, x1, y1)
} else {
    nINTS <- 0L
}
        ## the algo is not supposed to move backward:
        if (nINTS && DF$index < i) nINTS <- 0L
        ## update i and ARC (even if no vertex was added)
        if (nINTS) {
            ARC <- DF$arcORedge
            i <- DF$index
        } else {
            i <- j
            if (i == n + 1L) i <- 1L
            ARC <- !ARC
        }
###
        x0 <- x1
        y0 <- y1
        if (testBREAK()) break
        if (j < stop) READY <- TRUE
    }
    if (!quiet) cat("|\n")
    ## return the final matrix of vertices after transposing
    matrix(FINAL, ncol = 2, byrow = TRUE)
}

## arc[1|2] is a vector with: the x and y of the center, and the radius
## this version checks intersections of 2 circles -- not arcs!
circle2circle <- function(arc1, arc2, quiet = FALSE)
{
    Dx <- arc2[1] - arc1[1]
    Dy <- arc2[2] - arc1[2]
    if (!Dx && !Dy) {
        if (!quiet) warning("the circles have the same center")
        return(NULL)
    }
    d <- sqrt(Dx^2 + Dy^2)
    if (d >= arc1[3] + arc2[3]) {
        if (!quiet) warning("no overlap of the circles")
        return(NULL)
    }
    l <- (arc1[3]^2 - arc2[3]^2 + d^2)/(2 * d)
    if (l < 0 || !is.finite(l)) {
        if (!quiet) warning("no overlap of the circles")
        return(NULL)
    }
    h <- sqrt(arc1[3]^2 - l^2)
    q1 <- l / d
    q2 <- h / d
    x <- Dx * q1 + c(1, -1) * Dy * q2 + arc1[1]
    y <- Dy * q1 + c(-1, 1) * Dx * q2 + arc1[2]
    cbind(x, y)
}

## see comments in C code
betweenAngles <- function(x, B, A)
{
    if (x < 0) x <- x + 2 * pi
    if (A < 0) A <- A + 2 * pi
    if (B < 0) B <- B + 2 * pi
    if (B > A) {
	if (x < B && x > A) return(1L)
    } else {
	## if B < A, then the arc overlaps with 0
	if (x < B) return(1L)
	if (x > A) return(1L)
    }
    0L
}

intersArcArc <- function(arc1, arc2)
{
    res <- list(arc1, arc2, numeric(4), 0L)
    o <- circle2circle(arc1[1:3], arc2[1:3], TRUE)
    if (is.null(o)) return(res)
    ## get the angles of both pairs of points on each arc:
    ## (note A->B is in CWO)
    th1A <- atan2(arc1[5] - arc1[2], arc1[4] - arc1[1])
    th1B <- atan2(arc1[7] - arc1[2], arc1[6] - arc1[1])
    th2A <- atan2(arc2[5] - arc2[2], arc2[4] - arc2[1])
    th2B <- atan2(arc2[7] - arc2[2], arc2[6] - arc2[1])
    x1 <- atan2(o[1, 2] - arc1[2], o[1, 1] - arc1[1])
    x2 <- atan2(o[1, 2] - arc2[2], o[1, 1] - arc2[1])
    ## 1) is the 1st point on both arcs?
    if (betweenAngles(x1, th1A, th1B) &&
        betweenAngles(x2, th2A, th2B)) {
        res[[3L]][1:2] <- o[1, ]
        res[[4L]] <- res[[4L]] + 1L
    }
    ## 2) id. for the 2nd point?
    x1 <- atan2(o[2, 2] - arc1[2], o[2, 1] - arc1[1])
    x2 <- atan2(o[2, 2] - arc2[2], o[2, 1] - arc2[1])
    if (betweenAngles(x1, th1A, th1B) &&
        betweenAngles(x2, th2A, th2B)) {
        i <- if (res[[4L]] == 1) 3:4 else 1:2
        res[[3L]][i] <- o[2, ]
        res[[4L]] <- res[[4L]] + 1L
    }
    res
}
