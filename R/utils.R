## utils.R (2024-04-30)

##   Utilities for GIS data

## Copyright 2024 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../COPYING for licensing issues.

.check2cols <- function(x, y)
{
    if (is.null(y)) {
        ncx <- ncol(x)
        if (is.null(ncx) || ncx < 2)
            stop(paste(sQuote(deparse(substitute(x))),
                       "should have at least 2 columns if",
                       sQuote(deparse(substitute(y))), "is not given"),
                 call. = FALSE)
        y <- x[, 2L]
        x <- x[, 1L]
    } else {
        if (length(x) != length(y))
            stop(paste(sQuote(deparse(substitute(x))), "and",
                       sQuote(deparse(substitute(y))),
                       "should have the same length"), call. = FALSE)
    }
    list(x, y) # no names for faster computation
}

rose <- function(x, y, size = 1, width = size/4, cols = c("grey10", "white"),
                 labels = c("N", "S", "E", "W"), offset = 0, ...)
{
    w2 <- width / 2

    xx <- c(x, x, x + w2, NA, x, x + size, x + w2, NA,
            x, x, x - w2, NA, x, x - size, x - w2)
    yy <- c(y, y + size, y + w2, NA, y, y, y - w2, NA,
            y, y - size, y - w2, NA, y, y, y + w2)
    polygon(xx, yy, col = cols[2])
    xx <- c(x, x, x - w2, NA, x, x + size, x + w2, NA,
            x, x, x + w2, NA, x, x - size, x - w2)
    yy <- c(y, y + size, y + w2, NA, y, y, y + w2, NA,
            y, y - size, y - w2, NA, y, y, y - w2)
    polygon(xx, yy, col = cols[1])
    off <- strwidth("M")/2 + offset
    xx <- c(x, x, x + size + off, x - size - off)
    yy <- c(y + size + off, y - size - off, y, y)
    text(xx, yy, labels, ...)

}

axisMap <- function(latitude = FALSE, width = 0.05, len = 1, cols = c("black", "white"), ...)
{
    psr <- par("usr")
    at <- pretty(psr[if (latitude) 3:4 else 1:2])
    n <- length(at)
    if (latitude) {
        side <- 2
        hemisphere <- rep("N ", n)
        hemisphere[at < 0] <- "S "
        hemisphere[at == 0] <- ""
        labels <- paste0(hemisphere, abs(at), "\u00b0")
    } else {
        side <- 1
        half <- rep("E ", n)
        half[at < 0] <- "W "
        half[at == 0] <- ""
        labels <- paste0(half, abs(at), "\u00b0")
    }
    axis(side, at = at, labels = labels, las = 1,
         pos = psr[ifelse(latitude, 1L, 3L)], ...)
    selcol <- FALSE
    if (latitude) {
        w <- xinch(width)
        left <- c(psr[1] - w, psr[2])
        right <- c(psr[1], psr[2] + w)
        bottom <- psr[3]; top <- len * ceiling(ceiling(bottom)/len)
        while (bottom < psr[4]) {
            for (i in 1:2)
                rect(left[i], bottom, right[i], top,
                     col = cols[selcol + 1L], xpd = TRUE)
            bottom <- top
            top <- top + len
            if (top > psr[4]) top <- psr[4]
            selcol <- !selcol
        }
    } else {
        w <- yinch(width)
        top <- c(psr[3], psr[4] + w)
        bottom <- c(psr[3] - w, psr[4])
        left <- psr[1]; right <- len * ceiling(ceiling(left)/len)
        while (left < psr[2]) {
            rect(left, bottom, right, top,
                 col = cols[selcol + 1L], xpd = TRUE)
            left <- right
            right <- right + len
            if (right > psr[2]) right <- psr[2]
            selcol <- !selcol
        }
    }
}

## Source: http://www.physics.sfasu.edu/astro/color/spectra.html
wl2col <- function(x, gamma = 0.8, RGB = FALSE)
{
    out.of.range <- x < 380 | x > 780
    flag <- any(out.of.range)
    if (flag) {
        w <- which(!out.of.range)
        x <- x[w]
    }
    n <- as.integer(length(x))
    o <- .Fortran(wltocol, as.double(x), n, as.double(gamma), numeric(3 * n))
    res <- o[[4L]]
    dim(res) <- c(n, 3L)
    if (flag) {
        tmp <- matrix(0, length(out.of.range), 3L)
        tmp[w, ] <- res
        res <- tmp
    }
    if (RGB) {
        colnames(res) <- c("red", "green", "blue")
        if (!is.null(names(x))) rownames(res) <- names(x)
        return(res)
    }
    res <- rgb(res[, 1L], res[, 2L], res[, 3L])
    if (!is.null(names(x))) rownames(res) <- names(x)
    res
}

## Source: https://www.fourmilab.ch/documents/specrend/
spectrum2col <- function(spec, RGB = FALSE, no.warn = TRUE, color.system = 3)
{
    if (length(spec) != 81)
        stop("'spec' should have 81 values; see ?spectrum2col")
    o <- .C(C_specrend, spec, numeric(1), numeric(1), numeric(1), integer(1),
            as.integer(color.system))
    if (!no.warn) {
        if (o[[5]]) warning("computations were approximate")
    }
    if (RGB) return(c(red = o[[2]], green = o[[3]], blue = o[[4]]))
    rgb(o[[2]], o[[3]], o[[4]])
}

BlackBodySpectrum <- function(x, Temp = 300)
{
    wlm <- x * 1e-9 # nm -> m
    (3.74183e-16 / wlm^5) / (exp(1.4388e-2 / (wlm * Temp)) - 1)
}

fast2waytable <- function(x, y, levels = NULL)
{
    nulllevs <- is.null(levels)
    if (nulllevs) levels <- sort(unique(c(unique(x), unique(y))))
    K <- length(levels)
    transTab <- raw(256L)
    ## need to include one category for NA? maybe later...:
    transTab[levels + 1L] <- as.raw(0:(K - 1))
    res <- .Call(fast2waytable_Call, x, y, K, transTab)
    dimnames(res) <- list(levels, levels)
    w <- which(colSums(res) == 0)
    if (length(w)) res <- res[-w, -w, drop = FALSE]
    res
}
