## geo.R (2023-09-19)

##   Tools for Geographic Data

## Copyright 2023 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../DESCRIPTION for licensing issues.

geoTrans <- function(x, degsym = NULL, minsym = "'", secsym = "\"")
{
    if (is.null(degsym)) degsym <- "\u00b0"
    x <- as.character(x)
    sow <- grep("[SOW]", x)
    x <- gsub("[SNOWE]", "", x)
    x <- strsplit(x, paste0("[", degsym, minsym, secsym, "]"))
    foo <- function(x) {
        x <- as.numeric(x)
        sum(x[1:3] / c(1, 60, 3600), na.rm = TRUE)
    }
    res <- sapply(x, foo)
    if (length(sow)) res[sow] <- -res[sow]
    res
}

geoTrans2 <- function(lon, lat = NULL, degsym = NULL, minsym = "'",
                      secsym = "\"", dropzero = FALSE, digits = 3, latex = FALSE)
{
    if (is.null(lat)) {
        nc <- ncol(lon)
        if (is.null(nc) || nc < 2)
            stop("if 'lat' is not given, 'lon' should have at least 2 columns")
        lat <- lon[, 2]
        lon <- lon[, 1]
    } else {
        if (length(lon) != length(lat))
            stop("'lat' and 'lon' should have the same length")
    }
    if (is.null(degsym)) degsym <- "\u00b0"
    foo <- function(x) {
        d <- floor(round(x, 8))
        m2 <- (x - d) * 60
        m <- floor(round(m2, 8))
        s <- round(60 * (m2 - m), digits = digits)
        if (dropzero) {
            m <- ifelse(m == 0 & s == 0, "", paste0(m, minsym))
            s <- ifelse(s == 0, "", paste0(s, secsym))
        } else {
            m <- paste0(m, minsym)
            s <- paste0(s, secsym)
        }
        d <- paste0(d, degsym)
        paste0(d, m, s)
    }
    x <- foo(lon)
    y <- foo(lat)
    res <- paste0(ifelse(lat < 0, "S ", "N "), y, ", ", ifelse(lon < 0, "W ", "E "), x)
    if (latex) {
        res <- gsub(degsym, "\\\\textdegree ", res)
        res <- gsub("'", "$'$", res)
        res <- gsub("\"", "$''$", res)
    }
    res
}

geod <- function(lon, lat = NULL, R = 6371)
{
    deg2rad <- function(x) x * pi / 180
    if (is.null(lat)) {
        lat <- lon[, 2L]
        lon <- lon[, 1L]
    }
    n <- length(lat)
    lat <- deg2rad(lat)
    lon <- deg2rad(lon)
    absdiff <- function(x, y) abs(x - y)
    Delta_lon <- outer(lon, lon, absdiff)
    Delta_lat <- outer(lat, lat, absdiff)
    tmp <- cos(lat) # store the cosinus(lat) for all individuals
    A <- (sin(Delta_lat / 2))^2 + rep(tmp, n) * rep(tmp, each = n) * (sin(Delta_lon / 2))^2
    R * 2 * asin(sqrt(A))
### alternative:
### tmp <- sin(lat); tmp2 <- cos(lat)
### A <- acos(outer(tmp, tmp) + outer(tmp2, tmp2) * cos(Delta_lon))
### R * A
}

gcl <- function(x0, y0, x1, y1, linear = FALSE, npoints = 100)
{
    x0 <- x0 * pi/180
    y0 <- y0 * pi/180
    x1 <- x1 * pi/180
    y1 <- y1 * pi/180
    if (npoints < 2) stop("'npoints' should be at least 3")
    x <- seq(x0, x1, length.out = npoints)
    if (linear) {
        ## linear interpolation (Chamberlain and Duquette 2007, p. 9)
        foo <- function(X) y0 + (X - x0) * (y1 - y0) / (x1 - x0)
    } else {
        ## great circle line (Chamberlain and Duquette 2007, p. 10)
        foo <- function(X) {
            denom <- sin(x0 - x1)
            A <- tan(y0)
            B <- tan(y1)
            atan((A * sin(X - x1) - B * sin(X - x0)) / denom)
        }
    }
    cbind(x = x, y = foo(x)) * 180/pi
}

great_circle_line <- gcl

## https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Another_formula
## https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points (no need to calculate the slopes and intercepts)
dtl <- function(x, y = NULL, x0, y0, x1, y1, alpha = NULL, beta = NULL)
{
    if (!is.null(y)) x <- cbind(x, y)
    missx0 <- missing(x0)
    missy0 <- missing(y0)
    missx1 <- missing(x1)
    missy1 <- missing(y1)
    miss <- c(missx0, missy0, missx1, missy1)
    if (any(miss)) {
        if (is.null(alpha) || is.null(beta))
            stop("alpha and beta must be given if at least one of x0, y0, x1, and y1 is missing")
        return((alpha + beta * x[, 1L] - x[, 2L]) / sqrt(1 + beta^2))
    }
    if (sum(miss))
        stop("at least one of x0, y0, x1, and y1 is missing")
    Delta.x <- x1 - x0
    Delta.y <- y1 - y0
    num <- Delta.x * (y0 - x[, 2L]) - (x0 - x[, 1L]) * Delta.y
    denom <- sqrt(Delta.x^2 + Delta.y^2)
    num / denom
}

distance_to_line <- dtl

dta <- function(x, y = NULL, x0, y0, x1, y1, prec = 0.001)
{
    if (!is.null(y)) x <- cbind(x, y)
    foo <- function(x) {
        repeat {
            ## 1. trace an arc between the 2 points (0,1)
            XY <- gcl(x0, y0, x1, y1, npoints = 10L)
            ## 2. distances between the focus point and all points
            ## on the arc
            d <- geod(rbind(x, XY))[-1L, 1L]
            ## 3. find the shortest distance
            s <- which.min(d) # should not be 1 neither 10
            if (s == 1 || s == 10L) {
                if (!d[s]) return(0)
                stop("something wrong happened")
            }
            xys <- XY[s + c(-1, 1), ]
            ## 4. if the distance between the 2 points before and
            ## after the point found at 3. < precision, stop
            if (geod(xys)[2L] < prec) break
            ## 5. update the coordinates for the next iteration
            x0 <- xys[1L, 1L]; y0 <- xys[1L, 2L]
            x1 <- xys[2L, 1L]; y1 <- xys[2L, 2L]
        }
        d[s]
    }
    apply(x, 1, foo)
}

distance_to_arc <- dta
