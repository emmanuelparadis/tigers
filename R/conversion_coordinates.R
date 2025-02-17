## conversion_coordinates.R (2024-04-30)

##   Conversion Coordinates

## Copyright 2024 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../COPYING for licensing issues.

lonlat2UTM <- function(lon, lat = NULL, details = FALSE)
{
    DF <- .check2cols(lon, lat)
    res <- .Call("lonlat2UTM_Call", DF[[1]], DF[[2]])
    colnames(res) <- c("Easting", "Northing", "Zone")

    if (details) {
        res <- as.data.frame(res)
        res$Prefix <- paste0(res[, 3], c("S", "N")[as.integer(DF[[2]] > 0) + 1L])
    } else {
        res <- res[, -3]
    }
    res
}

UTM2lonlat <- function(x, y = NULL, zone = NULL, hemisphere = "N")
{
    if (is.data.frame(x) && length(x) == 4L) {
        y <- x[[2]]
        zone <- x[[3]]
        hemi <- ifelse(grepl("N", x[[4]]), 1L, -1L)
        x <- x[[1]] # en dernier ;)
    } else {
        xy <- .check2cols(x, y)
        x <- xy[[1]]
        y <- xy[[2]]
        n <- length(x)
        zone <- as.integer(zone)
        if (length(zone) < n) zone <- rep(zone, length.out = n)
        hemi <- setNames(c(1L, -1L), c("N", "S"))[hemisphere]
        if (length(hemi) < n) hemi <- rep(hemi, length.out = n)
    }
    .Call("UTM2lonlat_Call", x, y, zone, hemi)
}

lonlat2ECEF <- function(lon, lat = NULL, alt = 0, as.matrix = TRUE)
{
    xy <- .check2cols(lon, lat)
    lon <- xy[[1L]]
    lat <- xy[[2L]]
    a <- 6378137 # equatorial radius (m)
    b <- 6356752 # polar radius (m)
    b2a2 <- b^2/a^2
    e2 <- 1 - b2a2
    deg2rad <- pi/180
    lon <- lon * deg2rad
    lat <- lat * deg2rad
    clat <- cos(lat)
    slat <- sin(lat)
    clon <- cos(lon)
    slon <- sin(lon)

    N <- a / sqrt(1 - e2 * slat^2)

    x <- (N + alt) * clat * clon
    y <- (N + alt) * clat * slon
    z <- (b2a2 * N + alt) * slat
    if (as.matrix) cbind(x = x, y = y, z = z) else list(x = x, y = y, z = z)
}

ECEF2lonlat <- function(x, y = NULL, z = NULL)
{
    if (!is.null(y) && !is.null(z))
        x <- cbind(x, y, z)
    if (ncol(x) < 3) stop("'x' should have at least 3 columns")
    res <- .Call(ECEF2lonlat_Call, x)
    colnames(res) <- c("Lon", "Lat", "Alt")
    res
}
