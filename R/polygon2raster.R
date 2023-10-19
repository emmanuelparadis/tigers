## polygon2raster.R (2023-10-19)

##   Polygon Rasterisation

## Copyright 2023 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../DESCRIPTION for licensing issues.

polygon2mask <- function(XY, extent = NULL, k = 360,
                         value = 1L, backgrd = 0L)
{
    if (is.null(extent)) {
        west <- floor(min(XY[, 1L]))
        east <- ceiling(max(XY[, 1L]))
        south <- floor(min(XY[, 2L]))
        north <- ceiling(max(XY[, 2L]))
    } else {
        west <- extent[1L]
        east <- extent[2L]
        south <- extent[3L]
        north <- extent[4L]
    }
    NC <- (east - west) * k
    NR <- (north - south) * k
    z <- integer(NR * NC)
    if (!identical(backgrd, 0L)) z[] <- as.integer(backgrd)
    if (is.open(XY)) XY <- .close(XY)
    ## transform the coordinates (see C code for details):
    XY[, 1L] <- XY[, 1L] - west
    XY[, 2L] <- north - XY[, 2L]
    XY <- round(XY*k + 0.5/k)
    XY <- redundantVertices(XY, tol = 0)
    ## get the horizontal limits of polygon to avoid scanning
    ## through all the raster:
    north.lim <- floor(min(XY[, 2L]))
    south.lim <- ceiling(max(XY[, 2L]))
    if (north.lim == south.lim) {
        ik <- unique(XY[, 1L] + NC * XY[, 2L] + 1L)
        z[ik] <- value
        return(z)
    }
    PARS <- as.integer(c(NC, north.lim, south.lim, value))
    o <- .Call(singlePolygon2raster, XY, PARS, z)
    dim(z) <- c(NC, NR)
    z
    ##matrix(z, NR, NC, TRUE)
}
