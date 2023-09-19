## HausdorffDistance.R (2023-09-19)

##   Hausdorff Distance Between Two Polygons

## Copyright 2023 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../DESCRIPTION for licensing issues.

HausdorffDistance <- function(A, B, directed = FALSE)
{
    s <- 1:nrow(A)
    D <- as.matrix(dist(rbind(A, B)))
    D <- D[s, -s]
    h <- max(apply(D, 1, min))
    if (directed) return(h)
    max(h, apply(D, 2, min))
}
