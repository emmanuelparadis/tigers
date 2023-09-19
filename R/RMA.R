## RMA.R (2023-09-19)

##   Reduced Major Axis of a Polygon

## Copyright 2023 Emmanuel Paradis

## This file is part of the R-package `tigers'.
## See the file ../DESCRIPTION for licensing issues.

## https://mathworld.wolfram.com/LeastSquaresFittingPerpendicularOffsets.html

RMA <- function(x, y = NULL)
{
    res <- .Call(RMA_Call, x, y)
    colnames(res) <- c("alpha", "beta")
    res
}

#RMA <- function(x, y = NULL)
#{
#    if (is.null(y)) {
#        y <- x[, 2L]
#        x <- x[, 1L]
#    }
#    mx <- mean(x)
#    my <- mean(y)
#    vx <- var(x)
#    vy <- var(y)
#    cxy <- cov(x, y)
#    if (!cxy) {
#        beta <- c(0, Inf)
#        alpha <- c(my, mx)
#    } else {
#        B <- -0.5 * (vy - vx ) / cxy
#        beta <- -B + c(-1, 1) * sqrt(B^2 + 1)
#        alpha <- my - beta * mx
#    }
#    cbind(alpha = alpha, beta = beta)
#}
