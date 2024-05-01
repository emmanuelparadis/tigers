\name{rose}
\alias{rose}
\title{Draw a Compass}
\description{
  Draw a compass on a map or a plot.
}
\usage{
rose(x, y, size = 1, width = size/4, cols = c("grey10", "white"),
     labels = c("N", "S", "E", "W"), offset = 0, ...)
}
\arguments{
  \item{x, y}{coordinates of the center of the compass given as two
    vectors of length one (in user coordinates).}
  \item{size}{the length of the needles (in the same units than the user
    coordinates).}
  \item{width}{the width of the needles at the bottom; by default, one
    fourth of the previous.}
  \item{cols}{the colours for each side of the needles.}
  \item{labels}{the text printed at the tips of the needles.}
  \item{offset}{the space between the labels and the tips of the
    needles.}
  \item{\dots}{further arguments passed to \code{text} to format the
    labels.}
}
\details{
  To not print the directions, set \code{labels = rep("", 4)}.
}
\author{Emmanuel Paradis}
\keyword{aplot}
