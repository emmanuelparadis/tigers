\name{polygon2mask}
\alias{polygon2mask}
\alias{polygon2raster}
\title{Convert Polygon to a Raster Mask}
\description{
  Takes a polygon and returns a matrix with a mask that can be input
  into a raster.
}
\usage{
polygon2mask(XY, extent = NULL, k = 360,
             value = 1L, backgrd = 0L)
}
\arguments{
  \item{XY}{A two-column matrix giving the coordinates of a polygon.}
  \item{extent}{a vector with four numeric values giving the extent of
    the raster. By default, values are determined to minimally cover the
    polygon.}
  \item{k}{an integer value giving the number of pixels per unit (i.e.,
    the inverse of the resolution of the raster). The resolution is the
    same in both directions.}
  \item{value}{the value given to the pixels inside the polygon
    (converted to integer).}
  \item{backgrd}{idem for the pixels outside the polygon.}
}
\details{
  The mask is returned as a matrix which is filled rowwise (in agreement
  with the convention used in rasters) and can be input into functions
  in \pkg{terra} (e.g., \code{rast()}).

  \code{polygon2mask} does basically the same operation than
  \code{terra::rasterize} but is faster and can produce a vector for
  masking raster data.
}
\value{
  a matrix stored as integers; the dimensions of this matrix give the
  size of the raster.
}
\note{
  The code is still in development.
}
\author{Emmanuel Paradis}
\references{
  Nievergelt, J. and Preparata, F. P. (1982) Plane-sweep algorithms for
  intersecting geometric figures. \emph{Communications of the ACM},
  \bold{25}, 739--747. <doi:10.1145/358656.358681>.
}
\examples{
## from ?chullPolygon:
XY <- rbind(c(0, 0),
            c(1, 0),
            c(.25, .25),
            c(.5, .5),
            c(1.2, .8),
            c(1, .78),
            c(0, 1))

layout(matrix(1:9, 3, 3, TRUE))
k <- 2
for (i in 1:9) {
    msk <- polygon2mask(XY, k = k)
    d <- dim(msk)
    image(1:d[1], 1:d[2], msk)
    dm <- paste(d, collapse = "x")
    title(paste("k =", k, ", dim =", dm))
    k <- k * 2
}

layout(1)
}
\keyword{manip}