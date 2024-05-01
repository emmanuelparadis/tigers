\name{lonlat2UTM}
\alias{lonlat2UTM}
\alias{UTM2lonlat}
\title{Conversions of Coordinates}
\description{
  Functions to convert coordinates between angular (longitude, latitude)
  and UTM systems.
}
\usage{
lonlat2UTM(lon, lat = NULL, details = FALSE)
UTM2lonlat(x, y = NULL, zone = NULL, hemisphere = "N")
}
\arguments{
  \item{lon, lat}{coordinates given as two vectors or as a matrix (or
    a data frame).}
  \item{details}{a logical value indicating whether to return
    information on the UTM zones (see below).}
  \item{x, y}{either \code{x} is a data frame with four columns (as
    output from \code{lonlat2UTM}), or a matrix or a data frame giving
    the UTM coordinates in the standard way.}
  \item{zone, hemisphere}{an integer and a character string specifying the
    UTM square; can be single values (in which case they are recycled
    for all coordinates).}
}
\details{
  \code{lonlat2UTM} works for all UTM zones. If the coordinates cover
  several UTM zones and/or hemispheres, the option \code{details} is
  switched to \code{TRUE}.

  The formulas are in Karney (2011) originally from Kr\"uger (1912). The
  error is less than one micrometre. Karney (2011) presents algorithms
  that require arbitrary precision real numbers for errors of a few
  nanometres.
}
\value{
  a matrix or a data frame.
}
\author{Emmanuel Paradis}
\references{
  Karney, C. F. F. (2011) Transverse Mercator with an accuracy of a few
  nanometers. \emph{Journal of Geodesy}, \bold{85}, 475--485.

  Kr\"uger, J. H. L. (1912) Konforme Abbildung des Erdellipsoids in der
  Ebene. New Series 52. Royal Prussian Geodetic Institute, Potsdam.

}
\keyword{manip}
