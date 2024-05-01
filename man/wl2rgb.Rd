\name{wl2col}
\alias{wl2col}
\alias{spectrum2col}
\alias{BlackBodySpectrum}
\title{Wavelengths to Colours}
\description{
  Conversion from wavelengths to colours.
}
\usage{
wl2col(x, gamma = 0.8, RGB = FALSE)
spectrum2col(spec, RGB = FALSE, no.warn = TRUE, color.system = 3)
BlackBodySpectrum(x, Temp = 300)
}
\arguments{
  \item{x}{a numeric vector with wavelengths in nanometers (nm).}
  \item{gamma}{parameter for correcting the transitions. If \code{gamma
      = 1}, then the transitions are linear (see examples).}
  \item{RGB}{a logical value. By default, colours (in HTML code) are
    returned. If \code{RGB = TRUE}, a matrix with three columns is
    returned.}
  \item{spec}{a numeric vector with 81 values giving the (relative)
    intensity of the different wave lengths between 380 nm and 780 nm
    (see examples).}
  \item{no.warn}{a logical value. If \code{TRUE} and some approximate
    calculations were performed in the C routine, a warning message is
    issued.}
  \item{color.system}{a single integer between 1 and 6 specifying the
    colour system (see details).}
  \item{Temp}{temperature in Kelvins (K) of the black body.}
}
\details{
  Computations are mainly performed by C and Fortran codes (see
  References).

  The argument \code{spec} gives the (relative) intensity of visible
  light between 380 nm and 780 nm in intervals with a bandwith of 5 nm
  (i.e., [380--385], [385--390], \dots, [775-780]). The returned value is
  the perceived colour of the given spectrum. It could happen that some
  calculations were approximate which is done silently unless
  \code{no.warn = FALSE}.

  The six colour systems are: (1) NTSC, (2) EBU (PAL/SECAM), (3) SMPTE,
  (4) HDTV, (5) CIE, and (6) CIE REC 709.

  \code{BlackBodySpectrum} calculates the emittance at specified
  wavelength(s) of a black body of temperature \code{Temp} using
  Planck's law.
}
\value{
  \code{wl2col} and \code{spectrum2col} return by default a vector of
  mode character with colours in HTML code. If \code{RGB = TRUE}, they
  return a matrix with the values (between 0 and 1) of red, green, and
  blue arranged in a three-column matrix. If the input \code{x} has
  names, these are used in the returned object (as names or rownames).

  \code{BlackBodySpectrum} returns a numeric vector.
}
\author{Emmanuel Paradis, John Walker, Dan Bruton}
\references{
  Bruton, D. (1996) Approximate RGB values for visible wavelengths.
  \url{http://www.physics.sfasu.edu/astro/color/spectra.html}

  Planck, M. (1901) Ueber das Gesetz der Energieverteilung im
  Normalspectrum. \emph{Annalen der Physik}, \bold{309},
  553--563. (English translation:
  \url{http://web.ihep.su/dbserv/compas/src/planck01/eng.pdf})

  Walker, J. (1996) Color rendering of spectra.
  \url{https://www.fourmilab.ch/documents/specrend/}
}
\examples{
wl <- 370:790
COLS <- c("red", "green", "blue")
if (interactive()) layout(matrix(1:3, 3))
matplot(wl, wl2col(wl, , TRUE), "l", col = COLS, lty = 1, lwd = 3)
title("gamma = 0.8 (default)")
matplot(wl, wl2col(wl, 1, TRUE), "l", col = COLS, lty = 1, lwd = 3)
title("gamma = 1")
matplot(wl, wl2col(wl, 1/3, TRUE), "l", col = COLS, lty = 1, lwd = 3)
title("gamma = 1/3")
layout(1)

spec <- numeric(81)
spec[2] <- 1
names(spec) <- seq(380, 780, 5)
sapply(1:6, function(i) spectrum2col(spec, TRUE, color.system = i))

WL <- 380:780
xlab <- "Wavelength (nm)"
ylab <- expression("Emittance (W."*m^{-2}*")")
plot(WL, BlackBodySpectrum(WL, 306), type = "l", xlab = xlab,
     ylab =ylab, log = "")
lines(WL, BlackBodySpectrum(WL, 303), lty = 2)
legend("topleft", legend = paste(c(306, 303) - 273, "degrees C"), lty = 1:2)


## vector of wavelengths:
wl <- seq(382.5, by = 5, length.out = 81)
spectrum2col(BlackBodySpectrum(wl, 310), TRUE)
spectrum2col(BlackBodySpectrum(wl, 3100), TRUE)
spectrum2col(BlackBodySpectrum(wl, 31000), TRUE)

wl <- 10:1e5
col <- wl2col(wl)
plot(wl, BlackBodySpectrum(wl, 6000), "n", log = "xy", ylim = c(1, 1e14),
     xaxs = "i", xlab = xlab, ylab = ylab)
s <- col != "#000000" # do not show the black lines
abline(v = wl[s], col = col[s])
lines(wl, BlackBodySpectrum(wl, 300))
lines(wl, BlackBodySpectrum(wl, 3000))
lines(wl, BlackBodySpectrum(wl, 6000))
text(c(3000, 240, 200), c(5e6, 6e10, 5e13), paste(c(300, 3000, 6000), "K"))
}
\keyword{manip}
