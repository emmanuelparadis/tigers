c wavelengthtocol.f    2022-05-25
c
c Copyright 1996 Dan Bruton, ported to R by Emmanuel Paradis
c
c This file is part of the R-package `tigers'.
c See the file ../COPYING for licensing issues.
c
c
c      RGB VALUES FOR VISIBLE WAVELENGTHS   by Dan Bruton (astro@tamu.edu)
c
c      This program can be found at
c      http://www.physics.sfasu.edu/astro/color.html
c      and was last updated on February 20, 1996.
c
c      This program will create a ppm (portable pixmap) image of a spectrum.
c      The spectrum is generated using approximate RGB values for visible
c      wavelengths between 380 nm and 780 nm.
c      NetPBM's ppmtogif can be used to convert the ppm image
c      to a gif.  The red, green and blue values (RGB) are
c      assumed to vary linearly with wavelength (for GAMMA=1).
c      NetPBM Software: ftp://ftp.cs.ubc.ca/ftp/archive/netpbm/
c
      subroutine wltocol(WAVELEN, N, GAMMA, RES)
      DOUBLE PRECISION WAVELEN(N), GAMMA, RES(N,3), WL, SSS
      INTEGER N, I
      DO I=1,N
         WL = WAVELEN(I)
         IF ((WL.GE.380.).AND.(WL.LE.440.)) THEN
            R = -1.*(WL-440.)/(440.-380.)
            G = 0.
            B = 1.
         ENDIF
         IF ((WL.GE.440.).AND.(WL.LE.490.)) THEN
            R = 0.
            G = (WL-440.)/(490.-440.)
            B = 1.
         ENDIF
         IF ((WL.GE.490.).AND.(WL.LE.510.)) THEN
            R = 0.
            G = 1.
            B = -1.*(WL-510.)/(510.-490.)
         ENDIF
         IF ((WL.GE.510.).AND.(WL.LE.580.)) THEN
            R = (WL-510.)/(580.-510.)
            G = 1.
            B = 0.
         ENDIF
         IF ((WL.GE.580.).AND.(WL.LE.645.)) THEN
            R = 1.
            G = -1.*(WL-645.)/(645.-580.)
            B = 0.
         ENDIF
         IF ((WL.GE.645.).AND.(WL.LE.780.)) THEN
            R = 1.
            G = 0.
            B = 0.
         ENDIF
c
c      LET THE INTENSITY SSS FALL OFF NEAR THE VISION LIMITS
c
         IF (WL.GT.700.) THEN
            SSS=.3+.7* (780.-WL)/(780.-700.)
         ELSE IF (WL.LT.420.) THEN
            SSS=.3+.7*(WL-380.)/(420.-380.)
         ELSE
            SSS=1.
         ENDIF
c
c      GAMMA ADJUST
c
         RES(I,1)=(SSS*R)**GAMMA
         RES(I,2)=(SSS*G)**GAMMA
         RES(I,3)=(SSS*B)**GAMMA
      ENDDO
      RETURN
      END
