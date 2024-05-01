/* ECEF_lonlat.c    2024-01-28 */

/* Copyright 2024 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rinternals.h>

/* the constants (GRS80) */
#define a 6378137                   /* equatorial radius in meters */
#define f 0.003352810681183637418   /* flattening */
#define b 6356752.3                 /* polar radius = a * (1 - f) */
#define e2 0.0066943800229034       /* first eccentricity squared = f * (2 - f) */
#define ep2 0.0067394967754817      /* second eccentricity squared = 1 / pow(1 - f, 2) - 1 */

/* The notation is the same than in Blewitt (2024) */

/* additional constants for this function */
#define rad2deg 57.2957795130823229   /* 180/pi */
#define Const1 0.9966471893188164     /* 1 - f */
#define Const2 42697.6729161411276436 /* e2 * a */
#define Const3 12756274               /* 2 * a */
#define Const4 42841.3117236846446758 /* ep2 * b */

inline void ECEF2lonlat_(double x, double y, double z, double *res)
{
    double xy2, z2, p, r, P, R, C, S, T, lambda, phi, h;

    lambda = atan2(y, x);

    xy2 = pow(x, 2) + pow(y, 2);
    z2 = pow(z, 2);
    p = sqrt(xy2);
    r = sqrt(xy2 + z2);

    P = (p/Const1) * (1 - Const2/(r + f * pow(z/r, 2) * (Const3 - r)));
    R = sqrt(pow(P, 2) + z2);
    C = P / R;
    S = z / R;
    T = (z + Const4 * pow(S, 3)) / (p - Const2 * pow(C, 3));

    phi = atan(T);

    C = 1 / sqrt(1 + pow(T, 2));
    S = T * C;
    h = p*C + z*S - a*sqrt(1 - e2 * pow(S, 2));

    res[0] = lambda * rad2deg;
    res[1] = phi * rad2deg;
    res[2] = h;
}

SEXP ECEF2lonlat_Call(SEXP XYZ)
{
    int n, i;
    double *x, *y, *z, vec[3], *ptx, *pty, *ptz;
    SEXP res;
    PROTECT(XYZ = coerceVector(XYZ, REALSXP));
    n = nrows(XYZ);
    x = REAL(XYZ);
    y = x + n;
    z = y + n;

    PROTECT(res = allocMatrix(REALSXP, n, 3));
    ptx = REAL(res);
    pty = ptx + n;
    ptz = pty + n;

    for (i = 0; i < n; i++) {
	ECEF2lonlat_(x[i], y[i], z[i], vec);
	ptx[i] = vec[0];
	pty[i] = vec[1];
	ptz[i] = vec[2];
    }

    UNPROTECT(2);
    return res;
}
