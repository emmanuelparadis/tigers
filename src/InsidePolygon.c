/* InsidePolygon.c       2023-10-06 */

/* Copyright 2023 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../DESCRIPTION for licensing issues. */

#include <R.h>
#include <Rinternals.h>

/* Returns 1 if the point with coordinates (Px,Py) is stricly inside
   the polygon, 0 if it is stricly outside.

   If (Px,Py) coincides with a vertex or is exactly on an edge of the
   polygon, the results are unpredictable and these must be tested
   separately.

   The polygon must be open and can be in either clockwise or
   counterclockwise order.

   n: number of vertices = number of edges (a.k.a. pathlength) */

#define INTERSECT_POINT_EDGE				\
if (Py > fmin(y1, y2)) {				\
    if (Py <= fmax(y1, y2)) {				\
	if (Px <= fmax(x1, x2)) {			\
	    if (y1 != y2) {				\
		A = (Py - y1)*(x2 - x1)/(y2 - y1) + x1;	\
		if (x1 == x2 || Px <= A) test = !test;	\
	    }						\
	}						\
    }							\
 }

/* maybe check the code in convexPolygonOverlap.c to see some bits can
   be used here too */

int InsidePolygon(double *x, double *y, int n, double Px, double Py)
{
    int i, j, test = 0;
    double A, x1, x2, y1, y2;

    x1 = x[0]; y1 = y[0];
    for (i = 0, j = 1; j < n; i++, j++) {
	x2 = x[j]; y2 = y[j];
	INTERSECT_POINT_EDGE
	x1 = x2; y1 = y2;
    }
    /* handle the last edge (n-1 -> 0) separately */
    x2 = x[0]; y2 = y[0];
    INTERSECT_POINT_EDGE

    return test;
}

SEXP InsidePolygon_Call(SEXP XY, SEXP P)
{
    double *x, *y, *Px, *Py;
    int i, n, np, *ptr;
    SEXP res;

    PROTECT(XY = coerceVector(XY, REALSXP));
    PROTECT(P = coerceVector(P, REALSXP));

    n = nrows(XY);
    x = REAL(XY);
    y = x + n;

    np = nrows(P);
    Px = REAL(P);
    Py = Px + np;

    PROTECT(res = allocVector(LGLSXP, np));
    ptr = INTEGER(res);

    for (i = 0; i < np; i++)
	ptr[i] = InsidePolygon(x, y, n, Px[i], Py[i]);

    UNPROTECT(3);
    return res;
}
