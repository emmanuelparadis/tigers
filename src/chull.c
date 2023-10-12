/* chull.c       2023-10-06 */

/* Copyright 2023 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../DESCRIPTION for licensing issues. */

/* Reference: Graham and Yao 1983 J. Algo.

   This implementation requires the polygon to be open and in
   clockwise order (a crash will happen otherwise).

   It is possible to write a version for counterclockwise ordered by
   inverting the angles and directions. */

#include <R.h>
#include <Rinternals.h>
#include "tigers.h"

int findStartVertex(double *X, double *Y, int n, int *imax)
{
    int i, k = 0, m = 0;
    /* k: index of the lowest, leftest vertex
       m: index of the uppest, rightest vertex */

    for (i = 1; i < n; i++) {
	if (X[i] < X[k]) {
	    k = i;
	    continue;
	}
	if (X[i] > X[m]) {
	    m = i;
	    continue;
	}
	if (X[i] == X[k] && Y[i] < Y[k]) {
	    k = i;
	    continue;
	}
	if (X[i] == X[m] && Y[i] > Y[m]) {
	    m = i;
	    continue;
	}
    }

    *imax = m;
    return k;
}

/* returns 1, 0, or - 1, depending on whether c is to the right of,
   collinear with, or to the left of L[a, b] */
int angle_(double *X, double *Y, int a, int b, int c)
{
    double x;
    x = (X[c] - X[a]) * (Y[b] - Y[c]) - (X[c] - X[b]) * (Y[a] - Y[c]);
    if (x < 0) return -1;
    if (x > 0) return 1;
    return 0;
}

#define isRightOfLine(x, y, P) angle_(X, Y, x, y, P) > 0 ? 1 : 0

#define circ(x) x == -1 ? n - 1 : x

int LeftHull_(double *X, double *Y, int n, int *h, int start, int end)
{
    int i, t, A, B;

    /* Box I */
    h[0] = start;
    h[1] = end;
    i = circularIndex(end + 1, n);
    while (angle_(X, Y, h[0], h[1], i) <= 0)
	i = circularIndex(++i, n);
    h[2] = i;
    t = 2;

    for (;;) {
	if (i == start) break;

	i = circularIndex(++i, n); /* Box II */

	/* Box III */
	if (angle_(X, Y, h[t - 1], h[t], i) >= 0) {
	    /* below we need to consider that h[t] may be equal to 0 ('y' in Graham & Yao's paper) */
	    if (angle_(X, Y, circ(h[t] - 1), h[t], i) < 0) {
		A = h[t - 1];
		B = h[t];
	    } else {
		A = h[t];
		B = h[0];
	    }
	    while (isRightOfLine(A, B, i)) i = circularIndex(++i, n);
	}

	/* Box IV */
	while (angle_(X, Y, h[t - 1], h[t], i) <= 0) t--;
	t++;
	h[t] = i;
    }

    return t;
}

int convex_hull_GrahamYao(double *X, double *Y, int n, int *H, int isClockwise)
{
    int m, o, N, Nb, *buff;

    /* find the lowest (south), leftest (west) vertices, stored in o
       and m, respectively */
    o = findStartVertex(X, Y, n, &m);

    //    Rprintf("o = %d   m = %d\n", o, m);

    buff = (int*)R_alloc(n, sizeof(int));

    N = LeftHull_(X, Y, n, buff, o, m);
    N--;
    memcpy(H, buff + 1, N * sizeof(int));
    Nb = LeftHull_(X, Y, n, buff, m, o);
    Nb--;
    memcpy(H + N, buff + 1, Nb * sizeof(int));
    N += Nb;

    return N;
}

SEXP convex_hull_C(SEXP XY, SEXP isClockwise)
{
    int n, m, *H, *p;
    double *X, *Y;
    SEXP res;

    PROTECT(XY = coerceVector(XY, REALSXP));
    PROTECT(isClockwise = coerceVector(isClockwise, INTSXP));
    n = nrows(XY);
    X = REAL(XY);
    Y = X + n;

    H = (int*)R_alloc(n, sizeof(int));
    m = convex_hull_GrahamYao(X, Y, n, H, INTEGER(isClockwise)[0]);
    PROTECT(res = allocVector(INTSXP, m));
    p = INTEGER(res);
    memcpy(p, H, m * sizeof(int));

    UNPROTECT(3);
    return res;
}
