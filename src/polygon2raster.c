#include <R.h>
#include <Rinternals.h>

/* The coordinates of the polygon are converted beforehand (in R for
   the moment) so that (0,0) is the center of the top-left pixel of
   the raster; the X-coordinates increase eastward, and the
   Y-coordinates increase southward. The coordinates of the center of
   the bottom-right pixel of the raster is (NC - 1, NR - 1). The
   coordinates of all pixel centers are integers.

   The algorithm moves southward, tracing *horizontal lines* at the
   center of the pixels and finds intersections with the edges of the
   polygon. These intersections are then sorted in increasing order of
   their Y-coordinates (i.e., moving southward), then the pixels of
   the raster between two successive intersections are given the value
   given in 'v'.

   north: most northern line (so that all Y of the polygon are south to this line)
   south: most southern line (so that all Y of the polygon are north to this line)

   The polygon path must be closed (easier to scan the edges).
*/

typedef struct intersection {
    double x;
    double y;
} intersection;

static int comp__(const void *a, const void *b)
{
    double ya, yb;
    ya = ((struct intersection *)a)->y;
    yb = ((struct intersection *)b)->y;
    return (ya > yb) - (ya < yb);
}

int compdbl(const void *a, const void *b)
{
    int x = *((double*)a);
    int y = *((double*)b);
    return (x > y) - (x < y);
}

SEXP singlePolygon2raster(SEXP XY, SEXP PARS, SEXP raster)
{
    int b, i, j, k, n, h, l, nxt, NC, north, south, *z, y, v, vert;
    double *X, *Y, beta, alpha, *buffer_of_X, x, x1, x2;
    intersection *intsct;
    SEXP res;

    PROTECT(XY = coerceVector(XY, REALSXP));
    PROTECT(PARS = coerceVector(PARS, INTSXP));
    PROTECT(raster = coerceVector(raster, INTSXP));

    n = nrows(XY);
    X = REAL(XY);
    Y = X + n;
    NC = INTEGER(PARS)[0];
    north = INTEGER(PARS)[1];
    south = INTEGER(PARS)[2];
    v = INTEGER(PARS)[3];
    z = INTEGER(raster);

    /* try to allocate enough memory, maybe a better guess could be better... */
    intsct = (intersection*)R_alloc(1E6, sizeof(intersection));

    /* 1) Build the table of edges and intersections with horizontal lines
       h: the "highest" vertex (most north)
       l: the "lowest" vertex (most south) */
    k = 0;
    for (i = 0, j = 1; j < n; i++, j++) {
	if (Y[i] == Y[j]) continue; /* if the edge is horizontal, nothing to do */
	if (Y[i] < Y[j]) {
	    h = i; l = j;
	} else {
	    l = i; h = j;
	}
	nxt = ceil(Y[h]); /* the next horizontal line south of the "highest" vertex */
	/* if the "lowest" vertex is not south of this line, nothing to do */
	if (nxt > Y[l]) continue;
	vert = X[h] == X[l]; /* is the edge vertical? */
	if (!vert) { /* avoid calculating infinite slope (beta) */
	    beta = (Y[h] - Y[l]) / (X[h] - X[l]);
	    alpha = Y[h] - beta * X[h];
	}
	/* store as many intersections as necessary: */
	while (nxt < Y[l]) {
	    intsct[k].x = vert ? X[h] : (nxt - alpha) / beta;
	    intsct[k].y = nxt;
	    k++;
	    nxt++;
	}
    }

    /*
       There should be an even number of intersections for each horizontal line.
    */

    /* 2) sort the intersections in increasing order of the Y
       coordinates -- seems much more efficient than the code in
       getPixelsInPolygon() */
    qsort(intsct, k, sizeof(struct intersection), comp__);

    /* 3) move southward from the north limit of the raster */
    buffer_of_X = (double*)R_alloc(1000, sizeof(double));
    y = north;
    j = 0;
    while (intsct[j].y != y) j++;
    while (y < south) {
	b = 0;
	while (intsct[j].y == y) {
	    /* 3a) the values of Y are sorted, so we just move along 'intsct' */
	    buffer_of_X[b] = intsct[j].x;
	    j++;
	    b++;
	}
	/* 3b) sort the values of X in increasing order (i.e., from W to E) */
	qsort(buffer_of_X, b, sizeof(double), compdbl);
	for (i = 0; i < b; i += 2) {
	    x1 = floor(buffer_of_X[i]); /* the 1st pixel to be converted... */
	    x2 = ceil(buffer_of_X[i + 1]); /* ... and the last one */
	    for (x = x1; x < x2; x++) z[(int)x + NC * y] = v;
	}
    	y++;
    }

    PROTECT(res = allocVector(INTSXP, 1));
    INTEGER(res)[0] = 0;
    UNPROTECT(4);
    return res;
}
