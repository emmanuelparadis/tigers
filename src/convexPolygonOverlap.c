#include <R.h>
#include <Rinternals.h>
#include "tigers.h"

typedef struct DATA {
    double x;
    int o;
} DATA;

static int comp__(const void *a, const void *b)
{
    double x, y;
    x = ((struct DATA *)a)->x;
    y = ((struct DATA *)b)->x;;
    return (x > y) - (x < y);
}

static int compd__(const void *a, const void *b)
{
    double x, y;
    x = ((struct DATA *)a)->x;
    y = ((struct DATA *)b)->x;;
    return (x < y) - (x > y);
}

void order_(double *x, int n, int *o, int rev)
{
    int i;
    struct DATA *X;

    X = (DATA*)R_alloc(n, sizeof(DATA));

    for (i = 0; i < n; i++) {
	X[i].x = x[i];
	X[i].o = i;
    }

    if (rev) {
	qsort(X, n, sizeof(struct DATA), comp__);
    } else {
	qsort(X, n, sizeof(struct DATA), compd__);
    }

    for (i = 0; i < n; i++) o[i] = X[i].o;
}

/* returns 1 if the bounding boxes of two segments overlap
   (returns 0 if the boxes are only contiguous) */
int overlappingBbox(double xA0, double yA0, double xA1, double yA1,
		    double xB0, double yB0, double xB1, double yB1)
{
    if (fmin(xA0, xA1) >= fmax(xB0, xB1)) return 0;
    if (fmin(xB0, xB1) >= fmax(xA0, xA1)) return 0;
    if (fmin(yA0, yA1) >= fmax(yB0, yB1)) return 0;
    if (fmin(yB0, yB1) >= fmax(yA0, yA1)) return 0;
    return 1;
}

/* returns a single letter describing the shape of a segment:
   'V': vertical
   'H': horizontal
   'P': a (single) point
   'O': oblique (or ordinary) */
char segment_shape(double x0, double y0, double x1, double y1)
{
    if (x0 == x1) {
	if (y0 == y1) return 'P';
	return 'V';
    }
    if (y0 == y1) return 'H';
    return 'O';
}

/* returns the slope of a segment giving its 4 coordinates and its shape */
double get_slope(double x0, double y0, double x1, double y1, char shape)
{
    if (shape == 'V') return INFINITY;
    if (shape == 'H') return 0;
    return (y1 - y0) / (x1 - x0);
}

/* returns 1 if two segments intersect (returns 0 if "simple contact");
   the coordinates of the intersection are returned in (*OX, *OY) */
int segmentIntersection(double xA0, double yA0, double xA1, double yA1,
			double xB0, double yB0, double xB1, double yB1,
			double *OX, double *OY)
{
    char S1, S2;
    double beta1, beta2, alpha1, alpha2, xi, yi;

    if (!overlappingBbox(xA0, yA0, xA1, yA1, xB0, yB0, xB1, yB1)) return 0;

    S1 = segment_shape(xA0, yA0, xA1, yA1);
    S2 = segment_shape(xB0, yB0, xB1, yB1);

    if (S1 == 'P' || S2 == 'P') return 0;
    if (S1 == 'V' && S2 == 'V') return 0;
    if (S1 == 'H' && S2 == 'H') return 0;

    /* get the slopes */
    beta1 = get_slope(xA0, yA0, xA1, yA1, S1);
    beta2 = get_slope(xB0, yB0, xB1, yB1, S2);
    if (beta1 == beta2) return 0;

    if (S1 == 'O') alpha1 = yA0 - beta1 * xA0;
    if (S2 == 'O') alpha2 = yB0 - beta2 * xB0;

    /* get the intersection point of the lines defined by both
       segments, there are 7 cases to consider */
    for (;;) {
	if (S2 == 'O' && S1 == 'O') { /* the most general case */
	    xi = (alpha1 - alpha2) / (beta2 - beta1);
	    yi = beta1 * xi + alpha1;
	    break;
	}
	if (S1 == 'V') { /* 2 cases */
	    xi = xA0;
	    yi = S2 == 'H' ? yB0 : beta2 * xi + alpha2;
	    break;
	}
	if (S2 == 'V') { /* 2 cases */
	    xi = xB0;
	    yi = S1 == 'H' ? yA0 : beta1 * xi + alpha1;
	    break;
	}
	if (S1 == 'H' && S2 == 'O') { /* 1 case */
	    yi = yA0;
	    xi = (yi - alpha2) / beta2;
	    break;
	}
	if (S2 == 'H' && S1 == 'O') { /* 1 case */
	    yi = yB0;
	    xi = (yi - alpha1) / beta1;
	    break;
	}
    }
    /* check that the intersection point is on within both segments */
    if (xi > xA0 && xi > xA1) return 0;
    if (xi < xA0 && xi < xA1) return 0;
    if (xi > xB0 && xi > xB1) return 0;
    if (xi < xB0 && xi < xB1) return 0;
    if (yi > yA0 && yi > yA1) return 0;
    if (yi < yA0 && yi < yA1) return 0;
    if (yi > yB0 && yi > yB1) return 0;
    if (yi < yB0 && yi < yB1) return 0;
    *OX = xi;
    *OY = yi;
    return 1;
}


/* compute half the signed area of the triangle */
/* #define test_area_triangle(x1, y1, x2, y2, x3, y3) (x1 - x3)*(y2 - y1) - (x1 - x2)*(y3 - y1) == 0 */

/* see chull.c for details */
int isLeftOfLine(double x1, double y1, double x2, double y2,
		 double x3, double y3)
{
    double x;

    x = (x3 - x1) * (y2 - y3) - (x3 - x2) * (y1 - y3);
    if (x < 0) return 1;
    return 0;
}

/* Find the (convex) polygon which is the overlap of two _convex_ polygons.

  Input: polygons A, B and their number of vertices (nA and nB)
  Output: polygon O
  Return value: number of vertices in O

  A and B must be in clockwise order (e.g., output from chull.c).

  It is possible to change the code to work with polygons in CCWO by
  changing isLeftOfLine().
*/

#define OUTPUT_VERTEX(x, y) buff[w] = x; ybuff[w] = y; w++

/* avoid duplicate vertices (note we don't do the same when scanning A) */
#define SCAN_VERTICES_B							\
    for (j = 0; j < nB; j++) {						\
        if ((B[j] == x0 && yB[j] == y0) || (B[j] == x1 && yB[j] == y1)) { \
	    dropB[j] = 1;						\
	    continue;							\
	}								\
	if (isLeftOfLine(x0, y0, x1, y1, B[j], yB[j])) {		\
	    dropB[j] = 1;						\
	    ii = -1;							\
	    while (ii < 2) {						\
		l = circularIndex(j + ii, nB);				\
		if (segmentIntersection(x0, y0, x1, y1, B[j], yB[j], B[l], yB[l], &xi, &yi)) { \
		    OUTPUT_VERTEX(xi, yi);				\
		}							\
		ii += 2;						\
	    }								\
	}								\
    }									\

/* no need to check for edge crossing, so we can skip vertices already marked */
#define SCAN_VERTICES_A					\
    for (j = 0; j < nA; j++) {				\
	if (dropA[j]) continue;				\
	if (isLeftOfLine(x0, y0, x1, y1, A[j], yA[j]))	\
	    dropA[j] = 1;				\
    }


int convexPolygonOverlap(double *A, double *B, int nA, int nB, double *O, int N)
{
    int i, j, k, l, ii, w = 0, *dropA, *dropB, *o;
    double x0, y0, x1, y1, xi, yi, *yA, *yB, mX = 0, mY = 0, *angle, *buff, *ybuff;

    yA = A + nA;
    yB = B + nB;

    buff = (double*)R_alloc(2 * N, sizeof(double));
    ybuff = buff + N;

    /* two indicators to identify the vertices that will be included
       or not in the output */
    dropA = (int*)R_alloc(nA, sizeof(int));
    dropB = (int*)R_alloc(nB, sizeof(int));
    memset(dropA, 0, nA * sizeof(int));
    memset(dropB, 0, nB * sizeof(int));

    /* The 1st loop scans the edges of A and marks the vertices of B
       on its left: these won't be included in the output (vertices on
       the edge are included). It also scans for edge crossings and
       output the intersection points. */
    for (i = 0, k = 1; k < nA; i++, k++) {
	x0 = A[i]; y0 = yA[i]; x1 = A[k]; y1 = yA[k];
	SCAN_VERTICES_B
    }
    /* handle the last edge of A (i == nA - 1) to complete the 1st loop */
    x0 = A[i]; y0 = yA[i]; x1 = A[0]; y1 = yA[0];
    SCAN_VERTICES_B

    /* The 2nd loop scans the edges of B and marks the vertices of A
       on its left: these won't be included in the output (vertices on
       the edge are included) */
    for (i = 0, k = 1; k < nB; i++, k++) {
	x0 = B[i]; y0 = yB[i]; x1 = B[k]; y1 = yB[k];
	SCAN_VERTICES_A
    }
    /* do the last edge of B (i == nB - 1) to complete the 2nd loop */
    x0 = B[i]; y0 = yB[i]; x1 = B[0]; y1 = yB[0];
    SCAN_VERTICES_A

    /* output the unmarked vertices of both polygons */
    for (i = 0; i < nA; i++) {
	if (dropA[i]) continue;
	OUTPUT_VERTEX(A[i], yA[i]);
    }
    for (i = 0; i < nB; i++) {
	if (dropB[i]) continue;
	OUTPUT_VERTEX(B[i], yB[i]);
    }

    /* sort the vertices of the output in clockwise order (to sort in
       counterclockwise order: replace "maximum" by "minimum") */

    /* 1. get the mean coordinates */
    for (i = 0; i < w; i++) {
	mX += buff[i];
	mY += ybuff[i];
    }
    mX /= w;
    mY /= w;

    /* 2. get the angle of each vertex with respect to the mean point */
    angle = (double*)R_alloc(w, sizeof(double));
    for (i = 0; i < w; i++)
	angle[i] = atan2(ybuff[i] - mY, buff[i] - mX);

    /* 3. sort in decreasing order of the angles */
    o = (int*)R_alloc(w, sizeof(int));
    order_(angle, w, o, 1);

    /* 4. copy the coordinates in O */
    for (i = 0; i < w; i++) {
	k = o[i];
	O[i] = buff[k];
	O[i + w] = ybuff[k];
    }

    return w;
}

SEXP convexPolygonOverlap_Call(SEXP A, SEXP B)
{
    int nA, nB, N, n;
    double *O, *pA, *pB;
    SEXP res;

    PROTECT(A = coerceVector(A, REALSXP));
    PROTECT(B = coerceVector(B, REALSXP));
    nA = nrows(A);
    nB = nrows(B);
    pA = REAL(A);
    pB = REAL(B);
    N = 2 * (nA + nB);
    O = (double*)R_alloc(2 * N, sizeof(double));
    n = convexPolygonOverlap(pA, pB, nA, nB, O, N);
    PROTECT(res = allocMatrix(REALSXP, n, 2));
    memcpy(REAL(res), O, 2 * n * sizeof(double));

    UNPROTECT(3);
    return res;
}
