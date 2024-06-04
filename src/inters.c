/* inters.c       2024-05-31 */

/* Copyright 2024 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../DESCRIPTION for licensing issues. */

/*
  some code below is similar to some in tigers/src/convexPolygonOverlap.c
  but with a different parameterization
*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* returns 1 if the bounding boxes of two segments overlap
   (returns 0 if the boxes are only contiguous) */
int overlappingBboxBIS(double *A, double *B)
{
    if (fmin(A[0], A[2]) >= fmax(B[0], B[2])) return 0;
    if (fmin(B[0], B[2]) >= fmax(A[0], A[2])) return 0;
    if (fmin(A[1], A[3]) >= fmax(B[1], B[3])) return 0;
    if (fmin(B[1], B[3]) >= fmax(A[1], A[3])) return 0;
    return 1;
}

/* returns a single letter describing the shape of the segment:
   'V': vertical
   'H': horizontal
   'P': a (single) point
   'O': oblique (or ordinary) */
char segment_shapeBIS(double *x)
{
    if (x[0] == x[2]) {
	if (x[1] == x[3]) return 'P';
	return 'V';
    }
    if (x[1] == x[3]) return 'H';
    return 'O';
}

/* returns the slope of a segment giving its 4 coordinates and its shape */
double get_slopeBIS(double *x, char shape)
{
    if (shape == 'V') return INFINITY;
    if (shape == 'H') return 0;
    return (x[3] - x[1]) / (x[2] - x[0]);
}

/* returns 1 if two segments intersect (returns 0 if "simple contact");
   the coordinates of the intersection are returned in (*O) */
int segmentIntersectionBIS(double *A, double *B, double *O)
{
    char S1, S2;
    double beta1, beta2, alpha1, alpha2, xi, yi;

    if (!overlappingBboxBIS(A, B)) return 0;

    S1 = segment_shapeBIS(A);
    S2 = segment_shapeBIS(B);

    if (S1 == 'P' || S2 == 'P') return 0;
    if (S1 == 'V' && S2 == 'V') return 0;
    if (S1 == 'H' && S2 == 'H') return 0;

    /* get the slopes */
    beta1 = get_slopeBIS(A, S1);
    beta2 = get_slopeBIS(B, S2);
    if (beta1 == beta2) return 0;

    if (S1 == 'O') alpha1 = A[1] - beta1 * A[0];
    if (S2 == 'O') alpha2 = B[1] - beta2 * B[0];

    /* get the intersection point of the lines defined by both
       segments, there are 7 cases to consider */
    for (;;) {
	if (S2 == 'O' && S1 == 'O') { /* the most general case */
	    xi = (alpha1 - alpha2) / (beta2 - beta1);
	    yi = beta1 * xi + alpha1;
	    break;
	}
	if (S1 == 'V') { /* 2 cases */
	    xi = A[0];
	    yi = S2 == 'H' ? B[1] : beta2 * xi + alpha2;
	    break;
	}
	if (S2 == 'V') { /* 2 cases */
	    xi = B[0];
	    yi = S1 == 'H' ? A[1] : beta1 * xi + alpha1;
	    break;
	}
	if (S1 == 'H' && S2 == 'O') { /* 1 case */
	    yi = A[1];
	    xi = (yi - alpha2) / beta2;
	    break;
	}
	if (S2 == 'H' && S1 == 'O') { /* 1 case */
	    yi = B[1];
	    xi = (yi - alpha1) / beta1;
	    break;
	}
    }
    /* check that the intersection point is within both segments */
    if (xi > A[0] && xi > A[2]) return 0;
    if (xi < A[0] && xi < A[2]) return 0;
    if (xi > B[0] && xi > B[2]) return 0;
    if (xi < B[0] && xi < B[2]) return 0;
    if (yi > A[1] && yi > A[3]) return 0;
    if (yi < A[1] && yi < A[3]) return 0;
    if (yi > B[1] && yi > B[3]) return 0;
    if (yi < B[1] && yi < B[3]) return 0;
    O[0] = xi;
    O[1] = yi;
    return 1;
}

SEXP inters_seg_seg(SEXP SEGI, SEXP SEGII)
{
    SEXP res, INTERS, TEST;

    PROTECT(SEGI = coerceVector(SEGI, REALSXP));
    PROTECT(SEGII = coerceVector(SEGII, REALSXP));

    PROTECT(res = allocVector(VECSXP, 2));
    PROTECT(INTERS = allocVector(REALSXP, 2));
    PROTECT(TEST = allocVector(INTSXP, 1));

    INTEGER(TEST)[0] = segmentIntersectionBIS(REAL(SEGI), REAL(SEGII), REAL(INTERS));

    SET_VECTOR_ELT(res, 0, INTERS);
    SET_VECTOR_ELT(res, 1, TEST);

    UNPROTECT(5);
    return res;
}

int between(double x, double y, double *seg, char S)
{
    if (S == 'V') {
	if (y > seg[1] && y < seg[3]) return 1;
	if (y < seg[1] && y > seg[3]) return 1;
	return 0;
    }
    if (S == 'H') {
	if (x > seg[0] && x < seg[2]) return 1;
	if (x < seg[0] && x > seg[2]) return 1;
	return 0;
    }
    /* we know (x,y) is on the line, so just need to test either x or y */
    if (x > seg[0] && x < seg[2]) return 1;
    if (x < seg[0] && x > seg[2]) return 1;
    return 0;
}

/* A and B are two angles defining an arc with B->A in CWO.
   The test returns 1 if x is within the arc, 0 if not.
   The 3 input angles are outputs from atan2();
   they are first rescaled from [-pi,pi] to [0,2pi]. */
int between_angles(double x, double B, double A)
{
    if (x < 0) x += M_2PI;
    if (A < 0) A += M_2PI;
    if (B < 0) B += M_2PI;
    if (B > A) {
	if (x < B && x > A) return 1;
    } else {
	/* if B < A, then the arc overlaps with 0 */
	if (x < B) return 1;
	if (x > A) return 1;
    }
    return 0;
}

/* Returns 1|2 if a segment and an arc intersect (returns 0 if "simple
   contact"); the coordinates of the intersection(s) are returned in
   *O. The returned value is the number of intersections.
   The arc is parameterized in CWO. */
int segmentArcIntersection(double *seg, double *arc, double *O)
{
    int n = 0, test;
    char S;
    double beta, alpha, xi1, yi1, xi2, yi2, tmp, a, b, c, Delta,
	theta0, theta1, theta, circle_bbox[4], x_ctr, y_ctr, r;

    /* local copies of these 3 which are used several times */
    x_ctr = arc[0]; y_ctr = arc[1]; r = arc[2];

    /* 1) test if the bounding boxes overlap */
    circle_bbox[0] = x_ctr - r;
    circle_bbox[1] = y_ctr - r;
    circle_bbox[2] = x_ctr + r;
    circle_bbox[3] = y_ctr + r;
    if (!overlappingBboxBIS(seg, &circle_bbox[0])) return 0;

    S = segment_shapeBIS(seg);
    if (S == 'P') return 0;

    //    Rprintf("S = %c\n", S);

    for (;;) {
	/* 2a) the general case */
	if (S == 'O') {
	    /* get the slope and the intercept */
	    beta = get_slopeBIS(seg, S);
	    alpha = seg[1] - beta * seg[0];

	    /* the 3 coefficients of the polymial */
	    a = 1 + pow(beta, 2);
	    b = -2 * (x_ctr + beta * (y_ctr - alpha));
	    c = x_ctr*x_ctr + pow(y_ctr - alpha, 2) - r*r;

	    Delta = b*b - 4*a*c;
	    if (Delta <= 0) return 0;
	    tmp = sqrt(Delta);
	    xi1 = xi2 = -b;
	    xi1 -= tmp;
	    xi2 += tmp;
	    tmp = 2 * a;
	    xi1 /= tmp;
	    xi2 /= tmp;
	    yi1 = beta * xi1 + alpha;
	    yi2 = beta * xi2 + alpha;
	    break;
	}

	/* 2b) the special cases */
	if (S == 'V') {
	    if (fabs(seg[1] - y_ctr) >= r) return 0;
	    xi1 = xi2 = seg[0];
	    yi1 = yi2 = y_ctr;
	    tmp = sqrt(pow(r, 2) - pow(seg[0] - x_ctr, 2));
	    yi1 -= tmp;
	    yi2 += tmp;
	    break;
	}
	// if (S == 'H') {
	if (fabs(seg[0] - x_ctr) >= r) return 0;
	yi1 = yi2 = seg[1];
	xi1 = xi2 = x_ctr;
	tmp = sqrt(r*r - pow(seg[1] - y_ctr, 2));
	xi1 -= tmp;
	xi2 += tmp;
	break;
	//}
    }

    theta0 = atan2(arc[4] - y_ctr, arc[3] - x_ctr);
    theta1 = atan2(arc[6] - y_ctr, arc[5] - x_ctr);

    /* check that the 1st intersection point is on the segment... */
    if (between(xi1, yi1, seg, S)) {
	/* ... and on the arc */
	theta = atan2(yi1 - y_ctr, xi1 - x_ctr);
	test = between_angles(theta, theta0, theta1);
	if (test) {
	    O[0] = xi1;
	    O[1] = yi1;
	    n++;
	}
    }

    /* id. for the 2nd point */
    if (between(xi2, yi2, seg, S)) {
	theta = atan2(yi2 - y_ctr, xi2 - x_ctr);
	test = between_angles(theta, theta0, theta1);
	if (test) {
	    if (n == 1) {
		O[2] = xi2;
		O[3] = yi2;
	    } else  {
		O[0] = xi2;
		O[1] = yi2;
	    }
	    n++;
	}
    }

    return n;
}

SEXP inters_seg_arc(SEXP SEG, SEXP ARC)
{
    SEXP res, INTERS, TEST;

    PROTECT(SEG = coerceVector(SEG, REALSXP));
    PROTECT(ARC = coerceVector(ARC, REALSXP));

    PROTECT(res = allocVector(VECSXP, 2));
    PROTECT(INTERS = allocVector(REALSXP, 4));
    PROTECT(TEST = allocVector(INTSXP, 1));

    INTEGER(TEST)[0] = segmentArcIntersection(REAL(SEG), REAL(ARC), REAL(INTERS));

    SET_VECTOR_ELT(res, 0, INTERS);
    SET_VECTOR_ELT(res, 1, TEST);

    UNPROTECT(5);
    return res;
}
