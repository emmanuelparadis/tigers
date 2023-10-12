/* sim_polygons_landClasses.c       2023-10-06 */

/* Copyright 2023 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../DESCRIPTION for licensing issues. */

/* #include "sim_polygons_landClasses.h" */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "tigers.h"

#define BonTWO 247766915924930.75 /* B/2 */
#define Cspheroid 0.08209454909747739 /* C */
#define R 6371008 /* mean radius of the Earth (m) */
#define RSQDTWO 20294871468032 /* = 0.5 * R^2 (m^2) */
#define area_triangle(x1, y1, x2, y2, x3, y3) fabs((x1 - x3)*(y2 - y1) - (x1 - x2)*(y3 - y1))/2

int find_position_minimum(double *x, int n)
{
    int i, res = 0;
    double min = x[0];

    for (i = 1; i < n; i++) {
	if (x[i] < min) {
	    res = i;
	    min = x[i];
	}
    }

    return res;
}

int find_position_maximum(double *x, int n)
{
    int i, res = 0;
    double max = x[0];

    for (i = 1; i < n; i++) {
	if (x[i] > max) {
	    res = i;
	    max = x[i];
	}
    }

    return res;
}

/* is P inside the triangle (A,B,C)? */
int isInsideTriangle(int A, int B, int C, int P, double *x, double *y)
{
    int k, O, a, b, side, above, above2;
    long double tmp, tmp2, denom, numer, beta, alpha, alpha2;
    for (k = 1; k < 4; k++) {
        if (k == 1) {O = A; a = B; b = C;}
        if (k == 2) {O = B; a = A; b = C;}
        if (k == 3) {O = C; a = A; b = B;}
        /* O is the point opposite to the base (a,b) */
        denom = x[a] - x[b];
        if (!denom) {
            tmp = x[P] - x[a];
            side = tmp > 0;
            tmp2 = x[O] - x[a];
            if (side != (tmp2 > 0)) return 0;
            if (fabsl(tmp2) < fabsl(tmp)) return 0;
        } else {
            numer = y[a] - y[b];
            if (!numer) {
                tmp = y[P] - y[a];
                side = tmp > 0;
                tmp2 = y[O] - y[a];
                if (side != (tmp2 > 0)) return 0;
                if (fabsl(tmp2) < fabsl(tmp)) return 0;
            } else {
                beta = numer / denom;
                alpha = y[b] - beta * x[b];
                /* is P on the same side than O wrt (a,b)? */
                above = y[P] - (beta * x[P] + alpha) > 0;
                if (above != (y[O] - (beta * x[O] + alpha) > 0)) return 0;
		alpha2 = y[O] - beta * x[O];
                above2 = y[P] - (beta * x[P] + alpha2) > 0;
                if (above == above2) return 0;
	    }
        }
    }
    return 1;
}

/* the path is closed, so the first and the last vertex are the same point */
long double areaPolygon2(int pathlength, int *path, double *x, double *y)
{
    int i, a, b;
    long double s = 0;
    a = path[0];
    for (i = 1; i < pathlength; i++) {
	/* a = path[i - 1]; */
	b = path[i];
	s += x[a] * y[b] - x[b] * y[a];
	a = b;
    }
    //return fabsl(s)/2;
    return s;
}

/* P must be < pathlength - 2 (see below) */
long double areaPolygon_drop1(int pathlength, int *path, double *x, double *y, int P)
{
    int a, b, p;
    long double res;
    a = P == 0 ? path[pathlength - 2] : path[P - 1];
    b = path[P + 1]; /* this works in all cases because the path is closed */
    p = path[P];
    res = x[p]*(y[b] - y[a]) + y[p]*(x[a] - x[b]) - (x[a]*y[b] - y[a]*x[b]);
    return res;
}

/* same than areaPolygon2() but with angular (lon,lat) coordinates */
/* x (= lon) and y (= lat) must be in radians */
long double areaPolygon2_angular(int pathlength, int *path, double *x, double *y)
{
    int i, a, b, c;
    long double s = 0;
    a = path[0];
    b = path[1];
    c = path[2];
    for (i = 1; i < pathlength; i++) {
	s += (x[c] - x[a]) * sin(y[b]);
	a = b;
	b = c;
	c = path[i + 2];
    }
    //return fabsl(s * RSQDTWO);
    return s;
}

/* same than areaPolygon2_drop1() but with angular (lon,lat) coordinates */
/* P must be >= 0 AND < pathlength */
long double areaPolygon_drop1_angular(int pathlength, int *path, double *x, double *y, int P)
{
    int a, b, c, d, p, N = pathlength - 1;
    long double res;
    a = path[(P - 2) % N];
    b = path[(P - 1) % N];
    c = path[(P + 1) % N];
    d = path[(P + 2) % N];
    p = path[P];
    /* there are three terms to delete, and two new terms */
    res  = (x[p] - x[a]) * sin(y[b]);
    res += (x[c] - x[b]) * sin(y[p]);
    res += (x[d] - x[p]) * sin(y[c]);
    res -= (x[c] - x[a]) * sin(y[b]);
    res -= (x[d] - x[b]) * sin(y[c]);
    return res;
}

/* P must be < pathlength - 2 */
int isEar(int pathlength, int *path, double *x, double *y, int P, long double S, long double *T)
{
    int A, B, i;
    long double D;
    D = areaPolygon_drop1(pathlength, path, x, y, P);
    //    Rprintf("S = %1.20Lf\tD = %1.20Lf\t%d\n", S, D, fabsl(S - D) > fabsl(S));
    if (fabsl(S - D) > fabsl(S)) return 0;
    A = P == 0 ? pathlength - 2 : P - 1;
    B = P + 1;
    for (i = 0; i < pathlength; i++) {
	if (i == A || i == B || i == P) continue;
	if (isInsideTriangle(path[A], path[B], path[P], path[i], x, y)) {
	    //	    Rprintf("i=%d  A=%d  B=%d  P=%d\n", i, A, B, P);
	    return 0;
	}
    }
    T[0] = D;
    //    Rprintf("OREILLE !!!\n");
    return 1;
}

/* the path must be closed */
int * triangulate_polygon(int pathlength, int *path, double *x, double *y)
{
    int i, j, m, n, N, *newpath, *res;
    long double Parea, D;
    n = pathlength;
    Parea = areaPolygon2(n, path, x, y);
    N = n - 3; /* number of triangles  */
    res = (int*)R_alloc(3 * N, sizeof(int));
    newpath = (int*)R_alloc(n, sizeof(int));
    memcpy(newpath, path, n * sizeof(int));
    m = 0;
    while (n > 4) {
	//	Rprintf("n = %d\tParea = %Lf\n", n, Parea);
	for (i = 1; i < n; i++) { /* start at the second point */
	    if (isEar(n, newpath, x, y, i, Parea, &D)) {
		//		Rprintf("*** fabsl(D) = %0.20Lf\t %d ***\n", fabsl(D), fabsl(D) > 0);
		res[m] = newpath[i];
		res[m + N] = newpath[i - 1];
		res[m + 2*N] = newpath[i + 1];
		/* update the path by removing the vertex */
		for (j = i; j < n - 1; j++) newpath[j] = newpath[j + 1];
		m++;
		n--;
		Parea -= D;
		break;
	    }
	}
	/* if (n == 6) { */
	/*     for (int ii = 0; ii < 6; ii++) */
	/* 	Rprintf("x=%f  y=%f\n", x[newpath[ii]], y[newpath[ii]]); */
	/*     return res; */
	/* } */
    }
    /* export the last triangle */
    res[m] = newpath[0];
    res[m + N] = newpath[1];
    res[m + 2*N] = newpath[2];
    return res;
}

double angle_direction_change(int v1, int v2, int v3, double *x, double *y)
{
    double alpha, beta, res;

    alpha = atan2(y[v1] - y[v2], x[v1] - x[v2]);
    beta = atan2(y[v3] - y[v2], x[v3] - x[v2]);
    res = M_PI - fabs(alpha);
    if (alpha < 0) res = -res;
    res += beta;

    if (res > M_PI) {
	res -= M_2PI;
	return res;
    }

    if (res < -M_PI) res += M_2PI;

    return res;
}

/* the path of vertices must be open and in counterclockwise order*/
int * triangulate_polygon_thin(int pathlength, int *path, double *x, double *y)
{
    int i, j, m, n, N, a, b, c, *newpath, *res;
    double *angle;
    n = pathlength;
    N = n - 2; /* number of triangles */
    res = (int*)R_alloc(3 * N, sizeof(int));
    angle = (double*)R_alloc(n, sizeof(double));
    newpath = (int*)R_alloc(n, sizeof(int));
    memcpy(newpath, path, n * sizeof(int));

    /* compute all triangle angles */
    a = path[n - 1];
    b = path[0];
    for (i = 0; i < n - 1; i++) {
	c = path[i + 1];
	angle[i] = angle_direction_change(a, b, c, x, y);
	a = b;
	b = c;
    }
    c = path[0];
    angle[n - 1] = angle_direction_change(a, b, c, x, y);

    m = 0;
    while (n > 3) {
	i = find_position_maximum(angle, n);
	/* if (angle[i] < 0) */
	//Rprintf("i = %d, angle = %f, n = %d\n", i, angle[i], n);
	//for (int ii = 0; ii < n; ii++) Rprintf("%d ", newpath[ii]); Rprintf("\n");
	res[m] = newpath[i]; /* should be 0 <= i < n */
	res[m + N] = newpath[circularIndex(i - 1, n)]; /* in case i = 0 */
	res[m + 2*N] = newpath[circularIndex(i + 1, n)]; /* in case i = n - 1 */
	m++;
	/* update the path by removing vertex i */
	for (j = i; j < n - 1; j++) {
	    newpath[j] = newpath[j + 1];
	    angle[j] = angle[j + 1];
	}
	n--;

	/* update the two angles that are affected by the clipping */
	j = circularIndex(i - 1, n);
	//Rprintf("i = %d\tj = %d\tn = %d\n", i, j, n);
	a = newpath[circularIndex(i - 2, n)];
	b = newpath[j];
	c = newpath[i]; /* we should have 0 <= i < n */
	angle[j] = angle_direction_change(a, b, c, x, y);
	a = b;
	b = c;
	c = newpath[circularIndex(i + 1, n)];
	angle[i] = angle_direction_change(a, b, c, x, y);
    }

    /* export the last triangle */
    res[m] = newpath[0];
    res[m + N] = newpath[1];
    res[m + 2*N] = newpath[2];

    return res;
}

double calculate_fatness(int v1, int v2, int v3, double *x, double *y)
{
    double a, b, c, tmp;

    a = x[v1] - x[v2];
    a *= a;
    tmp = y[v1] - y[v2];
    a += tmp * tmp;
    a = sqrt(a);

    b = x[v1] - x[v3];
    b *= b;
    tmp = y[v1] - y[v3];
    b += tmp * tmp;
    b = sqrt(b);

    c = x[v2] - x[v3];
    c *= c;
    tmp = y[v2] - y[v3];
    c += tmp * tmp;
    c = sqrt(c);

    tmp = (a+b+c)/3;
    a -= tmp;
    b -= tmp;
    c -= tmp;
    return (a*a + b*b + c*c)/3;

    //    return a * (a/2 - b) + b * (b/2 - c) + c * (c/2 - a);
}

/* path of vertices must be open and in clockwise order */
int * triangulate_polygon_fat(int pathlength, int *path, double *x, double *y)
{
    int i, j, m, n, N, a, b, c, *newpath, *res;
    double *fatness;
    long double Parea, D;
    n = pathlength;
    Parea = areaPolygon2(n, path, x, y);
    N = n - 2; /* number of triangles  */
    res = (int*)R_alloc(3 * N, sizeof(int));
    fatness = (double*)R_alloc(n, sizeof(double));
    newpath = (int*)R_alloc(n, sizeof(int));
    memcpy(newpath, path, n * sizeof(int));

    /* compute all triangle fatness */
    a = path[n - 1];
    b = path[0];
    for (i = 0; i < n - 1; i++) {
	c = path[i + 1];
	fatness[i] = isEar(n, path, x, y, i, Parea, &D) ? calculate_fatness(a, b, c, x, y) : R_PosInf;
	a = b;
	b = c;
    }
    c = path[0];
    fatness[n - 1] = isEar(n, path, x, y, n - 1, Parea, &D) ? calculate_fatness(a, b, c, x, y) : R_PosInf;

    m = 0;
    while (n > 3) {
	i = find_position_minimum(fatness, n);
	//	Rprintf("i = %d   fatness[i] = %f\n", i, fatness[i]);
	res[m] = newpath[i];
	res[m + N] = newpath[circularIndex(i - 1, n)];
	res[m + 2*N] = newpath[circularIndex(i + 1, n)];
	/* update the path by removing the vertex i */
	for (j = i; j < n - 1; j++) {
	    newpath[j] = newpath[j + 1];
	}

	m++;
	n--;

	Parea = areaPolygon2(n, path, x, y);

	a = path[n - 1];
	b = path[0];
	for (i = 0; i < n - 1; i++) {
	    c = path[i + 1];
	    fatness[i] = isEar(n, newpath, x, y, i, Parea, &D) ? calculate_fatness(a, b, c, x, y) : R_PosInf;
	    a = b;
	    b = c;
	}
	c = path[0];
	fatness[n - 1] = isEar(n, newpath, x, y, n - 1, Parea, &D) ? calculate_fatness(a, b, c, x, y) : R_PosInf;
    }

    /* export the last triangle */
    res[m] = newpath[0];
    res[m + N] = newpath[1];
    res[m + 2*N] = newpath[2];

    return res;
}

int * triangulate_polygon_det(int pathlength, int *path, double *x, double *y)
{
    int i, j, m, n, N, a, b, c, *newpath, *res;
    double *det;
    n = pathlength;
    N = n - 2; /* number of triangles  */
    res = (int*)R_alloc(3 * N, sizeof(int));
    det = (double*)R_alloc(n, sizeof(double));
    newpath = (int*)R_alloc(n, sizeof(int));
    memcpy(newpath, path, n * sizeof(int));

    /* compute all triangle determinants (= signed areas) */
    a = path[n - 1];
    b = path[0];
    for (i = 0; i < n - 1; i++) {
	c = path[i + 1];
	det[i] = (x[b] - x[a])*(y[c] - y[a]) - (x[c] - x[a])*(y[b] - y[a]);
	a = b;
	b = c;
    }
    c = path[0];
    det[n - 1] = (x[b] - x[a])*(y[c] - y[a]) - (x[c] - x[a])*(y[b] - y[a]);

    m = 0;
    while (n > 3) {
	i = find_position_maximum(det, n);
	res[m] = newpath[i];
	res[m + N] = newpath[circularIndex(i - 1, n)];
	res[m + 2*N] = newpath[circularIndex(i + 1, n)];
	m++;
	/* update the path by removing the vertex i */
	for (j = i; j < n - 1; j++) {
	    newpath[j] = newpath[j + 1];
	    det[j] = det[j + 1];
	}
	n--;

	/* update the two angles that are affected */
	j = circularIndex(i - 1, n);
	a = newpath[circularIndex(i - 2, n)];
	b = newpath[j];
	c = newpath[i]; /* we should have 0 <= i < n */
	det[j] = (x[b] - x[a])*(y[c] - y[a]) - (x[c] - x[a])*(y[b] - y[a]);
	a = b;
	b = c;
	c = newpath[circularIndex(i + 1, n)];
	det[i] = (x[b] - x[a])*(y[c] - y[a]) - (x[c] - x[a])*(y[b] - y[a]);
    }

    /* export the last triangle */
    res[m] = newpath[0];
    res[m + N] = newpath[1];
    res[m + 2*N] = newpath[2];

    return res;
}

SEXP triangulate_Call(SEXP XY, SEXP METHOD)
{
    double *x, *y;
    int i, pathlength, *path, *res, method, S;
    SEXP v;

    /* S=3 if the path is closed, S=2 if open */

    PROTECT(XY = coerceVector(XY, REALSXP));
    PROTECT(METHOD = coerceVector(METHOD, INTSXP));
    pathlength = nrows(XY);
    path = (int*)R_alloc(pathlength, sizeof(int));
    for (i = 0; i < pathlength; i++) path[i] = i;
    x = REAL(XY);
    y = x + pathlength;
    method = INTEGER(METHOD)[0];
    switch(method) {
    case 1: {
	res = triangulate_polygon(pathlength, path, x, y);
	S = 3;
	break;
    }
    case 2: {
	res = triangulate_polygon_thin(pathlength, path, x, y);
	S = 2;
	break;
    }
    case 3: {
	res = triangulate_polygon_fat(pathlength, path, x, y);
	S = 2;
	break;
    }
    case 4: {
	res = triangulate_polygon_det(pathlength, path, x, y);
	S = 2;
	break;
    }
    }

    PROTECT(v = allocMatrix(INTSXP, pathlength - S, 3));
    memcpy(INTEGER(v), res, 3 * (pathlength - S) * sizeof(int));
    UNPROTECT(3);
    return v;
}

SEXP haveOverlapTwoPolygons(SEXP P, SEXP Q)
{
    int *tri, *path, i, j, n, m, N, done = 0;
    double *x, *y, *XX, *YY, *q;
    SEXP res;

    PROTECT(P = coerceVector(P, REALSXP));
    PROTECT(Q = coerceVector(Q, REALSXP));

    n = nrows(P);
    m = nrows(Q);

    //Rprintf("n = %d\tm = %d\n", n, m);

    x = (double*)R_alloc(n, sizeof(double));
    y = (double*)R_alloc(n, sizeof(double));
    XX = (double*)R_alloc(4, sizeof(double));
    YY = (double*)R_alloc(4, sizeof(double));

    N = n * sizeof(double); /* temporary use */
    memcpy(x, REAL(P), N);
    memcpy(y, REAL(P) + n, N);
    //for (i =0; i<n;i++) Rprintf("%f\t%f\n", x[i], y[i]);
    q = REAL(Q);

    /* check the paths are not "closed" */
    /* if (x[0] == x[n - 1] && y[0] == y[n - 1]) n--; */
    /* if (q[0] == q[m - 1] && q[m] == q[2*m - 1]) m--; */

    /* triangulate the first polygon */
    if (n == 3) {
	N = 1; /* number of triangles  */
	tri = (int*)R_alloc(3, sizeof(int));
	tri[0] = 0; tri[1] = 1; tri[2] = 2;
    } else {
	//Rprintf("start triangulation...");
	N = n - 3; /* number of triangles  */
	path = (int*)R_alloc(n, sizeof(int));
	for (i = 0; i < n; i++) path[i] = i;
	//Rprintf(" calling triangulate_polygon()...");
	tri = triangulate_polygon(n, path, x, y);
	//Rprintf(" done.\n");
    }

    for (i = 0; i < N; i++) {
	//	Rprintf("N = %d\ti = %d\n", N, i);
	XX[0] = x[tri[i]];
	YY[0] = y[tri[i]];
	XX[1] = x[tri[i + N]];
	YY[1] = y[tri[i + N]];
	XX[2] = x[tri[i + 2*N]];
	YY[2] = y[tri[i + 2*N]];
	for (j = 0; j < m; j++) {
	    XX[3] = q[j]; YY[3] = q[j + m];
	    if (isInsideTriangle(0, 1, 2, 3, XX, YY)) {
		done = 1;
		break;
	    }
	}
	if (done) break;
    }

    PROTECT(res = allocVector(LGLSXP, 1));
    INTEGER(res)[0] = done;
    UNPROTECT(3);
    return res;
}

int check_identical_vertices(double *x, int n, int *red, int check)
{
    int i = 0, j, N = 0;

    for (;;) {
	j = i + 1;
	if (j >= n) break;
	while (x[i] == x[j] && x[i + n] == x[j + n]) {
	    red[j] = 1;
	    N++;
	    j++;
	    if (j >= n) break;
	}
	i = j;
    }

    if (!check) return N;

    if (N) {
	Rprintf("Found %d redundant ", N);
	N == 1 ? Rprintf("vertex") : Rprintf("vertices");
	Rprintf("\n(identical vertices are on the same line; indices are 0...n-1):\n");
	int flag = 0;
	for (i = 1; i < n - 1; i++) {
	    if (red[i]) {
		if (!flag) {
		    Rprintf("\n%d", i - 1);
		    flag = 1;
		}
		Rprintf(", %d", i);
	    } else {
		flag = 0;
	    }
	}
	Rprintf("\n");
    } else {
	Rprintf("No identical vertices.\n");
    }

    return N;
}

int check_close_vertices(double *x, int n, double tol, int *red, int check)
{
    int i = 0, j, N = 0;

    for (;;) {
	j = i + 1;
	if (j >= n) break;
	while (sqrt(pow(x[i] - x[j], 2) + pow(x[i + n] - x[j + n], 2)) <= tol) {
	    red[j] = 1;
	    N++;
	    j++;
	    if (j >= n) break;
	}
	i = j;
    }

    if (!check) return N;

    if (N) {
	Rprintf("Found %d close ", N);
	N == 1 ? Rprintf("vertex") : Rprintf("vertices");
	Rprintf(" (tolerance = %e)", tol);
	Rprintf("\n(close vertices are on the same line; indices are 0...n-1):\n");
	int flag = 0;
	for (i = 1; i < n - 1; i++) {
	    if (red[i]) {
		if (!flag) {
		    Rprintf("\n%d", i - 1);
		    flag = 1;
		}
		Rprintf(", %d", i);
	    } else {
		flag = 0;
	    }
	}
	Rprintf("\n");
    } else {
	Rprintf("No close vertices.\n");
    }

    return N;
}

void remove_vertices(double *x, int n, double *y, int *red)
{
    int i, j = 0;
    for (i = 0; i < n; i++) {
	if (red[i]) continue;
	y[j] = x[i];
	j++;
    }
    for (i = n; i < 2*n; i++) {
	if (red[i - n]) continue;
	y[j] = x[i];
	j++;
    }
}

SEXP redundant_vertices(SEXP POLYGON, SEXP TOL, SEXP CHECK_ONLY)
{
    double *x, *ptr, tol, N;
    int n, /* method, */ check, *red;
    SEXP res;

    PROTECT(POLYGON = coerceVector(POLYGON, REALSXP));
    /* PROTECT(METHOD = coerceVector(METHOD, INTSXP)); */
    PROTECT(TOL = coerceVector(TOL, REALSXP));
    PROTECT(CHECK_ONLY = coerceVector(CHECK_ONLY, LGLSXP));

    n = nrows(POLYGON);

    x = REAL(POLYGON);
    /* method = INTEGER(METHOD); */
    tol = REAL(TOL)[0];
    check = INTEGER(CHECK_ONLY)[0];

    red = (int*)R_alloc(n, sizeof(int));
    memset(red, 0, n * sizeof(int));

    N = tol == 0 ? check_identical_vertices(x, n, red, check) : check_close_vertices(x, n, tol, red, check);

    if (check) {
	PROTECT(res = allocVector(INTSXP, 1));
	INTEGER(res)[0] = 1;
    } else {
	if (N) {
	    PROTECT(res = allocMatrix(REALSXP, n - N, 2));
	    ptr = REAL(res);
	    remove_vertices(x, n, ptr, red);
	} else {
	    PROTECT(res = allocMatrix(REALSXP, n, 2));
	    ptr = REAL(res);
	    memcpy(ptr, x, 2 * n * sizeof(double));
	}
    }

    UNPROTECT(4);
    return res;
}

/* the path is open */
double areaPolygon2open(double *x, double *y, int n)
{
    int a, b;
    double s = 0;
    for (a = 0, b = 1; b < n; a++, b++)
	s += x[a] * y[b] - x[b] * y[a];
    b--; /* b = n - 1 */
    s += x[b] * y[0] - x[0] * y[b];
    return fabs(s)/2;
}

/* same than areaPolygon2open() but with angular (lon,lat) coordinates */
/* x (= lon) and y (= lat) must be in radians */
double areaPolygon2open_angular(double *x, double *y, int n)
{
    int a = 0, b = 1, c = 2;
    double s = 0;
    while (c < n) {
	s += (x[c] - x[a]) * sin(y[b]);
	a++; b++; c++;
    }
    /* handle the last two triplets */
    s += (x[0] - x[n - 2]) * sin(y[n - 1]);
    s += (x[1] - x[n - 1]) * sin(y[0]);
    return fabs(s * RSQDTWO);
}

double areaPolygon2open_spheroid(double *x, double *y, int n)
{
    int a = 0, b = 1, c = 2;
    double s = 0;
    while (c < n) {
	s += (x[c] - x[a]) * atan(Cspheroid * sin(y[b]));
	a++; b++; c++;
    }
    /* handle the last two triplets */
    s += (x[0] - x[n - 2]) * atan(Cspheroid * sin(y[n - 1]));
    s += (x[1] - x[n - 1]) * atan(Cspheroid * sin(y[0]));
    return fabs(s * BonTWO);
}

SEXP area_Call(SEXP XY, SEXP METHOD)
{
    SEXP res;
    double *x, *y;
    int n, *method;
    double (*FUN)(double *x, double *y, int n);

    PROTECT(XY = coerceVector(XY, REALSXP));
    PROTECT(METHOD = coerceVector(METHOD, INTSXP));
    method = INTEGER(METHOD);
    n = nrows(XY);
    x = REAL(XY);
    y = x + n;

    if (method[0]) {
	FUN = method[1] ? &areaPolygon2open_spheroid : &areaPolygon2open_angular;
    } else {
	FUN = areaPolygon2open;
    }

    PROTECT(res = allocVector(REALSXP, 1));
    REAL(res)[0] = FUN(x, y, n);

    UNPROTECT(3);
    return res;
}
