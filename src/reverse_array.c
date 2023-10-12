/* reverse_array.c       2023-10-06 */

/* Copyright 2023 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../DESCRIPTION for licensing issues. */

#include <R.h>
#include <Rinternals.h>

void rev_double(double *x, int n)
{
    double c;
    int i, j;

    for (i = 0, j = n - 1; i < n / 2; i++, j--) {
        c = x[i];
        x[i] = x[j];
        x[j] = c;
    }
}

void rev_2cols_double(double *x, int n)
{
    rev_double(x, n);
    rev_double(x + n, n);
}

void rev_copy_double(double *x, double *z, int n)
{
    int i, j;

    for (i = 0, j = n - 1; i < n; i++, j--) z[i] = x[j];
}

void rev_copy_2cols_double(double *x, double *z, int n)
{
    rev_copy_double(x, z, n);
    rev_copy_double(x + n, z + n, n);
}

SEXP rev_Call(SEXP x, SEXP copy)
{
    int n;
    SEXP res;
    double *xp;

    PROTECT(x = coerceVector(x, REALSXP));
    PROTECT(copy = coerceVector(copy, INTSXP));

    n = LENGTH(x);
    xp = REAL(x);

    if (INTEGER(copy)[0]) {
        double *z;
        PROTECT(res = allocVector(REALSXP, n));
        z = REAL(res);
	rev_copy_double(xp, z, n);
	/* alternative (~30% slower):
	   memcpy(z, xp, n * sizeof(double));
	   rev_double(z, n); */
    } else {
	PROTECT(res = allocVector(INTSXP, 1));
	rev_double(xp, n);
	INTEGER(res)[0] = 0;
    }

    UNPROTECT(3);
    return res;
}

SEXP rev_2cols_Call(SEXP x, SEXP copy)
{
    int n;
    SEXP res;
    double *xp;

    PROTECT(x = coerceVector(x, REALSXP));
    PROTECT(copy = coerceVector(copy, INTSXP));

    n = nrows(x);
    xp = REAL(x);

    if (INTEGER(copy)[0]) {
        double *z;
        PROTECT(res = allocMatrix(REALSXP, n, 2));
        z = REAL(res);
	rev_copy_2cols_double(xp, z, n);
    } else {
	PROTECT(res = allocVector(INTSXP, 1));
	rev_2cols_double(xp, n);
	INTEGER(res)[0] = 0;
    }

    UNPROTECT(3);
    return res;
}
