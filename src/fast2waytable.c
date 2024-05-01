/* fast2waytable.c    2024-01-05 */

/* Copyright 2021-2024 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rinternals.h>

SEXP fast2waytable_Call(SEXP X, SEXP Y, SEXP NCAT, SEXP TT)
{
    int i, K, n, *p, *x, *y;
    unsigned char *trans_tab;
    SEXP res;

    PROTECT(X = coerceVector(X, INTSXP));
    PROTECT(Y = coerceVector(Y, INTSXP));
    PROTECT(NCAT = coerceVector(NCAT, INTSXP));
    PROTECT(TT = coerceVector(TT, RAWSXP));

    x = INTEGER(X);
    y = INTEGER(Y);
    trans_tab = RAW(TT);

    n = LENGTH(X);
    K = INTEGER(NCAT)[0];

    PROTECT(res = allocMatrix(INTSXP, K, K));
    p = INTEGER(res);
    memset(p, 0, K * K * sizeof(int));

    for (i = 0; i < n; i++) (p[trans_tab[x[i]] + K * trans_tab[y[i]]])++;

    UNPROTECT(5);
    return res;
}

