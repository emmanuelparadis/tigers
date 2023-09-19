#include <R.h>
#include <Rinternals.h>

/* how to handle NA's?? if any... */

int RMA_(double *x, double *y, int n, double *out)
{
    int i;
    long double A, B;
    long double mx = 0, my = 0, vx = 0, vy = 0, cxy = 0;

    for (i = 0; i < n; i++) {
	mx += x[i];
	my += y[i];
    }
    mx /= n;
    my /= n;
    for (i = 0; i < n; i++) {
	A = x[i] - mx;
	B = y[i] - my;
	vx += A * A;
	vy += B * B;
	cxy += A * B;
    }
    vx /= n - 1;
    vy /= n - 1;
    cxy /= n - 1;

    if (!cxy) {
	out[0] = (double)my;
	out[1] = (double)mx;
	out[2] = 0;
	out[3] = R_PosInf;
    } else {
        B = 0.5 * (vy - vx ) / cxy;
	A = sqrtl(B * B + 1);
	out[2] = (double)B - A;
	out[3] = (double)B + A;
	out[0] = (double)my - out[2] * mx;
	out[1] = (double)my - out[3] * mx;
    }

    return 0;
}

SEXP RMA_Call(SEXP X, SEXP Y)
{
    int n, y_is_null;
    double *x, *y, *o;
    SEXP res;

    PROTECT(X = coerceVector(X, REALSXP));
    x = REAL(X);

    y_is_null = isNull(Y);

    if (y_is_null) {
	n = nrows(X);
	y = x + n;
    } else {
	PROTECT(Y = coerceVector(Y, REALSXP));
	n = LENGTH(X);
	if (LENGTH(Y) != n)
	    error("both vectors must have the same length\n");
	y = REAL(Y);
    }

    PROTECT(res = allocMatrix(REALSXP, 2, 2));
    o = REAL(res);

    RMA_(x, y, n, o);

    if (y_is_null) UNPROTECT(2); else UNPROTECT(3);
    return res;
}
