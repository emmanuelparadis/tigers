/* UTM_lonlat.c    2024-04-30 */

/* Copyright 2024 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rinternals.h>

/* the constants (GRS80) */
#define a 6378137                  /* equatorial radius in meters */
#define f 0.003352810681183637418  /* flattening */
#define b a * (1 - f)              /* polar radius = a * (1 - f) */
#define e2 f * (2 - f)             /* first eccentricity squared = f * (2 - f) */
#define n f / (2 - f)              /* third flattening */
#define A (a / (1 + n)) * (1 + pow(n, 2) / 4 + pow(n, 4) / 64)

#define rad2deg 57.2957795130823229 /* 180/pi */
#define deg2rad 0.0174532925199433  /* pi/180 */

/* The code below shows that there's no need to go beyond n=4 with double precision
#include <stdio.h>
#include <math.h>

void main()
{
    printf("b = %f\te2 = %f\tn = %f\n", b, e2, n);

    double left_term = a / (1 + n), right_term = 1;

    printf("A = %.16f\n", left_term * right_term);

    right_term += pow(n, 2) / 4;
    printf("A = %.16f\n", left_term * right_term);

    right_term += pow(n, 4) / 64;
    printf("A = %.16f\n", left_term * right_term);

    right_term += pow(n, 6) / 256;
    printf("A = %.16f\n", left_term * right_term);

    right_term += pow(n, 8) / 16384;
    printf("A = %.16f\n", left_term * right_term);
}
*/

static double alpha[] = {0,
    n/2 - 2*pow(n,2)/3 + 5*pow(n,3)/16,
    13*pow(n,2)/48 - 3*pow(n,3)/5,
    61*pow(n,3)/240};
static double beta[] = {0,
    n/2 - 2*pow(n,2)/3 + 37*pow(n,3)/96,
    pow(n,2)/48 + pow(n,3)/15,
    17*pow(n,3)/480};
static double delta[] = {0,
    2*n - 2*pow(n,2)/3 - 2*pow(n,3),
    7*pow(n,2)/3 - 8*pow(n,3)/5,
    56*pow(n,3)/15};

#define k0 0.9996
#define E0 500000 /* in km */

/* other constants */
#define const1 2 * sqrt(n) / (1 + n)
#define const2 (1 - n) / (1 + n)
#define const3 k0 * A

inline int getUTMzone(double lon, double lat, double *lambda0)
{
    if (lat >= 72 && lat < 84) {
	if (lon >= 0 && lon < 9) {
	    *lambda0 = 4.5;
	    return 31;
	}
	if (lon >= 9 && lon < 21) {
	    *lambda0 = 15;
	    return 33;
	}
	if (lon >= 21 && lon < 33) {
	    *lambda0 = 27;
	    return 35;
	}
	if (lon >= 33 && lon < 42) {
	    *lambda0 = 37.5;
	    return 37;
	}
    }
    if (lat >= 56 && lat < 64) {
	if (lon >= 0 && lon < 3) {
	    *lambda0 = 1.5;
	    return 31;
	}
	if (lon >= 3 && lon < 12) {
	    *lambda0 = 7.5;
	    return 32;
	}
    }
    *lambda0 = floor(lon / 6) * 6 + 3;
    return (int) floor((180 + lon)/6) + 1;
}

void lonlat2UTM_(double *lambda, double *phi, double *E, double *N, double *zone, int sample_size)
{
    int i, j;
    double N0, lambda0, sin_phi, diff_lambda, t, xi_prime, zeta_prime, sigma, tau;
    double inter1, inter2, inter3;
    double cos_inter2[4], sin_inter2[4], sinh_inter3[4], cosh_inter3[4];

    for (i = 0; i < sample_size; i++) {
	if (ISNA(lambda[i]) || ISNA(phi[i])) {
	    E[i] = N[i] = zone[i] = NA_REAL;
	    continue;
	}
	if (phi[i] >= 84 || phi[i] < -80) {
	    E[i] = N[i] = zone[i] = NA_REAL;
	    continue;
	}
	N0 = phi[i] > 0 ? 0 : 10000000;
	//	lambda0 = floor(lambda[i] / 6) * 6 + 3;
	zone[i] = getUTMzone(lambda[i], phi[i], &lambda0);
	//	Rprintf("lambda0 = %f\n", lambda0);
	sin_phi = sin(deg2rad * phi[i]);
	diff_lambda = deg2rad * (lambda[i] - lambda0);
	t = sinh(atanh(sin_phi) - const1 * atanh(const1 * sin_phi));
	xi_prime = atan(t / cos(diff_lambda));
	zeta_prime = atanh(sin(diff_lambda) / sqrt(1 + t * t));
	sigma = 1;
	tau = 0;
	for (j = 1; j <= 3; j++) {
	    inter1 = 2 * j * alpha[j];
	    inter2 = 2 * j * xi_prime;
	    cos_inter2[j] = cos(inter2);
	    sin_inter2[j] = sin(inter2);
	    inter3 = 2 * j * zeta_prime;
	    sinh_inter3[j] = sinh(inter3);
	    cosh_inter3[j] = cosh(inter3);
	    sigma += inter1 * cos_inter2[j] * cosh_inter3[j];
	    tau += inter1 * sin(inter2) * sinh_inter3[j];
	}
	E[i] = zeta_prime;
	N[i] = xi_prime;
	for (j = 1; j <= 3; j++) {
	    E[i] += alpha[j] * cos_inter2[j] * sinh_inter3[j];
	    N[i] += alpha[j] * sin_inter2[j] * cosh_inter3[j];
	}
	E[i] *= k0 * A;
	E[i] += E0;
	N[i] *= k0 * A;
	N[i] += N0;
    }
}

SEXP lonlat2UTM_Call(SEXP LON, SEXP LAT)
{
    int sample_size;
    double *E, *N, *Z;
    SEXP res;

    PROTECT(LON = coerceVector(LON, REALSXP));
    PROTECT(LAT = coerceVector(LAT, REALSXP));
    sample_size = LENGTH(LON);
    PROTECT(res = allocMatrix(REALSXP, sample_size, 3));

    E = REAL(res);
    N = E + sample_size;
    Z = N + sample_size;
    lonlat2UTM_(REAL(LON), REAL(LAT), E, N, Z, sample_size);

    /* the version with a data frame as output
    PROTECT(res = allocVector(VECSXP, 3));
    PROTECT(easting = allocVector(REALSXP, sample_size));
    PROTECT(northing = allocVector(REALSXP, sample_size));
    PROTECT(zone = allocVector(INTSXP, sample_size));
    lonlat2UTM_(REAL(LON), REAL(LAT), REAL(easting), REAL(northing), INTEGER(zone), sample_size);
    SET_VECTOR_ELT(res, 0, easting);
    SET_VECTOR_ELT(res, 1, northing);
    SET_VECTOR_ELT(res, 2, zone);
    UNPROTECT(6); */

    UNPROTECT(3);
    return res;
}

inline double get_mid_meridian(int zone, double lat)
{
    if (lat >= 72 && lat < 84) {
	if (zone == 31) return 4.5;
	if (zone == 33) return 15;
	if (zone == 35) return 27;
	if (zone == 37) return 37.5;
    }
    if (lat >= 56 && lat < 64) {
	if (zone == 31) return 1.5;
	if (zone == 32) return 7.5;
    }
    return zone * 6 - 183;
}

void UTM2lonlat_(double *E, double *N, int *zone, int *hemi, double *lambda, double *phi, int sample_size)
{
    int i, j;
    double N0, xi, zeta, xi_prime, zeta_prime, sigma_prime,
	tau_prime, inter1, inter2, inter3, sin_inter1, cosh_inter2,
	cos_inter1, sinh_inter2, chi, tmp, lambda0;

    for (i = 0; i < sample_size; i++) {
	if (ISNA(E[i]) || ISNA(N[i])) {
	    lambda[i] = phi[i] = NA_REAL;
	    continue;
	}
	N0 = hemi[i] > 0 ? 0 : 10000000;
	xi = (N[i] - N0) / const3;
	zeta = (E[i] - E0) / const3;

	xi_prime = xi;
	zeta_prime = zeta;
	sigma_prime = 1;
	tau_prime = 0;

	for (j = 1; j <= 3; j++) {
	    inter1 = 2 * j * xi;
	    inter2 = 2 * j * zeta;
	    inter3 = 2 * j * beta[j];
	    sin_inter1 = sin(inter1);
	    cosh_inter2 = cosh(inter2);
	    cos_inter1 = cos(inter1);
	    sinh_inter2 = sinh(inter2);

	    xi_prime -= beta[j] * sin_inter1 * cosh_inter2;
	    zeta_prime -= beta[j] * cos_inter1 * sinh_inter2;
	    sigma_prime -= inter3 * cos_inter1 * cosh_inter2;
	    tau_prime += inter3 * sin_inter1 * sinh_inter2;
	}

	chi = asin(sin(xi_prime) / cosh(zeta_prime));

	tmp = xi;
	for (j = 1; j <= 3; j++) tmp += delta[j] * sin(2 * j * chi);
	phi[i] = rad2deg * tmp;

	//	lambda0 = zone[i] * 6 - 183;
	lambda0 = get_mid_meridian(zone[i], phi[i]);
	lambda[i] = lambda0 + rad2deg * atan(sinh(zeta_prime) / cos(xi_prime));
    }
}

SEXP UTM2lonlat_Call(SEXP EASTING, SEXP NORTHING, SEXP ZONE, SEXP HEMI)
{
    int sample_size;
    double *lon, *lat;
    SEXP res;

    PROTECT(EASTING = coerceVector(EASTING, REALSXP));
    PROTECT(NORTHING = coerceVector(NORTHING, REALSXP));
    PROTECT(ZONE = coerceVector(ZONE, INTSXP));
    PROTECT(HEMI = coerceVector(HEMI, INTSXP));
    sample_size = LENGTH(EASTING);
    PROTECT(res = allocMatrix(REALSXP, sample_size, 2));

    lon = REAL(res);
    lat = lon + sample_size;
    UTM2lonlat_(REAL(EASTING), REAL(NORTHING), INTEGER(ZONE),
		INTEGER(HEMI), lon, lat, sample_size);

    UNPROTECT(5);
    return res;
}
