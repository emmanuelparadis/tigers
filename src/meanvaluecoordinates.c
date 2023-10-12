/* meanvaluecoordinates.c       2023-10-06 */

/* Copyright 2023 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../DESCRIPTION for licensing issues. */

#include <R.h>

/* i is looped between 0 and n - 1 */
/* #define i_minus_one(i) i == 0 ? n - 1 : i - 1 */
#define i_plus_one(i) (i + 1) % n

/*
  Kai Hormann and Michael S. Floater (2006) Mean value coordinates for
  arbitrary planar polygons. ACM Transactions on Graphics 25(4):
  1424-1441. <doi:10.1145/1183287.1183295>

  The polygon must be *open*.

  The algo seems to work with both clockwise and counterclockwise
  polygons.

  xy: double array of size 2*n
  Pxy: double array of size 2
  bary_coords: double array of size n
  n: number of vertices (= number of edges)

  Return value:
  0: Pxy is inside the polygon
  1: Pxy coincides with a vertex of the polygon
  2: Pxy is on an edge of the polygon
*/
int mean_value_coordinates_Hormann_Floater(double *xy, double *Pxy, double *bary_coords, int n)
{
    int i, ip;
    double *dxy, ri, rip, rim, Ai, Aim, Di, Dim, denom, eps = 1e-8, w, SW;

    /* signed distances in and x and y ($s_i$ in the paper) */
    dxy = (double*)R_alloc(2 * n, sizeof(double));

    for (i = 0; i < n; i++) {
        dxy[i] = xy[i] - Pxy[0];
        dxy[i + n] = xy[i + n] - Pxy[1];
    }

    memset(bary_coords, 0, n * sizeof(double));

    /* check and handle special cases */
    for (i = 0; i < n; i++) {
        ip = i_plus_one(i);
        ri = sqrt(pow(dxy[i], 2) + pow(dxy[i + n], 2)); /* ||s_i|| */
        if (ri <= eps) {
            bary_coords[i] = 1;
            return 1;
        } else {
	    Ai = (dxy[i] * dxy[ip + n] - dxy[ip] * dxy[i + n]) / 2;
	    Di = dxy[i] * dxy[ip] + dxy[i + n] * dxy[ip + n];
	    if (Ai == 0 && Di < 0) {
		rip = sqrt(pow(dxy[ip], 2) + pow(dxy[ip + n], 2)); /* ||s_i+|| */
		denom = ri + rip;
		bary_coords[i]  = rip / denom;
		bary_coords[ip]  = ri / denom;
		return 2;
	    }
        }
    }

    /* Note: it is *not* faster to store the Ai's and ri's in arrays,
       even for large values of 'n' */

    SW = 0;
    /* start with im = -1 (= n - 1), i = 0, ip = 1: */
    ri = sqrt(pow(dxy[0], 2) + pow(dxy[n], 2));
    rim = sqrt(pow(dxy[n - 1], 2) + pow(dxy[2 * n - 1], 2));
    Aim = (dxy[n - 1] * dxy[n] - dxy[0] * dxy[2 * n - 1]) / 2;
    Dim = dxy[n - 1] * dxy[0] + dxy[2 * n - 1] * dxy[n];
    i = 0;
    ip = 1;
    /* 'im' is not used below because if the update 'i' -> 'i-'
       so only to use 'i' and 'i+' (ip)*/
    while (i < n) {
	w = 0;
        rip = sqrt(pow(dxy[ip], 2) + pow(dxy[ip + n], 2));
	if (Aim) w += (rim - Dim / ri) / Aim;
        Ai = (dxy[i] * dxy[ip + n] - dxy[ip] * dxy[i + n]) / 2;
	Di = dxy[i] * dxy[ip] + dxy[i + n] * dxy[ip + n];
	if (Ai) w += (rip - Di / ri) / Ai;
	SW += w;
	bary_coords[i] = w;
	/* update for next loop*/
	rim = ri;
	ri = rip;
	Aim = Ai;
	Dim = Di;
	i++;
	ip++;
	if (ip >= n) ip %= n;
    }

    for (i = 0; i < n; i++) bary_coords[i] /= SW;

    return 0;
}

void test_Hormann_Floater(double *xy, double *Pxy, double *bary_coords, int *pathlength)
{
    int n = pathlength[0]/* , res */;

    /* res =  */mean_value_coordinates_Hormann_Floater(xy, Pxy, bary_coords, n);
    /* Rprintf("res = %d\n", res); */
}

