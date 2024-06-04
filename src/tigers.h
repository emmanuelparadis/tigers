/* tigers.h       2024-05-31 */

/* Copyright 2023-2024 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../DESCRIPTION for licensing issues. */

#include <R.h>
#include <Rinternals.h>

int circularIndex(int i, int n);

SEXP convex_hull_C(SEXP XY, SEXP isClockwise);
SEXP convexPolygonOverlap_Call(SEXP A, SEXP B);
SEXP InsidePolygon_Call(SEXP XY, SEXP P);
void test_Hormann_Floater(double *xy, double *Pxy, double *bary_coords, int *pathlength);
SEXP triangulate_Call(SEXP XY, SEXP METHOD);
SEXP haveOverlapTwoPolygons(SEXP P, SEXP Q);
SEXP redundant_vertices(SEXP POLYGON, SEXP TOL, SEXP PARS);
SEXP area_Call(SEXP XY);
SEXP rev_Call(SEXP x, SEXP copy);
SEXP rev_2cols_Call(SEXP x, SEXP copy);
SEXP RMA_Call(SEXP X, SEXP Y);
SEXP singlePolygon2raster(SEXP XY, SEXP PARS, SEXP raster);
void C_specrend(double *spec_intens, double *R, double *G, double *B, int *approx, int *color_system);
int F77_NAME(wltocol)(double *WAVELEN, int *N, double *GAMMA, double *RES);
SEXP fast2waytable_Call(SEXP X, SEXP Y, SEXP NCAT, SEXP TT);
SEXP ECEF2lonlat_Call(SEXP XYZ);
SEXP lonlat2UTM_Call(SEXP LON, SEXP LAT);
SEXP UTM2lonlat_Call(SEXP EASTING, SEXP NORTHING, SEXP ZONE, SEXP HEMI);
SEXP inters_seg_seg(SEXP SEGI, SEXP SEGII);
SEXP inters_seg_arc(SEXP SEG, SEXP ARC);

