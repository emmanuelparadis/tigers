/* tigers.h       2024-01-26 */

/* Copyright 2023-2024 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../DESCRIPTION for licensing issues. */

#include <Rinternals.h>

int circularIndex(int i, int n);

SEXP convex_hull_C(SEXP XY, SEXP isClockwise);
SEXP convexPolygonOverlap_Call(SEXP A, SEXP B);
SEXP InsidePolygon_Call(SEXP XY, SEXP P);
void test_Hormann_Floater(double *xy, double *Pxy, double *bary_coords, int *pathlength);
SEXP triangulate_Call(SEXP XY, SEXP METHOD);
SEXP haveOverlapTwoPolygons(SEXP P, SEXP Q);
SEXP redundant_vertices(SEXP POLYGON, SEXP TOL, SEXP CHECK_ONLY);
SEXP area_Call(SEXP XY);
SEXP rev_Call(SEXP x, SEXP copy);
SEXP rev_2cols_Call(SEXP x, SEXP copy);
SEXP RMA_Call(SEXP X, SEXP Y);
SEXP singlePolygon2raster(SEXP XY, SEXP PARS, SEXP raster);

