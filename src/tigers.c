/* tigers.c       2024-04-30 */

/* Copyright 2023-2024 Emmanuel Paradis */

/* This file is part of the R-package `tigers'. */
/* See the file ../DESCRIPTION for licensing issues. */

#include <R_ext/Rdynload.h>
#include "tigers.h"

/* from sim_polygons_landClasses.c */
int circularIndex(int i, int n)
{
    if (i > n - 1) return i % n;
    while (i < 0) i += n;
    return i;
}

static R_CMethodDef C_entries[] = {
    {"test_Hormann_Floater", (DL_FUNC) &test_Hormann_Floater, 4},
    {"C_specrend", (DL_FUNC) &C_specrend, 6},
    {NULL, NULL, 0}
};

static R_CallMethodDef Call_entries[] = {
    {"convex_hull_C", (DL_FUNC) &convex_hull_C, 2},
    {"convexPolygonOverlap_Call", (DL_FUNC) &convexPolygonOverlap_Call, 2},
    {"InsidePolygon_Call", (DL_FUNC) &InsidePolygon_Call, 2},
    {"triangulate_Call", (DL_FUNC) &triangulate_Call, 2},
    {"haveOverlapTwoPolygons", (DL_FUNC) &haveOverlapTwoPolygons, 2},
    {"redundant_vertices", (DL_FUNC) &redundant_vertices, 3},
    {"area_Call", (DL_FUNC) &area_Call, 1},
    {"rev_Call", (DL_FUNC) &rev_Call, 2},
    {"rev_2cols_Call", (DL_FUNC) &rev_2cols_Call, 2},
    {"RMA_Call", (DL_FUNC) &RMA_Call, 2},
    {"singlePolygon2raster", (DL_FUNC) &singlePolygon2raster, 3},
    {"fast2waytable_Call", (DL_FUNC) &fast2waytable_Call, 4},
    {"ECEF2lonlat_Call", (DL_FUNC) &ECEF2lonlat_Call, 1},
    {"lonlat2UTM_Call", (DL_FUNC) &lonlat2UTM_Call, 2},
    {"UTM2lonlat_Call", (DL_FUNC) &UTM2lonlat_Call, 4},
    {NULL, NULL, 0}
};

static R_FortranMethodDef Fortran_entries[] = {
    {"wltocol", (DL_FUNC) &F77_SUB(wltocol), 4},
    {NULL, NULL, 0}
};

void R_init_tigers(DllInfo *info)
{
    R_registerRoutines(info, C_entries, Call_entries, Fortran_entries, NULL);
    R_useDynamicSymbols(info, FALSE);
}
