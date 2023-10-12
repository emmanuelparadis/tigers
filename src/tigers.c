/* tigers.c       2023-10-06 */

/* Copyright 2023 Emmanuel Paradis */

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
    {NULL, NULL, 0}
};

static R_CallMethodDef Call_entries[] = {
    {"convex_hull_C", (DL_FUNC) &convex_hull_C, 2},
    {"convexPolygonOverlap_Call", (DL_FUNC) &convexPolygonOverlap_Call, 2},
    {"InsidePolygon_Call", (DL_FUNC) &InsidePolygon_Call, 2},
    {"triangulate_Call", (DL_FUNC) &triangulate_Call, 2},
    {"haveOverlapTwoPolygons", (DL_FUNC) &haveOverlapTwoPolygons, 2},
    {"redundant_vertices", (DL_FUNC) &redundant_vertices, 3},
    {"area_Call", (DL_FUNC) &area_Call, 2},
    {"rev_Call", (DL_FUNC) &rev_Call, 2},
    {"rev_2cols_Call", (DL_FUNC) &rev_2cols_Call, 2},
    {"RMA_Call", (DL_FUNC) &RMA_Call, 2},
    {"singlePolygon2raster", (DL_FUNC) &singlePolygon2raster, 3},
    {NULL, NULL, 0}
};

void R_init_tigers(DllInfo *info)
{
    R_registerRoutines(info, C_entries, Call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
