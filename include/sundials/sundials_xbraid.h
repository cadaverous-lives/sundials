/* --------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * --------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * --------------------------------------------------------------------------
 * This is the header file for SUNDIALS + XBraid interface base class and
 * NVector interface.
 * -------------------------------------------------------------------------- */

#ifndef _SUNDIALS_XBRAID_H
#define _SUNDIALS_XBRAID_H

#include "sundials/sundials_types.h"
#include "sundials/sundials_nvector.h"
#include "braid.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* -----------------------
 * XBraid vector structure
 * ----------------------- */


struct _braid_Vector_struct
{
  N_Vector  y; /* the solution vector */
  void *vdata; /* auxilliary time-dependent data, handled by SUNBraidOps */
};

/* Pointer to vector wrapper (same as braid_Vector) */
typedef struct _braid_Vector_struct *SUNBraidVector;


/* -----------------------------
 * XBraid ops and app structures
 * ----------------------------- */


/* Structure containing function pointers to operations */
struct _SUNBraidOps
{
  int (*getvectmpl)(braid_App app, N_Vector *tmpl);
  // 1. Add new “methods” to the SUNBraidOps structure (vtable)
  int (*getbufsize)(braid_App app, braid_Int *size_ptr);
  int (*bufpack)(braid_App app, void* buffer, void *vdata_ptr);
  int (*bufunpack)(braid_App app, void* buffer, void **vdata_ptr);
  int (*freevecdata)(braid_App app, void *vdata);
};

/* Pointer to operations structure */
typedef struct _SUNBraidOps *SUNBraidOps;


/* Define XBraid App structure */
struct _braid_App_struct
{
  void        *content;
  SUNBraidOps ops;
};

/* Pointer to the interface object (same as braid_App) */
typedef struct _braid_App_struct *SUNBraidApp;


/* -----------------------
 * SUNBraid app operations
 * ----------------------- */


SUNDIALS_EXPORT int SUNBraidApp_NewEmpty(braid_App *app);

SUNDIALS_EXPORT int SUNBraidApp_FreeEmpty(braid_App *app);

SUNDIALS_EXPORT int SUNBraidApp_GetVecTmpl(braid_App app, N_Vector *tmpl);

// 2. Add generic version of these functions (declare in include/sundials/sundials_xbraid.h and define in src/sundials/sundials_xbraid.c) e.g.,
SUNDIALS_EXPORT int SUNBraidApp_GetBufSize(braid_App app, braid_Int *size_ptr);

SUNDIALS_EXPORT int SUNBraidApp_BufPack(braid_App app, void* buffer, void* data_ptr);

SUNDIALS_EXPORT int SUNBraidApp_BufUnpack(braid_App app, void* buffer, void** data_ptr);

SUNDIALS_EXPORT int SUNBraidApp_FreeVecData(braid_App app, void* data_ptr);

/* -------------------------
 * SUNBraid vector functions
 * ------------------------- */


SUNDIALS_EXPORT int SUNBraidVector_New(N_Vector y, SUNBraidVector *u);

SUNDIALS_EXPORT int SUNBraidVector_GetNVector(SUNBraidVector u, N_Vector *y);

SUNDIALS_EXPORT int SUNBraidVector_Clone(braid_App app, braid_Vector u,
                                         braid_Vector *v_ptr);

SUNDIALS_EXPORT int SUNBraidVector_Free(braid_App app, braid_Vector u);

SUNDIALS_EXPORT int SUNBraidVector_Sum(braid_App app,
                                       braid_Real alpha, braid_Vector x,
                                       braid_Real beta, braid_Vector y);

SUNDIALS_EXPORT int SUNBraidVector_SpatialNorm(braid_App app, braid_Vector u,
                                               braid_Real *norm_ptr);

SUNDIALS_EXPORT int SUNBraidVector_BufSize(braid_App app, braid_Int *size_ptr,
                                           braid_BufferStatus bstatus);

SUNDIALS_EXPORT int SUNBraidVector_BufPack(braid_App app, braid_Vector u,
                                           void *buffer,
                                           braid_BufferStatus bstatus);

SUNDIALS_EXPORT int SUNBraidVector_BufUnpack(braid_App app, void *buffer,
                                             braid_Vector *u_ptr,
                                             braid_BufferStatus bstatus);


/* ----------------------
 * SUNBraid return values
 * ---------------------- */


#define SUNBRAID_SUCCESS       0  /* call/operation was successful  */

#define SUNBRAID_ALLOCFAIL    -1  /* a memory allocation failed     */
#define SUNBRAID_MEMFAIL      -2  /* a memory access fail           */
#define SUNBRAID_OPNULL       -3  /* the SUNBraid operation is NULL */
#define SUNBRAID_ILLINPUT     -4  /* an invalid input was provided  */
#define SUNBRAID_BRAIDFAIL    -5  /* an XBraid function failed      */
#define SUNBRAID_SUNFAIL      -6  /* a SUNDIALS function failed     */


#ifdef __cplusplus
}
#endif

#endif
