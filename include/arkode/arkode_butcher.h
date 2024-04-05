/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for ARKode Butcher table structures.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_BUTCHER_H
#define _ARKODE_BUTCHER_H

#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*---------------------------------------------------------------
  Types : struct ARKodeButcherTableMem, ARKodeButcherTable
  ---------------------------------------------------------------*/
struct ARKodeButcherTableMem {

  int q;           /* method order of accuracy       */
  int p;           /* embedding order of accuracy    */
  int stages;      /* number of stages               */
  realtype **A;    /* Butcher table coefficients     */
  realtype *c;     /* canopy node coefficients       */
  realtype *b;     /* root node coefficients         */
  realtype *d;     /* embedding coefficients         */

};


typedef _SUNDIALS_STRUCT_ ARKodeButcherTableMem *ARKodeButcherTable;


/* Utility routines to allocate/free/output Butcher table structures */
SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_Alloc(int stages,
                                                            booleantype embedded);
SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_Create(int s, int q,
                                                             int p,
                                                             realtype *c,
                                                             realtype *A,
                                                             realtype *b,
                                                             realtype *d);

SUNDIALS_EXPORT void ARKodeButcherTable_Fill(ARKodeButcherTable B,
                                             realtype *c, realtype *A,
                                                          realtype *b, 
                                                          realtype *d);

SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_Copy(ARKodeButcherTable B);
SUNDIALS_EXPORT void ARKodeButcherTable_Space(ARKodeButcherTable B,
                                              sunindextype *liw,
                                              sunindextype *lrw);
SUNDIALS_EXPORT void ARKodeButcherTable_Free(ARKodeButcherTable B);
SUNDIALS_EXPORT void ARKodeButcherTable_Write(ARKodeButcherTable B,
                                              FILE *outfile);

SUNDIALS_EXPORT int ARKodeButcherTable_CheckOrder(ARKodeButcherTable B, int *q,
                                                  int *p, FILE *outfile);
SUNDIALS_EXPORT int ARKodeButcherTable_CheckARKOrder(ARKodeButcherTable B1,
                                                     ARKodeButcherTable B2,
                                                     int *q, int *p,
                                                     FILE *outfile);

/* utility functions for computing elementary weights */
int _phi_order2(realtype *b, realtype *c, int s, realtype *bc);

int _phi_order3a(realtype *b, realtype *c1, realtype *c2, int s, realtype *bcc);
int _phi_order3b(realtype *b, realtype **A, realtype *c, int s, realtype *bAc);

int _phi_order4a(realtype *b, realtype *c1, realtype *c2, realtype *c3, int s, realtype *bccc);
int _phi_order4b(realtype *b, realtype *c1, realtype **A, realtype *c2, int s, realtype *bcAc);
int _phi_order4c(realtype *b, realtype **A, realtype *c1, realtype *c2, int s, realtype *bAcc);
int _phi_order4d(realtype *b, realtype **A1, realtype **A2, realtype *c, int s, realtype*bAAc);

int _phi_order5a(realtype *b, realtype *c1, realtype *c2, realtype *c3, realtype *c4, int s, realtype *bcccc);
int _phi_order5b(realtype *b, realtype *c1, realtype *c2, realtype **A, realtype *c3, int s, realtype *bccAc);
int _phi_order5c(realtype *b, realtype **A1, realtype *c1, realtype **A2, realtype *c2, int s, realtype *bAcAc);
int _phi_order5d(realtype *b, realtype *c1, realtype **A, realtype *c2, realtype *c3, int s, realtype *bcAcc);
int _phi_order5e(realtype *b, realtype **A, realtype *c1, realtype *c2, realtype *c3, int s, realtype *bAccc);
int _phi_order5f(realtype *b, realtype *c1, realtype **A1, realtype **A2, realtype *c2, int s, realtype *bcAAc);
int _phi_order5g(realtype *b, realtype **A1, realtype *c1, realtype **A2, realtype *c2, int s, realtype *bAcAc);
int _phi_order5h(realtype *b, realtype **A1, realtype **A2, realtype *c1, realtype *c2, int s, realtype *bAAcc);
int _phi_order5i(realtype *b, realtype **A1, realtype **A2, realtype **A3, realtype *c, int s, realtype *bAAAc);

#ifdef __cplusplus
}
#endif

#endif
