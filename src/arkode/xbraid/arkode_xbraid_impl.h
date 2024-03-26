/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation header file for the ARKode + XBraid interface.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKSTEP_XBRAID_IMP_H
#define _ARKSTEP_XBRAID_IMP_H

#include "sundials/sundials_types.h"
#include "sundials/sundials_nonlinearsolver.h"
#include "arkode_impl.h"
#include "braid.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* --------------
 * Utility macros
 * -------------- */


#define CHECK_BRAID_RETURN(last_flag, flag)                             \
  do { (last_flag) = (flag); if ((flag) != 0) return SUNBRAID_BRAIDFAIL; } while(0)

#define CHECK_ARKODE_RETURN(last_flag, flag)                            \
  do { (last_flag) = (flag); if ((flag) != 0) return SUNBRAID_SUNFAIL; } while(0)


/* --------------------------
 * SUNBraid private constants
 * -------------------------- */


/* TakeSetup step result flags */
#define STEP_FAILED   -1
#define STEP_SUCCESS   0
#define STEP_ADAPT     1


/* ------------------------
 * Nonlinear solver memory 
 * ------------------------ */


/* Define ARKBraidNlsData, stores nonlinear solver data */
struct _ARKBraidNlsMem
{
  int      order;
  int      nconds;
  N_Vector rhs;
  N_Vector th0;
  N_Vector thcur;
  N_Vector thcor;
  N_Vector weight;
  N_Vector x;
  SUNMatrix J;
  SUNLinearSolver LS;
  void (*res)(sunrealtype *r, const realtype *th, const realtype *rhs);
  void (*jac)(sunrealtype *J, const realtype *th);
  void (*init)(sunrealtype *th);
  int (*alt_init)(sunrealtype *th, const sunindextype resets);
};

typedef struct _ARKBraidNlsMem *ARKBraidNlsMem;


/* ------------------------------
 * ARKBraid app structure content
 * ------------------------------ */


/* Define SUNBraidApp content */
struct _ARKBraidContent
{
  /* ARKODE memory structure */
  ARKodeMem ark_mem;

  /* Options */
  int rfac_limit;  /* refinement factor limit           */
  int rfac_fail;   /* refinement factor for failed step */

  /* Functions provided to XBraid (user may override) */
  braid_PtFcnStep        step;    /* take time step       */
  braid_PtFcnInit        init;    /* initialize vector    */
  braid_PtFcnSpatialNorm snorm;   /* norm over space      */
  braid_PtFcnAccess      access;  /* user access function */

  /* Saved return flags */
  braid_Int last_flag_braid;
  int       last_flag_arkode;

  /* flags for theta method */
  int flag_refine_downcycle;
  int flag_skip_downcycle;

  /* Fine and coarse grid method orders */
  int order_fine;
  int order_coarse;
  int num_order_conditions;

  /* ARKODE memory for coarse grid */
  ARKodeMem ark_mem_coarse;

  /* Butcher tables */
  int num_levels;  /* number of levels in MGRIT hierarchy */
  int *num_tables; /* number of Butcher tables for each level */
  ARKodeButcherTable **coarse_btables; /* array of Butcher tables for each level */
  ARKodeButcherTable fine_btable;

  /* Nonlinear solver */
  SUNNonlinearSolver NLS;
  ARKBraidNlsMem NLS_mem;

  /* Output time and state */
  realtype tout;
  N_Vector yout;
};

typedef struct _ARKBraidContent *ARKBraidContent;


/* ---------------------
 * ARKBraid vector data
 * --------------------- */


/* Define ARKBraidThetaVecData content */
struct _ARKBraidThetaVecData
{
  /* Store time value at previous C-point */
  realtype tprior;

  /* Store time-step normalization factor */
  realtype etascale;

  /* Store order condition rhs */
  realtype *Phi;
};

typedef struct _ARKBraidThetaVecData *ARKBraidThetaVecData;


/* ---------------------
 * AKBraidTheta
 * --------------------- */


booleantype _ARKBraid_IsCPoint(int tindex, int cfactor);

int _ARKBraidTheta_GetNumOrderConditions(int fine_order, int coarse_order);

int ARKBraidTheta_InitVecData(braid_App app, void** vdata_ptr);

int ARKBraidTheta_FreeVecData(braid_App app, void* vdata_ptr);

int ARKBraidTheta_GetBufSize(braid_App app, braid_Int* size_ptr);

int ARKBraidTheta_BufPack(braid_App app, void* buffer, void* vdata_ptr);

int ARKBraidTheta_BufUnpack(braid_App app, void* buffer, void** vdata_ptr);

int ARKBraidTheta_Sync(braid_App app, braid_SyncStatus sstatus);

int ARKBraidTheta_StepElemWeights(ARKBraidContent content,
                                  braid_StepStatus status, braid_Vector u);

int _ARKBraidTheta_SetBtable(ARKodeButcherTable B, ARKBraidContent content,
                             realtype* theta);

int _ARKBraidTheta_GetBTable(ARKBraidContent content, braid_StepStatus status,
                             braid_Int level, braid_Int ti, ARKodeButcherTable* B);

int _ARKBraidTheta_AllocCGBtables(ARKBraidContent content,
                                         braid_SyncStatus sstatus);

int _ARKBraidTheta_FreeCGBtables(ARKBraidContent content);

int ARKBraidTheta_NlsResidual(N_Vector thcor, N_Vector r, void* mem);

int ARKBraidTheta_NlsLSetup(booleantype jbad, booleantype* jcur, void* mem);

int ARKBraidTheta_NlsLSolve(N_Vector b, void* mem);

int ARKBraidTheta_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                              realtype tol, N_Vector ewt, void* mem);

int ARKBraidTheta_NlsSetup(ARKodeMem ark_mem, ARKBraidNlsMem nls_mem,
                           SUNNonlinearSolver* NLS_ptr);

int ARKBraidTheta_NlsSolve(SUNNonlinearSolver NLS, ARKBraidNlsMem nls_mem,
                           realtype* rhs);

int ARKBraidTheta_NlsMem_Create(ARKBraidContent content, ARKBraidNlsMem* nlsmem);

int ARKBraidTheta_NlsMem_Free(ARKBraidNlsMem nlsmem);

#ifdef __cplusplus
}
#endif

#endif
