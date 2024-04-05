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


/* Define ARKBraidThetaOrdCondsMem, 
stores nonlinear solver data for theta order conditions */

struct _ARKBraidThetaOrdCondsMem
{
  int      order;
  int      nconds;
  int      nstages;

  /* Functions defining the butcher table */
  void (*init)(sunrealtype *th);
  int  (*alt_init)(sunrealtype *th, const sunindextype resets);
  void (*btable_A)(sunrealtype *A, const sunrealtype *th);
  void (*btable_b)(sunrealtype *b, const sunrealtype *th);
  void (*btable_c)(sunrealtype *c, const sunrealtype *th);

  /* Residual and Jacobian functions for order conditions solver */
  void (*res)(sunrealtype *r, const realtype *th, const realtype *rhs);
  void (*jac)(sunrealtype *J, const realtype *th);

  /* Workspace variables for computing order conditions */
  realtype *phi1;
  realtype *phi2;

  /* Workspace variables for nonlinear solver */
  N_Vector rhs;
  N_Vector th0;
  N_Vector thcur;
  N_Vector thcor;
  N_Vector weight;
  N_Vector x;
  SUNMatrix J;
  SUNLinearSolver LS;
};

typedef struct _ARKBraidThetaOrdCondsMem *ARKBraidThetaOrdCondsMem;


/* This stores time-dependent but local information for a particular 
 * level of the MGRIT hierarchy. Data contained here is not communicated
 * between processors.
 */
struct _ARKBraidGridData
{
  int level;
  int num_steps_stored;

  /* Storage for improved initial guesses for stage values */
  N_Vector **stage_zs;   /* stage values stored per step */
  int       *num_stages; /* */

  /* Storage for coarse grid theta method Butcher tables */
  ARKodeButcherTable *coarse_btables;
  
};

typedef struct _ARKBraidGridData *ARKBraidGridData;

int ARKBraidGridData_Create(braid_Int level, braid_Int ntpoints, ARKBraidGridData *grid_ptr);
int ARKBraidGridData_Free(ARKBraidGridData *grid_ptr);


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
  int storage;               /* storage level: 
                               (enables improved initial guess for implicit stages)
                               -1: minimum storage, only C-points 
                                0: full storage on all levels
                           x >= 1: full storage on all levels >= 1 */
  booleantype stage_storage; /* if true, stage values will also be stored for all steps 
                                on all levels where full storage is set */

  /* Functions provided to XBraid (user may override) */
  braid_PtFcnStep        step;    /* take time step       */
  braid_PtFcnInit        init;    /* initialize vector    */
  braid_PtFcnSpatialNorm snorm;   /* norm over space      */
  braid_PtFcnAccess      access;  /* user access function */

  /* Saved return flags */
  braid_Int last_flag_braid;
  int       last_flag_arkode;

  /* Current number of levels in MGRIT hierarchy */
  int  num_levels_alloc;
  int *num_steps_stored; /* number of steps per level which are stored */

  /* Grid data storage */
  ARKBraidGridData *grids;

  /* Storage for implicit stage values on each level */
  int       **num_stages;       /* number of stages stored per step on each level */
  N_Vector ***stage_zs;         /* storage for stage values */

  /* Butcher tables */
  int *num_tables;  /* number of Butcher tables for each level */
  ARKodeButcherTable **coarse_btables; /* array of Butcher tables for each level */
  ARKodeButcherTable fine_btable;

  /* Flags for theta method */
  int flag_refine_downcycle;
  int flag_skip_downcycle;

  /* Fine and coarse grid method orders */
  int order_fine;
  int order_coarse;
  int num_order_conditions;

  /* ARKODE memory for coarse grid */
  ARKodeMem ark_mem_coarse;

  /* Nonlinear solver for theta order conditions */
  SUNNonlinearSolver NLS;
  ARKBraidThetaOrdCondsMem theta_mem;

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
  realtype  tprior;   /* time value at previous C-point */
  realtype  etascale; /* time-step normalization factor */
  realtype *Phi;      /* order condition rhs */
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

int ARKBraidTheta_InitHierarchy(braid_App app, braid_SyncStatus sstatus);

int ARKBraidTheta_StepElemWeights(ARKBraidContent content,
                                  braid_StepStatus status, braid_Vector u);

int _ARKBraidTheta_SetBtable(ARKodeButcherTable B, ARKBraidContent content,
                             realtype* theta);

int _ARKBraidTheta_GetBTable(ARKBraidContent content, braid_StepStatus status,
                             braid_Int level, braid_Int ti, ARKodeButcherTable* B);

int _ARKBraidTheta_AllocCGBtables(ARKBraidContent content);

int ARKBraidTheta_NlsResidual(N_Vector thcor, N_Vector r, void* mem);

int ARKBraidTheta_NlsLSetup(booleantype jbad, booleantype* jcur, void* mem);

int ARKBraidTheta_NlsLSolve(N_Vector b, void* mem);

int ARKBraidTheta_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                              realtype tol, N_Vector ewt, void* mem);

int ARKBraidTheta_NlsSetup(ARKodeMem ark_mem, ARKBraidThetaOrdCondsMem nls_mem,
                           SUNNonlinearSolver* NLS_ptr);

int ARKBraidTheta_NlsSolve(SUNNonlinearSolver NLS, ARKBraidThetaOrdCondsMem nls_mem,
                           realtype* rhs);

int ARKBraidTheta_NlsMem_Create(ARKBraidContent content, ARKBraidThetaOrdCondsMem* nlsmem);

int ARKBraidTheta_NlsMem_Free(ARKBraidThetaOrdCondsMem nlsmem);

#ifdef __cplusplus
}
#endif

#endif
