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


/* Define ARKBraidVecData content */
struct _ARKBraidVecData
{
  /* Store time value at previous C-point */
  realtype tprior;

  /* Store time-step normalization factor */
  realtype etascale;

  /* Store order condition rhs */
  realtype *Phi;
};

typedef struct _ARKBraidVecData *ARKBraidVecData;




#ifdef __cplusplus
}
#endif

#endif
