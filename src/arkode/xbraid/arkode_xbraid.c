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
 * This is the implementation file for the ARKODE + XBraid interface.
 * -------------------------------------------------------------------------- */

#include "arkode/arkode_xbraid.h"

#include "arkode/arkode.h"
#include "arkode_arkstep_impl.h"
#include "arkode_xbraid_impl.h"
#include "sundials/sundials_math.h"

#define ONE RCONST(1.0)

/* -------------------------------
 * Construct, initialize, and free
 * ------------------------------- */

/* Create XBraid app strucutre */
int ARKBraid_Create(void* arkode_mem, braid_App* app)
{
  int flag;
  ARKBraidContent content;
  ARKodeARKStepMem step_mem;

  /* Check input */
  if (arkode_mem == NULL) return SUNBRAID_ILLINPUT;

  /* Create XBraid interface object */
  flag = SUNBraidApp_NewEmpty(app);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Set operations */
  (*app)->ops->getvectmpl  = ARKBraid_GetVecTmpl;
  (*app)->ops->getbufsize  = ARKBraidTheta_GetBufSize;
  (*app)->ops->bufpack     = ARKBraidTheta_BufPack;
  (*app)->ops->bufunpack   = ARKBraidTheta_BufUnpack;
  (*app)->ops->initvecdata = ARKBraidTheta_InitVecData;
  (*app)->ops->freevecdata = ARKBraidTheta_FreeVecData;

  /* Create ARKODE interface content */
  content = NULL;
  content = (ARKBraidContent)malloc(sizeof(*content));
  if (content == NULL)
  {
    (void)SUNBraidApp_FreeEmpty(app);
    return SUNBRAID_ALLOCFAIL;
  }

  /* Initialize content */

  /* Attach ARKODE memory */
  content->ark_mem = (ARKodeMem)arkode_mem;

  /* Interface functions */
  content->step   = ARKBraid_Step;
  content->init   = ARKBraid_Init;
  content->snorm  = SUNBraidVector_SpatialNorm;
  content->access = ARKBraid_Access;

  /* Options */
  content->storage         = -1;
  content->stage_storage   =  0;
  content->use_theta       =  0;
  content->loose_tol_fac   =  1.;   /* default is same tolerance on every level */
  content->loose_fine_rtol =  1e-3; /* for the coarse grid correction to have at least three    */
  content->tight_fine_rtol =  1e-3; /* digits of accuracy, the tau correction and the coarse    */
  content->coarse_rtol     =  1e-3; /* grid equation both need to gain three digits of accuracy */
                                    /* relative to the previous iteration */

  /* Saved return flags */
  content->last_flag_braid  = SUNBRAID_SUCCESS;
  content->last_flag_arkode = SUNBRAID_SUCCESS;

  /* Output time and solution (allocated in access if necessary) */
  content->tout = content->ark_mem->tn;
  content->yout = NULL;

  /* flags for theta method */
  content->flag_refine_downcycle = SUNFALSE;
  content->flag_skip_downcycle   = SUNFALSE;

  /* Newton solver for theta method order conditions */
  content->NLS = NULL;
  content->theta_mem = NULL;
  
  /* Fine and coarse grid method orders */

  if (content->ark_mem->step_mem == NULL) return SUNBRAID_ILLINPUT;
  step_mem = (ARKodeARKStepMem)content->ark_mem->step_mem;

  content->order_fine = step_mem->q;
  content->order_coarse = step_mem->q; // TODO: should this default to q + 1?

  content->num_order_conditions =
    _ARKBraidTheta_GetNumOrderConditions(content->order_fine,
                                         content->order_coarse);

  content->num_levels_alloc = 0;
  content->grids          = NULL;
  content->ark_mem_coarse = NULL;
  content->fine_btable    = NULL;

  /* Attach content */
  (*app)->content = content;

  return SUNBRAID_SUCCESS;
}

/* Initialize XBraid, attach interface functions */
int ARKBraid_BraidInit(MPI_Comm comm_w, MPI_Comm comm_t, realtype tstart,
                       realtype tstop, sunindextype ntime, braid_App app,
                       braid_Core* core)
{
  braid_Int braid_flag;
  ARKBraidContent content;

  /* Check inputs */
  if (comm_w == MPI_COMM_NULL || comm_t == MPI_COMM_NULL || ntime < 2 ||
      app == NULL)
    return SUNBRAID_ILLINPUT;

  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  /* Shortcut to content */
  content = (ARKBraidContent)app->content;

  /* Initialize XBraid */
  braid_flag = braid_Init(comm_w, comm_t, tstart, tstop, ntime, app,
                          content->step, content->init, SUNBraidVector_Clone,
                          SUNBraidVector_Free, SUNBraidVector_Sum,
                          content->snorm, content->access,
                          SUNBraidVector_BufSize, SUNBraidVector_BufPack,
                          SUNBraidVector_BufUnpack, core);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  /* Set sync function */
  braid_flag = braid_SetSync(*core, ARKBraid_Sync);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  /* Set storage */
  braid_flag = braid_SetStorage(*core, content->storage);

  return SUNBRAID_SUCCESS;
}

/* Deallocate XBraid app structure */
int ARKBraid_Free(braid_App* app)
{
  ARKBraidContent content; /* ARKBraid app content  */

  if (*app == NULL) return SUNBRAID_SUCCESS;

  if ((*app)->content != NULL)
  {
    content = (ARKBraidContent)(*app)->content;

    if (content->yout != NULL)
    {
      arkFreeVec(content->ark_mem, &(content->yout));
      content->yout = NULL;
    }

    if (content->ark_mem_coarse != NULL)
    {
      ARKStepFree((void**)&content->ark_mem_coarse); // Free integrator memory
      content->ark_mem_coarse = NULL;
    }

    if (content->fine_btable != NULL)
    {
      ARKodeButcherTable_Free(content->fine_btable);
      content->fine_btable = NULL;
    }
    
    if (content->grids != NULL)
    {
      for (int lvl = 0; lvl < content->num_levels_alloc; lvl++)
      {
        ARKBraidGridData_Free(&content->grids[lvl]);
      }
      free(content->grids);
      content->grids = NULL;
    }

    ARKBraidTheta_NlsMem_Free(content->theta_mem);
    SUNNonlinSolFree(content->NLS);

    free((*app)->content);
    (*app)->content = NULL;
  }
  return SUNBraidApp_FreeEmpty(app);
}

/* ----------------------
 * ARKBraid Set Functions
 * ---------------------- */

int ARKBraid_SetTheta(braid_App app, booleantype theta)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent)app->content;

  content->use_theta = theta;

  return SUNBRAID_SUCCESS;
}

int ARKBraid_SetCoarseOrder(braid_App app, sunindextype order)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent)app->content;

  /* Check order */
  if (order < 0 || order > 4) return SUNBRAID_ILLINPUT;

  /* Restore default or set */
  if (order == 0) content->order_coarse = content->order_fine;
  else content->order_coarse = order;

  /* Get new number of order conditions */
  content->num_order_conditions =
    _ARKBraidTheta_GetNumOrderConditions(content->order_fine,
                                         content->order_coarse);

  return SUNBRAID_SUCCESS;
}

int ARKBraid_SetFullStorage(braid_App app, braid_Int storage, braid_Int stage_storage)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent)app->content;
  content->storage       = storage;
  content->stage_storage = stage_storage;

  return SUNBRAID_SUCCESS;
}

int ARKBraid_SetLooseTolFactor(braid_App app, braid_Real loose_tol_fac)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent)app->content;
  content->loose_tol_fac = loose_tol_fac;

  return SUNBRAID_SUCCESS;
}

int ARKBraid_SetFineRTol(braid_App app, braid_Real tight_rtol, braid_Real loose_rtol)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent)app->content;
  content->tight_fine_rtol = tight_rtol;
  content->loose_fine_rtol = loose_rtol;

  return SUNBRAID_SUCCESS;
}

int ARKBraid_SetCoarseRTol(braid_App app, braid_Real rtol)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent)app->content;
  content->coarse_rtol = rtol;

  return SUNBRAID_SUCCESS;
}

int ARKBraid_SetStepFn(braid_App app, braid_PtFcnStep step)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent)app->content;

  /* Restore default or set function pointer */
  if (step == NULL) content->step = ARKBraid_Step;
  else content->step = step;

  return SUNBRAID_SUCCESS;
}

int ARKBraid_SetInitFn(braid_App app, braid_PtFcnInit init)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent)app->content;

  /* Restore default or set function pointer */
  if (init == NULL) content->init = ARKBraid_Init;
  else content->init = init;

  return SUNBRAID_SUCCESS;
}

int ARKBraid_SetSpatialNormFn(braid_App app, braid_PtFcnSpatialNorm snorm)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent)app->content;

  /* Restore default or set function pointer */
  if (snorm == NULL) content->snorm = SUNBraidVector_SpatialNorm;
  else content->snorm = snorm;

  return SUNBRAID_SUCCESS;
}

int ARKBraid_SetAccessFn(braid_App app, braid_PtFcnAccess access)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent)app->content;

  /* Restore default or set function pointer */
  if (access == NULL) content->access = ARKBraid_Access;
  else content->access = access;

  return SUNBRAID_SUCCESS;
}

/* ----------------------
 * ARKBraid Get Functions
 * ---------------------- */

int ARKBraid_GetVecTmpl(braid_App app, N_Vector* tmpl)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent)app->content;
  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;
  *tmpl = content->ark_mem->yn;
  return SUNBRAID_SUCCESS;
}

int ARKBraid_GetARKStepMem(braid_App app, void** arkode_mem)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent)app->content;
  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;
  *arkode_mem = (void*)content->ark_mem;
  return SUNBRAID_SUCCESS;
}

int ARKBraid_GetUserData(braid_App app, void** user_data)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent)app->content;
  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;
  *user_data = content->ark_mem->user_data;
  return SUNBRAID_SUCCESS;
}

int ARKBraid_GetLastBraidFlag(braid_App app, int* last_flag)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content    = (ARKBraidContent)app->content;
  *last_flag = content->last_flag_braid;
  return SUNBRAID_SUCCESS;
}

int ARKBraid_GetLastARKStepFlag(braid_App app, int* last_flag)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content    = (ARKBraidContent)app->content;
  *last_flag = content->last_flag_arkode;
  return SUNBRAID_SUCCESS;
}

int ARKBraid_GetSolution(braid_App app, realtype* tout, N_Vector yout)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent)app->content;
  if (content->yout == NULL) return SUNBRAID_MEMFAIL;
  *tout = content->tout;
  N_VScale(ONE, content->yout, yout);
  return SUNBRAID_SUCCESS;
}

/* --------------------------
 * XBraid Interface Functions
 * -------------------------- */

/* Take a time step */
int ARKBraid_Step(braid_App app, braid_Vector ustop, braid_Vector fstop,
                  braid_Vector u, braid_StepStatus status)
{
  braid_Int braid_flag;              /* braid function return flag  */
  int ark_flag;                      /* arkode step return flag     */
  int flag;                          /* arkode function return flag */
  int level;                         /* current level               */
  int rfac;                          /* refinement factor           */
  int fixedstep;                     /* flag for fixed step size    */
  int iter;                          /* MGRIT iteration             */
  int iu, il;                        /* lowest and highest time indices on this level    */
  int ti, tir;                       /* time index and relative time index for this step */
  int gotzs, gotustop;
  realtype tstart;                   /* current time                */
  realtype tstop;                    /* evolve to this time         */
  realtype hacc;                     /* accuracy based step size    */
  ARKBraidContent content;           /* ARKBraid app content        */

  N_Vector *z = NULL; /* stage initial guess */
  ARKodeARKStepMem step_mem;
  ARKBraidGridData grid; /* grid data for this step */

  /* Check input */
  if (app == NULL || status == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL || u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Access app content */
  content = (ARKBraidContent)app->content;
  step_mem = (ARKodeARKStepMem)content->ark_mem->step_mem;

  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;

  /* Remember if we are using fixed time-stepping */
  fixedstep = (content->ark_mem->fixedstep);

  /* Get info from XBraid */
  braid_flag = braid_StepStatusGetLevel(status, &level);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);
  braid_flag = braid_StepStatusGetIter(status, &iter);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);
  braid_flag = braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);
  braid_flag = braid_StepStatusGetTIndex(status, &ti);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);
  braid_flag = braid_StepStatusGetTIUL(status, &iu, &il, level);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);
  tir = ti - il + 1; /* relative time index */

  /* Get Grid storage */
  grid = content->grids[level];

  /* Compute theta method order conditions and set Butcher table */
  ARKBraidTheta_StepElemWeights(content, status, u);

  // TODO: Improve spatial accuracy routine
  braid_Real step_atol, step_rtol; /* abs and rel solver tolerances for this step   */
  sunrealtype atol, rtol;          /* tolerances set by user in ARKStepSStolerances */
  atol = content->ark_mem->Sabstol;
  rtol = content->ark_mem->reltol;
  step_rtol = rtol;

  /* Adjust tolerances */
  // braid_GetSpatialAccuracy(status, content->loose_tol_fac * atol, atol, &step_atol);
  ARKBraid_GetSpatialAccuracy(status, content->loose_tol_fac * atol, atol, &step_atol);
  ARKBraid_GetSpatialAccuracy(status, content->loose_tol_fac * rtol, rtol, &step_rtol);
  // if (ti == 0)
  // {
  //   printf("atol: %f, rtol: %f\n", step_atol, step_rtol);
  // }

  /* Get stored stage initial guess if available */
  if (content->stage_storage && grid->stage_zs)
  {
    /* This allows using stored intermediate stage values 
    from previous calls to step at this time index */
    /* ARKODE will allocate this if it is NULL */
    flag = arkSetFullStorage((void*)content->ark_mem, SUNTRUE);
    CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

    flag = arkSetStageZs((void*)content->ark_mem, grid->stage_zs[tir], grid->num_stages[tir]);
    CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);
  }
  gotzs    = (grid->num_stages != NULL && grid->num_stages[tir] > 0);
  gotustop = (u != ustop);

  /* TODO: Store linear solvers/preconditioners either per level or per step */

  /* Set step tolerances */
  flag = ARKStepSStolerances(content->ark_mem, step_rtol, step_atol);
  CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

  /* Turn off error estimation on coarse grids */
  if (!fixedstep && level > 0)
    arkSetFixedStep(content->ark_mem, tstop - tstart);

  /* Get number of linear iterations the step took */
  int liters_start, liters_stop;
  ARKLsMem lsmem;
  lsmem = (ARKLsMem)content->ark_mem->step_getlinmem((void*)content->ark_mem);

  liters_start = lsmem->nli;

  flag = ARKBraid_TakeStep((void*)(content->ark_mem), tstart, tstop, u->y, ustop->y, &ark_flag);

  CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

  liters_stop = lsmem->nli;
  // printf("Level: %d, ti: %d, rtol=%e, atol=%e, z: %d, liters: %d\n", level, ti, step_rtol, step_atol, gotzs, liters_stop-liters_start);

  /* Restore fixedstep value */
  if (!fixedstep && level > 0)
    /* Setting zero here turns adaptivity back on */
    arkSetFixedStep(content->ark_mem, ZERO);

  /* Restore step tolerances */
  flag = ARKStepSStolerances(content->ark_mem, rtol, atol);
  CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

  /* Retrieve updated intermediate stages */
  if (grid->stage_zs)
  {
    arkGetStageZs(content->ark_mem, &grid->stage_zs[tir], &grid->num_stages[tir]);
    CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);
  }

  /* Refine grid (XBraid will ignore if refinement is disabled) */

  /* Compute refinement factor */
  if (iter > 0 && braid_StepStatusAcceptsRFactor(status))
  {
    // TODO: only refine when residual sufficiently small?
    /* Default to no refinement */
    rfac = 1;

    /* The step failed due to solver failure or too much error */
    if (ark_flag != 0)
    {
      /* Get the suggested step size. The rfac value is given by ETACF on a
         solver failure and limited by ETAMIN on an error test failure */
      flag = ARKStepGetCurrentStep((void*)(content->ark_mem), &hacc);
      CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

      /* Set the refinement factor */
      rfac = (int)(SUNRceil((tstop - tstart) / hacc));

      /* Limit the refinement factor */
      rfac = (rfac < 1) ? 1 : rfac;
    }

    /* set the refinement factor */
    // if (rfac > 1)
    //   printf("Refine: ti: %d, rfactor=%d\n", ti, rfac);
    braid_flag = braid_StepStatusSetRFactor(status, rfac);
    CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);
  }

  return SUNBRAID_SUCCESS;
}

/* Create and initialize vectors */
int ARKBraid_Init(braid_App app, realtype t, braid_Vector* u_ptr)
{
  int flag;                   /* return flag          */
  N_Vector y;                 /* output N_Vector      */
  ARKBraidThetaVecData vdata; /* output ARKBraid vector data 
                                 (used for Theta method) */
  ARKBraidContent content;    /* ARKBraid app content */

  /* Check input */
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  /* Access app content */
  content = (ARKBraidContent)app->content;

  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;

  /* Create new NVector */
  y = NULL;
  if (!arkAllocVec(content->ark_mem, content->ark_mem->yn, &y))
    return SUNBRAID_ALLOCFAIL;

  /* Create new XBraid vector */
  flag = SUNBraidVector_New(app, y, u_ptr);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Set initial solution at all time points */
  N_VScale(ONE, content->ark_mem->yn, y);

  return SUNBRAID_SUCCESS;
}

/* User access function */
int ARKBraid_Access(braid_App app, braid_Vector u, braid_AccessStatus astatus)
{
  braid_Int braid_flag;    /* braid return flag    */
  braid_Int done;          /* braid finished flag  */
  braid_Int ntpoints;      /* num pts on fine grid */
  braid_Int idx;           /* time index for u     */
  braid_Real time;         /* time value for u     */
  ARKBraidContent content; /* ARKBraid app content */

  /* Check input */
  if (app == NULL || u == NULL || astatus == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL || u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Access app content */
  content = (ARKBraidContent)app->content;

  if (content->ark_mem) return SUNBRAID_MEMFAIL;

  /* Check if XBraid is done with the current simulation */
  braid_flag = braid_AccessStatusGetDone(astatus, &done);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  if (done)
  {
    /* Get global number of points on the fine grid */
    braid_flag = braid_AccessStatusGetNTPoints(astatus, &ntpoints);
    CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

    /* Get the time index for the vector u */
    braid_flag = braid_AccessStatusGetTIndex(astatus, &idx);
    CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

    /* Get the time for the vector u */
    braid_flag = braid_AccessStatusGetT(astatus, &time);
    CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

    /* Check if this is the last time point */
    if (idx == ntpoints - 1)
    {
      /* Allocate yout if necessary */
      if (content->yout == NULL)
      {
        if (!arkAllocVec(content->ark_mem, content->ark_mem->yn, &(content->yout)))
          return SUNBRAID_ALLOCFAIL;
      }

      /* Save solution for output to user */
      content->tout = time;
      N_VScale(ONE, u->y, content->yout);
    }
  }

  return SUNBRAID_SUCCESS;
}

int ARKBraid_Sync(braid_App app, braid_SyncStatus sstatus)
{
  int flag;           /* return flag          */
  int caller;         /* XBraid calling function */
  int level, nlevels; /* XBraid level, number of levels */
  int skip;           /* XBraid skip option */
  int iter;           /* XBraid iteration */
  int proc;           /* ID of current proc */
  int iu, il;         /* upper and lower time-indices on this proc */
  int ntpoints;       /* XBraid number of time points */
  ARKBraidContent content; /* ARKBraid app content */
  ARKodeARKStepMem step_mem; /* ARKStep memory */
  ARKBraidGridData grid;

  /* Check input */

  if (app == NULL || sstatus == NULL) return SUNBRAID_ILLINPUT;

  /* Access app content */
  content = (ARKBraidContent)app->content;

  /* Store the fine-grid butcher table */
  if (content->fine_btable == NULL)
  {
    /* Get the ARKStep memory */
    step_mem = (ARKodeARKStepMem)content->ark_mem->step_mem;

    if ((step_mem->Be == NULL) && (step_mem->Bi == NULL))
    {
      /* Initialize arkode_mem so that there is a valid Butcher table */
      flag = arkStep_SetButcherTables(content->ark_mem);
      CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);
    }

    /* TODO: support for ImEx methods */
    /* Get the fine-grid butcher table */
    if (step_mem->Be != NULL)
      content->fine_btable = ARKodeButcherTable_Copy(step_mem->Be);
    else if (step_mem->Bi != NULL)
      content->fine_btable = ARKodeButcherTable_Copy(step_mem->Bi);
    else return SUNBRAID_ILLINPUT;
  }

  /* Get information from XBraid status */
  flag = braid_SyncStatusGetCallingFunction(sstatus, &caller);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_SyncStatusGetSkip(sstatus, &skip);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_SyncStatusGetLevel(sstatus, &level);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_SyncStatusGetNLevels(sstatus, &nlevels);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_SyncStatusGetIter(sstatus, &iter);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_SyncStatusGetTIUL(sstatus, &iu, &il, 0);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  
  /* proc 0 will own the zero time index (il == 0) */
  proc = il;

  /* Check if this is a new hierarchy */
  if (caller == braid_ASCaller_Drive_AfterInit ||
      caller == braid_ASCaller_FRefine_AfterInitHier)
  {
    /* Destroy old grid data */
    if (content->grids)
    {
      for (int lvl = 0; lvl < content->num_levels_alloc; lvl++)
      {
        ARKBraidGridData_Free(&content->grids[lvl]);
      }
    }

    /* Store new number of levels */
    content->num_levels_alloc = nlevels;

    /* Allocate new grid data objects */
    content->grids = (ARKBraidGridData*)realloc(content->grids, nlevels * sizeof(ARKBraidGridData));
    if (content->grids == NULL) return SUNBRAID_ALLOCFAIL;

    /* Initialize grid data storage */
    for (int lvl = 0; lvl < nlevels; lvl++)
    {
      content->grids[lvl] = NULL;

      flag = braid_SyncStatusGetTIUL(sstatus, &iu, &il, lvl);
      CHECK_BRAID_RETURN(content->last_flag_braid, flag);
      ntpoints = iu - il + 1;

      flag = ARKBraidGridData_Create(lvl, ntpoints, &content->grids[lvl]);
      if (flag != SUNBRAID_SUCCESS) return flag;

      /* Allocate stage storage */
      if (content->stage_storage)
      {
        grid = content->grids[lvl];
        grid->stage_zs = (N_Vector**)calloc(ntpoints, sizeof(N_Vector*));
        if (grid->stage_zs == NULL) return SUNBRAID_ALLOCFAIL;
        grid->num_stages = (sunindextype*)calloc(ntpoints, sizeof(sunindextype));
        if (grid->num_stages == NULL) return SUNBRAID_ALLOCFAIL;
      }
    }

    /* Check that theta method needs computing */
    if (nlevels > 1 && content->use_theta)
    {
      /* Set flag if the first down-cycle will be skipped */
      content->flag_skip_downcycle = (skip && caller == braid_ASCaller_Drive_AfterInit);

      /* Initialize coarse grid theta method Butcher tables */
      flag = ARKBraidTheta_InitHierarchy(app, sstatus);
      if (flag != SUNBRAID_SUCCESS) return flag;
    }
  }
  return SUNBRAID_SUCCESS;
}

/* -----------------
 * Utility Functions
 * ----------------- */

/* Force a single step with ARKEvolve */
int ARKBraid_TakeStep(void* arkode_mem, realtype tstart, realtype tstop,
                      N_Vector y, N_Vector ystop, int* ark_flag)
{
  int flag;      /* generic return flag      */
  int tmp_flag;  /* evolve return flag       */
  realtype tret; /* return time              */

  /* Check inputs */
  if (arkode_mem == NULL) return ARK_MEM_NULL;
  if (y == NULL) return ARK_ILL_INPUT;

  /* Reset ARKStep state */
  flag = ARKStepReset(arkode_mem, tstart, y);
  if (flag != ARK_SUCCESS) return flag;

  /* Set the time step size */
  flag = ARKStepSetInitStep(arkode_mem, tstop - tstart);
  if (flag != ARK_SUCCESS) return flag;

  /* Ignore temporal error test result and force step to pass */
  flag = arkSetForcePass(arkode_mem, SUNTRUE);
  if (flag != ARK_SUCCESS) return flag;

  /* Set improved initial guess, if available */
  if (ystop != y)
  {
    /* This allows linear/cubic hermite interpolation */
    flag = ARKStepSetStepGuess(arkode_mem, tstop, ystop);
    if (flag != ARK_SUCCESS && flag != ARK_ILL_INPUT) return flag;
  }

  /* Take step, check flag below */
  tmp_flag = ARKStepEvolve(arkode_mem, tstop, y, &tret, ARK_ONE_STEP);

  /* Re-enable temporal error test check */
  flag = arkSetForcePass(arkode_mem, SUNFALSE);
  if (flag != ARK_SUCCESS) return flag;

  flag = arkSetFullStorage(arkode_mem, SUNFALSE);
  if (flag != ARK_SUCCESS) return flag;

  /* Check if evolve call failed */
  if (tmp_flag < 0)
  {
    *ark_flag = STEP_FAILED;
    return ARK_SUCCESS;
  }

  /* Check if temporal error test failed */
  flag = arkGetLastKFlag(arkode_mem, &tmp_flag);
  if (flag != ARK_SUCCESS) return flag;

  if (tmp_flag > 0)
  {
    *ark_flag = STEP_ADAPT;
    return ARK_SUCCESS;
  }

  /* Step was successful and passed the error test */
  *ark_flag = STEP_SUCCESS;
  return ARK_SUCCESS;
}

int ARKBraidGridData_Create(braid_Int level, braid_Int ntpoints, ARKBraidGridData *grid_ptr)
{
  ARKBraidGridData grid;
  braid_Int        iu, il;

  grid = (ARKBraidGridData)malloc(sizeof(*grid));
  if (grid == NULL) return SUNBRAID_ALLOCFAIL;

  grid->level            = level;
  grid->num_steps_stored = ntpoints;
  grid->num_stages       = NULL;
  grid->stage_zs         = NULL;
  grid->coarse_btables   = NULL;

  *grid_ptr = grid;

  return SUNBRAID_SUCCESS;
}

/* Free allocations for stage values and coarse-grid Butcher tables */
int ARKBraidGridData_Free(ARKBraidGridData *grid_ptr)
{
  ARKBraidGridData grid = *grid_ptr; 

  /* Check input */
  if (grid == NULL) return SUNBRAID_SUCCESS;

  if (grid->coarse_btables != NULL)
  {
      for (int i = 0; i < grid->num_steps_stored; i++)
      {
        ARKodeButcherTable_Free(grid->coarse_btables[i]);
      }
      free(grid->coarse_btables);
      grid->coarse_btables = NULL;
  }

  /* Destroy stored stage values */
  if (grid->stage_zs)
  {
      for (int step = 0; step < grid->num_steps_stored; step++)
      {
        if (grid->stage_zs[step] == NULL) continue;
        for (int is = 0; is < grid->num_stages[step]; is++)
        {
          N_VDestroy(grid->stage_zs[step][is]);
        }
        free(grid->stage_zs[step]);
        grid->stage_zs[step] = NULL;
      }
    free(grid->stage_zs);
    free(grid->num_stages);
    grid->stage_zs   = NULL;
    grid->num_stages = NULL;
  }
  free(grid);
  *grid_ptr = NULL;

  return SUNBRAID_SUCCESS;
}

int ARKBraid_GetSpatialAccuracy(braid_StepStatus  status, braid_Real loose_tol,
                             braid_Real tight_tol, braid_Real *tol_ptr)
{
   braid_Int nrequest   = 2;
   braid_Real stol, tol, rnorm, rnorm0, old_fine_tolx;
   braid_Int level;
   braid_Real l_rnorm, l_ltol, l_ttol, l_tol;
   braid_Real *rnorms = (braid_Real *) malloc( 2*sizeof(braid_Real) );

   braid_StepStatusGetTol(status, &tol);
   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetOldFineTolx(status, &old_fine_tolx);

   /* Get the first and then the current residual norms */
   rnorms[0] = -1.0; rnorms[1] = -1.0;
   braid_StepStatusGetRNorms(status, &nrequest, rnorms);
   if((rnorms[0] == -1.0) && (rnorms[1] != -1.0)){
      rnorm0 = rnorms[1];
   }
   else{
      rnorm0 = rnorms[0];
   }
   nrequest = -2;
   braid_StepStatusGetRNorms(status, &nrequest, rnorms);
   if((rnorms[1] == -1.0) && (rnorms[0] != -1.0)){
      rnorm = rnorms[0];
   }
   else{
      rnorm = rnorms[1];
   }


   if (level == 0)
   {
      /* Always return tight tolerance here to ensure refinement is done correctly */
      *tol_ptr = tight_tol;

      free(rnorms);
      return SUNBRAID_SUCCESS;
   }
   else if ( (nrequest == 0) || (rnorm0 == -1.0) )
   {
      /* Always return the loose tolerance, if there is no residual history yet 
       * (this is the first Braid iteration with skip turned on) */
      *tol_ptr = loose_tol;
   }
   else
   {
      /* Else, do a variable tolerance for the fine grid */
      // l_rnorm = -log10(rnorm / rnorm0);
      // l_tol   = -log10(tol / rnorm0);
      l_rnorm = -log10(rnorm);
      l_tol   = -log10(tol);
      l_ltol  = -log10(loose_tol);
      l_ttol  = -log10(tight_tol);

      if ( l_rnorm >= (7.0/8.0)*l_tol )
      {
         /* Close to convergence, return tight_tol */
         *tol_ptr = tight_tol;
      }
      else
      {
         /* Do linear interpolation between loose_tol and tight_tol (but with respect to log10) */
         stol = (l_rnorm / l_tol) * (l_ttol - l_ltol) + l_ltol;
         *tol_ptr = pow(10, -stol);

         /* The fine grid tolerance MUST never decrease */
         if ( ((*tol_ptr) > old_fine_tolx) && (old_fine_tolx > 0) )
         {
            *tol_ptr = old_fine_tolx;
         }
      }
   }

   /* Store this tolerance */
   braid_StepStatusSetOldFineTolx(status, (*tol_ptr));

   free(rnorms);
   /* printf( "lev: %d, accuracy: %1.2e, nreq: %d, rnorm: %1.2e, rnorm0: %1.2e, loose: %1.2e, tight: %1.2e, old: %1.2e, braid_tol: %1.2e \n", level, *tol_ptr, nrequest, rnorm, rnorm0, loose_tol, tight_tol, old_fine_tolx, tol); */
   return SUNBRAID_SUCCESS;
}