/* --------------------------------------------------------------------------
 * Programmer(s): David A. Vargas @ Univ. of New Mexico
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

#include "arkode/arkode.h"
#include "arkode/arkode_xbraid.h"
#include "arkode_arkstep_impl.h"
#include "arkode_xbraid_impl.h"
#include "arkode_xbraid_theta_methods.c"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"

/* -----------------
 * Helper functions
 * ----------------- */
booleantype _ARKBraid_IsCPoint(int tindex, int cfactor)
{
  if (tindex % cfactor == 0) return SUNTRUE;
  return SUNFALSE;
}

int _ARKBraidTheta_GetNumOrderConditions(int fine_order, int coarse_order)
{
  int num_order_conditions = 0;

  /* Check if theta method needs computing */
  if (coarse_order <= fine_order) return num_order_conditions;

  /* Compute number of order conditions */
  if (coarse_order >= 2) num_order_conditions += 1;
  if (coarse_order >= 3) num_order_conditions += 2;
  if (coarse_order >= 4) num_order_conditions += 4;
  if (coarse_order >= 5) num_order_conditions += 9;

  return num_order_conditions;
}

int _ARKBraidTheta_CloneARKodeMem(braid_App app, ARKodeMem ark_mem,
                                  ARKodeMem* ark_mem_clone)
{
  int flag;
  N_Vector ytmp;
  ARKodeARKStepMem step_mem, step_mem_new;
  ARKodeMem ark_mem_new = NULL;

  /* Check input */
  if (ark_mem == NULL || ark_mem->step_mem == NULL) return SUNBRAID_ILLINPUT;

  /* Get template vector */
  ARKBraid_GetVecTmpl(app, &ytmp);

  /* Create new ARKodeMem */
  ark_mem_new = ARKStepCreate(step_mem->fe, step_mem->fi, 0., ytmp,
                              ark_mem->sunctx);

  /* Copy ops */
  ark_mem_new->step_attachlinsol   = ark_mem->step_attachlinsol;
  ark_mem_new->step_attachmasssol  = ark_mem->step_attachmasssol;
  ark_mem_new->step_disablelsetup  = ark_mem->step_disablelsetup;
  ark_mem_new->step_disablemsetup  = ark_mem->step_disablemsetup;
  ark_mem_new->step_getlinmem      = ark_mem->step_getlinmem;
  ark_mem_new->step_getmassmem     = ark_mem->step_getmassmem;
  ark_mem_new->step_getimplicitrhs = ark_mem->step_getimplicitrhs;
  ark_mem_new->step_mmult          = ark_mem->step_mmult;
  ark_mem_new->step_getgammas      = ark_mem->step_getgammas;
  ark_mem_new->step_init           = ark_mem->step_init;
  ark_mem_new->step_fullrhs        = ark_mem->step_fullrhs;
  ark_mem_new->step                = ark_mem->step;

  /* Copy over settings */
  ark_mem_new->user_data = ark_mem->user_data;
  ark_mem_new->itol      = ark_mem->itol;
  ark_mem_new->ritol     = ark_mem->ritol;
  ark_mem_new->reltol    = ark_mem->reltol;
  ark_mem_new->Sabstol   = ark_mem->Sabstol;
  if (ark_mem->Vabstol != NULL)
    ark_mem_new->Vabstol = N_VClone(ark_mem->Vabstol);
  ark_mem_new->atolmin0 = ark_mem->atolmin0;
  ark_mem_new->SRabstol = ark_mem->SRabstol;
  if (ark_mem->VRabstol != NULL)
    ark_mem_new->VRabstol = N_VClone(ark_mem->VRabstol);
  ark_mem_new->Ratolmin0      = ark_mem->Ratolmin0;
  ark_mem_new->user_efun      = ark_mem->user_efun;
  ark_mem_new->efun           = ark_mem->efun;
  ark_mem_new->e_data         = ark_mem->e_data;
  ark_mem_new->user_rfun      = ark_mem->user_rfun;
  ark_mem_new->rfun           = ark_mem->rfun;
  ark_mem_new->r_data         = ark_mem->r_data;
  ark_mem_new->constraintsSet = ark_mem->constraintsSet;
  if (ark_mem->constraints != NULL)
    ark_mem_new->constraints = N_VClone(ark_mem->constraints);

  ark_mem_new->mxstep         = ark_mem->mxstep;
  ark_mem_new->mxhnil         = ark_mem->mxhnil;
  ark_mem_new->maxconstrfails = ark_mem->maxconstrfails;
  ark_mem_new->maxnef         = ark_mem->maxnef;
  ark_mem_new->maxncf         = ark_mem->maxncf;

  ark_mem_new->report = ark_mem->report;
  ark_mem_new->diagfp = ark_mem->diagfp;

  ark_mem_new->ProcessStep  = ark_mem->ProcessStep;
  ark_mem_new->ps_data      = ark_mem->ps_data;
  ark_mem_new->ProcessStage = ark_mem->ProcessStage;

  /* Copy step_mem */
  step_mem     = (ARKodeARKStepMem)ark_mem->step_mem;
  step_mem_new = (ARKodeARKStepMem)ark_mem_new->step_mem;

  /* Linear Solver Data */
  step_mem_new->linit       = step_mem->linit;
  step_mem_new->lsetup      = step_mem->lsetup;
  step_mem_new->lsolve      = step_mem->lsolve;
  step_mem_new->lfree       = step_mem->lfree;
  step_mem_new->lmem        = step_mem->lmem;
  step_mem_new->lsolve_type = step_mem->lsolve_type;

  /* Mass matrix solver data */
  step_mem_new->minit       = step_mem->minit;
  step_mem_new->msetup      = step_mem->msetup;
  step_mem_new->mmult       = step_mem->mmult;
  step_mem_new->msolve      = step_mem->msolve;
  step_mem_new->mfree       = step_mem->mfree;
  step_mem_new->mass_mem    = step_mem->mass_mem;
  step_mem_new->mass_type   = step_mem->mass_type;
  step_mem_new->msolve_type = step_mem->msolve_type;

  /* Set output pointer */
  *ark_mem_clone = ark_mem_new;

  return SUNBRAID_SUCCESS;
}

/*-------------------------------------
 * VecData structure
 *-------------------------------------*/

int ARKBraidTheta_InitVecData(braid_App app, void** vdata_ptr)
{
  int flag;
  ARKBraidThetaVecData vdata;

  /* Check input */
  if (app == NULL) return SUNBRAID_ILLINPUT;

  /* Access content */
  ARKBraidContent content = (ARKBraidContent)app->content;
  if (content == NULL) return SUNBRAID_ILLINPUT;

  /* Allocate vector data struct */
  vdata = NULL;
  vdata = (ARKBraidThetaVecData)malloc(sizeof(struct _ARKBraidThetaVecData));
  if (vdata == NULL) return SUNBRAID_ALLOCFAIL;

  /* Create vector data */
  vdata->tprior   = RCONST(0.);
  vdata->etascale = RCONST(1.);
  vdata->Phi      = NULL;

  /* Allocate for order conditions (Phi) */
  int num_conditions = _ARKBraidTheta_GetNumOrderConditions(content->order_fine,
                                                            content->order_coarse);
  if (num_conditions > 0)
  {
    vdata->Phi = (realtype*)malloc(num_conditions * sizeof(realtype));
    if (vdata->Phi == NULL) return SUNBRAID_ALLOCFAIL;
    for (size_t i = 0; i < num_conditions; i++) {
      vdata->Phi[i] = 0.;
    }
  }

  /* Attach vector data */
  *vdata_ptr = (void*)vdata;

  return SUNBRAID_SUCCESS;
}

int ARKBraidTheta_FreeVecData(braid_App app, void* vdata_ptr)
{
  ARKBraidThetaVecData vdata;

  if (vdata_ptr == NULL) return SUNBRAID_SUCCESS;

  vdata = (ARKBraidThetaVecData)vdata_ptr;

  if (vdata->Phi != NULL)
  {
    free(vdata->Phi);
    vdata->Phi = NULL;
  }
  free(vdata);
  vdata_ptr = NULL;

  return SUNBRAID_SUCCESS;
}

/* return the buffer space needed, this will get added to the buffer space
 * needed for the vector data in SUNBraidVector_BufSize */
int ARKBraidTheta_GetBufSize(braid_App app, braid_Int* size_ptr)
{
  ARKBraidContent content; /* ARKBraid app content */

  /* Check input */
  if (app == NULL || size_ptr == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  /* Access app content */
  content = (ARKBraidContent)app->content;

  /* Set buffer size */
  *size_ptr = (content->num_order_conditions + 2) * sizeof(realtype);

  return SUNBRAID_SUCCESS;
}

/* pack/unpack data in/out of the buffer, SUNBraidVector_BufPack will offset the
 * input pointer past the vector data */
int ARKBraidTheta_BufPack(braid_App app, void* buffer, void* vdata_ptr)
{
  /* Check input */
  if (app == NULL || buffer == NULL || vdata_ptr == NULL)
    return SUNBRAID_ILLINPUT;

  /* Access vector data */
  ARKBraidThetaVecData vdata = (ARKBraidThetaVecData)vdata_ptr;

  /* Access app content */
  ARKBraidContent content = (ARKBraidContent)app->content;

  /* Copy data into buffer */
  realtype* buf = (realtype*)buffer;
  buf[0]        = vdata->tprior;
  buf[1]        = vdata->etascale;
  for (int i = 0; i < content->num_order_conditions; i++) {
    buf[i + 2] = vdata->Phi[i];
  }

  return SUNBRAID_SUCCESS;
}

int ARKBraidTheta_BufUnpack(braid_App app, void* buffer, void** vdata_ptr)
{
  int flag;                   /* return flag                  */
  ARKBraidThetaVecData vdata; /* ARKBraid vector data, output */

  /* Check input */
  if (app == NULL || buffer == NULL || vdata_ptr == NULL)
    return SUNBRAID_ILLINPUT;

  /* Access app content */
  ARKBraidContent content = (ARKBraidContent)app->content;

  /* Allocate vdata */
  flag = ARKBraidTheta_InitVecData(app, vdata_ptr);
  if (flag != SUNBRAID_SUCCESS) return flag;

  vdata = (ARKBraidThetaVecData)*vdata_ptr;

  /* Copy data from buffer */
  realtype* buf   = (realtype*)buffer;
  vdata->tprior   = buf[0];
  vdata->etascale = buf[1];
  for (int i = 2; i < content->num_order_conditions; i++) {
    vdata->Phi[i] = buf[i];
  }

  /* Return vector data */
  *vdata_ptr = vdata;

  return SUNBRAID_SUCCESS;
}

/*--------------------------------
 * XBraid wrapper functions
 *--------------------------------*/

/* Sync is used by theta method to allocate memory for coarse grid butcher
 * tables */
int ARKBraidTheta_Sync(braid_App app, braid_SyncStatus sstatus)
{
  int flag;                /* return flag          */
  ARKBraidContent content; /* ARKBraid app content */
  int caller;              /* XBraid calling function */
  int level, nlevels;      /* XBraid level, number of levels */
  int iter;                /* XBraid iteration */
  int iu, il;   /* XBraid index of lowest and highest time indices owned by this
                   processor */
  int ncpoints; /* XBraid number of coarse grid points */
  ARKodeARKStepMem step_mem; /* ARKStep memory */

  /* Check input */

  if (app == NULL || sstatus == NULL) return SUNBRAID_ILLINPUT;

  /* Access app content */
  content = (ARKBraidContent)app->content;

  /* Need to store the fine-grid butcher table */
  if (content->fine_btable == NULL)
  {
    // ARKodeButcherTable B;
    // /* Get default coarse btable */
    // B = ARKodeButcherTable_Alloc(content->order_coarse, &B);
    // if (B == NULL) return SUNBRAID_ALLOCFAIL;
    // flag = _ARKBraidTheta_SetBtable(B, content, NULL);
    // if (flag != SUNBRAID_SUCCESS) return flag;

    /* Trick ARKODE into allocating enough memory for the coarse method */
    // flag = ARKStepSetTables(content->ark_mem, content->order_coarse,
    //                         content->order_coarse, B, NULL);
    // CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);
    // flag = arkInitialSetup(content->ark_mem, ONE);
    // CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

    /* Initialize arkode_mem so that there is a valid Butcher table */
    flag = arkStep_SetButcherTables(content->ark_mem);
    CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

    /* Get the ARKStep memory */
    step_mem = (ARKodeARKStepMem)content->ark_mem->step_mem;

    /* Get the fine-grid butcher table */
    if (step_mem->Be != NULL)
      content->fine_btable = ARKodeButcherTable_Copy(step_mem->Be);
    else if (step_mem->Bi != NULL)
      content->fine_btable = ARKodeButcherTable_Copy(step_mem->Bi);
    else return SUNBRAID_ILLINPUT;

    // ARKodeButcherTable_Free(B);
  }

  /* Check that theta method needs computing */
  if (content->order_fine >= content->order_coarse) return SUNBRAID_SUCCESS;

  // /* Create the coarse grid ARKodeMem */
  // flag = _ARKBraidTheta_CloneARKodeMem(app, content->ark_mem,
  // &content->ark_mem_coarse); if (flag != SUNBRAID_SUCCESS) return flag;

  // /* Set the coarse grid to always use implicit methods */
  // flag = ARKStepSetImplicit(content->ark_mem_coarse);
  // CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

  // /* Set coarse grid to not use temporal adaptivity? */
  // flag = ARKStepSetFixedStep((void*)content->ark_mem_coarse, RCONST(0.1));
  // CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

  /* Get information from XBraid status */
  flag = braid_SyncStatusGetCallingFunction(sstatus, &caller);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_SyncStatusGetLevel(sstatus, &level);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_SyncStatusGetNLevels(sstatus, &nlevels);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_SyncStatusGetIter(sstatus, &iter);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);

  /* Check if this is a new hierarchy */
  if (caller == braid_ASCaller_Drive_AfterInit ||
      caller == braid_ASCaller_FRefine_AfterInitHier)
  {
    _ARKBraidTheta_FreeCGBtables(content);

    /* Allocate memory for coarse grid butcher tables */
    flag = _ARKBraidTheta_AllocCGBtables(content, sstatus);
    if (flag != SUNBRAID_SUCCESS) return flag;

    /* Setup nonlinear solver */
    flag = ARKBraidTheta_NlsMem_Create(content, &content->NLS_mem);
    if (flag != SUNBRAID_SUCCESS) return flag;
    flag = ARKBraidTheta_NlsSetup(content->ark_mem, content->NLS_mem,
                                  &content->NLS);
    if (flag != SUNBRAID_SUCCESS) return flag;

    /* Turn on computation of coarse grid Butcher tables */
    if (nlevels > 1)
    {
      content->flag_refine_downcycle = SUNTRUE;
      printf("ARKBraid: Computing coarse grid Butcher tables\n");
    }
  }

  /* make sure this is still computed even when we skip the first downcycle */
  /* TODO: implement this correctly */
  // if (caller == braid_ASCaller_Drive_TopCycle && iter == 0)
  // {
  //   content->flag_refine_downcycle = SUNTRUE;
  //   printf("ARKBraid: Skipping first downcycle, but still computing coarse "
  //          "grid Butcher tables\n");
  // }

  return SUNBRAID_SUCCESS;
}

int ARKBraidTheta_StepElemWeights(ARKBraidContent content,
                                  braid_StepStatus status, braid_Vector u)
{
  int flag;                   /* return flag from XBraid calls */
  int caller;                 /* XBraid calling function */
  int level;                  /* level of hierarchy */
  int ti;                     /* time index */
  booleantype frelax;         /* F-relaxation flag */
  ARKBraidThetaVecData vdata; /* vector data (contains elementary weights) */
  ARKodeButcherTable B;       /* Butcher table */

  /* Check input */
  if (content == NULL || content->ark_mem == NULL || status == NULL || u == NULL)
    return SUNBRAID_ILLINPUT;

  flag = braid_StepStatusGetCallingFunction(status, &caller);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);

  /* Turn off computation of elementary weights once upcycle starts */
  if (content->flag_refine_downcycle && caller == braid_ASCaller_FInterp)
  {
    printf("Turning off elementary weight computation\n");
    content->flag_refine_downcycle = SUNFALSE;
  }

  flag = braid_StepStatusGetLevel(status, &level);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_StepStatusGetTIndex(status, &ti);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);

  /* Get the butcher table for this step */
  _ARKBraidTheta_GetBTable(content, status, level, ti, &B);

  /* F-relax */
  frelax = (caller == braid_ASCaller_FRestrict ||
            caller == braid_ASCaller_FASResidual);
  if (content->flag_refine_downcycle && frelax && u->vdata != NULL)
  {
    int cfactor;            /* coarsening factor */
    realtype tstart, tstop; /* start and stop time of current step */
    realtype eta, eta_pr;   /* normalized time step sizes       */
    realtype eta_c;         /* predicted normalized coarse step */

    /* Get information from XBraid status */
    flag = braid_StepStatusGetCFactor(status, &cfactor, level);
    CHECK_BRAID_RETURN(content->last_flag_braid, flag);
    flag = braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
    CHECK_BRAID_RETURN(content->last_flag_braid, flag);

    /* Access vector data */
    vdata = (ARKBraidThetaVecData)u->vdata;

    /* First step of f-interval */
    if (_ARKBraid_IsCPoint(ti, cfactor))
    {
      /* Reset elementary weights */
      for (int i = 0; i < content->num_order_conditions; i++)
        vdata->Phi[i] = RCONST(0.0);

      vdata->tprior   = tstart;
      vdata->etascale = cfactor * (tstop - tstart);
    }

    /* etascale is a prediction of the total length of the coarse interval
     * we normalize by this for numerical stability (scalars are O(1/m) rather
     * than O(h))
     */
    eta    = (tstop - tstart) / vdata->etascale;
    eta_pr = (tstart - vdata->tprior) / vdata->etascale;
    realtype *phi_step, *phi_pr;

    /* Allocate temporary arrays TODO: preallocate in app->content? */
    phi_step = (realtype*)malloc(content->num_order_conditions * sizeof(realtype));
    if (phi_step == NULL) return SUNBRAID_MEMFAIL;
    phi_pr = (realtype*)malloc(content->num_order_conditions * sizeof(realtype));
    if (phi_pr == NULL) return SUNBRAID_MEMFAIL;

    /* Compute elementary weights */

    /* Temporarily copy prior weights */
    for (int i = 0; i < content->num_order_conditions; i++)
      phi_pr[i] = vdata->Phi[i];

    /* These recurrences explained in (upcoming) theta methods paper... */

    /* Second order */
    _phi_order2(B->b, B->c, B->stages, &phi_step[0]);
    vdata->Phi[0] += eta * eta * phi_step[0] + eta * eta_pr;

    /* Third order */
    if (content->order_coarse >= 3)
    {
      _phi_order3a(B->b, B->c, B->c, B->stages, &phi_step[1]);
      _phi_order3b(B->b, B->A, B->c, B->stages, &phi_step[2]);
      vdata->Phi[1] += SUNRpowerI(eta, 3) * phi_step[1] +
                       RCONST(2.) * eta * eta * eta_pr * phi_step[0] +
                       eta * eta_pr * eta_pr;
      vdata->Phi[2] += SUNRpowerI(eta, 3) * phi_step[2] + eta * phi_pr[0] +
                       eta * eta * eta_pr * phi_step[0];
    }

    /* Last step of f-interval */
    if (_ARKBraid_IsCPoint(ti + 1, cfactor))
    {
      int iu, il; /* highest and lowest time indices on this processor */
      int tic;    /* time index on next coarsest level */

      /* Normalize weights */
      eta_c = (tstop - vdata->tprior) / vdata->etascale;

      /* Second order */
      vdata->Phi[0] /= SUNRpowerI(eta_c, 2);

      /* Third order */
      if (content->order_coarse >= 3)
      {
        vdata->Phi[1] /= SUNRpowerI(eta_c, 3);
        vdata->Phi[2] /= SUNRpowerI(eta_c, 3);
      }

      /* Solve order conditions (solution in NLS_mem->thcur)*/
      ARKBraidTheta_NlsSolve(content->NLS, content->NLS_mem, vdata->Phi);

      /* Set new Butcher table */

      /* Computing on level for level + 1 */
      flag = braid_StepStatusGetTIUL(status, &iu, &il, level + 1);
      CHECK_BRAID_RETURN(content->last_flag_braid, flag);

      tic = (ti+1) / cfactor - il;

      // printf("(This proc il=%d, iu=%d) SetBtable: level=%d, tic=%d\n", il, iu, level, tic); 

      _ARKBraidTheta_SetBtable(content->coarse_btables[level + 1][tic], content,
                               NV_DATA_S(content->NLS_mem->thcur));
    }

    free(phi_step);
    free(phi_pr);
  }

  /* Set the Butcher table for this step TODO: support for explicit methods */
  // if (app->explicit && level == 0)
  //   flag = ARKStepSetTables(ark_mem, B->q, B->p, NULL, B);
  // else
  // flag = ARKStepSetTables(content->ark_mem_coarse, B->q, B->p, B, NULL);
  flag = ARKStepSetTables(content->ark_mem, B->q, B->p, B, NULL);
  CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

  return SUNBRAID_SUCCESS;
}

/*-------------------------
 * Butcher table functions
 *-------------------------*/

/* Set coarse grid Butcher table given theta parameters */
int _ARKBraidTheta_SetBtable(ARKodeButcherTable B, ARKBraidContent content,
                             realtype* theta)
{
  /* TODO: get/store somewhere the actual number of stages */
  int s;             /* Number of stages                                 */
  booleantype guess; /* If we don't know theta, provide an initial guess */
  realtype* A;       /* Butcher table entries                            */

  /* Check inputs */
  if (B == NULL || content == NULL) return SUNBRAID_ILLINPUT;
  guess = (theta == NULL);

  s = B->stages;

  A = (realtype*)malloc(s * s * sizeof(realtype));
  if (A == NULL) return SUNBRAID_ALLOCFAIL;

  if (guess)
  {
    theta = (realtype*)malloc(content->num_order_conditions * sizeof(realtype));
    if (theta == NULL)
    {
      free(A);
      return SUNBRAID_ALLOCFAIL;
    }
  }

  /* Set Butcher table */
  /* TODO: store the guess and btable functions in the content struct? */
  switch (content->order_coarse)
  {
  case 2:
    if (guess) _theta_sdirk2_guess(theta);
    _theta_sdirk2_btable_A(A, theta);
    _theta_sdirk2_btable_b(B->b, theta);
    _theta_sdirk2_btable_c(B->c, theta);
    B->q = 1;
    B->p = 1;
    break;
  case 3:
    if (guess) _theta_sdirk3_guess(theta);
    _theta_sdirk3_btable_A(A, theta);
    _theta_sdirk3_btable_b(B->b, theta);
    _theta_sdirk3_btable_c(B->c, theta);
    B->q = 2;
    B->p = 2;
    break;
  default:
    free(A);
    if (guess) free(theta);
    return SUNBRAID_ILLINPUT;
  }

  for (int i = 0; i < s; i++)
  {
    for (int j = 0; j < s; j++) { B->A[i][j] = A[i * s + j]; }
  }

  free(A);
  if (guess) free(theta);
  return SUNBRAID_SUCCESS;
}

/* returns the Butcher table to be used for index ti on the given level */
int _ARKBraidTheta_GetBTable(ARKBraidContent content, braid_StepStatus status,
                             braid_Int level, braid_Int ti, ARKodeButcherTable* B)
{
  int flag;   /* return flag */
  int iu, il; /* upper and lower indices for the current level */
  int tir; /* time index, relative to first time point owned by current proc */

  /* Get the coarse (relative) time index */
  flag = braid_StepStatusGetTIUL(status, &iu, &il, level);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  tir = ti - il + 1;

  // printf("GetBTable: il=%d, iu=%d, level=%d, ti=%d\n", il, iu, level, ti);
  if (level == 0 || content->order_coarse <= content->order_fine ||
      content->coarse_btables == NULL || content->coarse_btables[level] == NULL ||
      content->coarse_btables[level][tir] == NULL)
  {
    /* Use the fine table */
    *B = content->fine_btable;
  }
  else
  {
    /* Use the coarse table */
    *B = content->coarse_btables[level][tir];
  }
  return SUNBRAID_SUCCESS;
}

int _ARKBraidTheta_AllocCGBtables(ARKBraidContent content,
                                  braid_SyncStatus sstatus)
{
  int flag;     /* Return flag */
  int ntpoints; /* Number of time-points on this level */
  int ns;       /* Number of theta method stages */
  int nlevels;  /* Number of levels in MGRIT hierarchy */
  int iu, il;   /* Highest and lowest time points owned by this proc */

  /* Check input */
  if (content == NULL || sstatus == NULL) return SUNBRAID_ILLINPUT;

  /* Query XBraid status */
  flag = braid_SyncStatusGetNLevels(sstatus, &nlevels);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);

  content->num_levels = nlevels;

  content->coarse_btables = (ARKodeButcherTable**)malloc((nlevels - 1) * sizeof(ARKodeButcherTable*));
  if (content->coarse_btables == NULL) return SUNBRAID_ALLOCFAIL;
  content->coarse_btables -= 1; /* offset array by one, since we are storing
                                   nothing for lvl 0 */

  content->num_tables = (int*)malloc((nlevels - 1) * sizeof(int));
  if (content->num_tables == NULL) return SUNBRAID_ALLOCFAIL;
  content->num_tables -= 1; /* offset array by one, since we are storing nothing
                               for lvl 0 */

  /* TODO: This is obviously not going to generalize past order 3 */
  ns = content->order_coarse;

  /* Preallocate Butcher tables for each coarse level */

  for (int lvl = 1; lvl < nlevels; lvl++)
  {
    content->coarse_btables[lvl] = NULL;

    flag = braid_SyncStatusGetTIUL(sstatus, &iu, &il, lvl);
    CHECK_BRAID_RETURN(content->last_flag_braid, flag);

    ntpoints                 = iu - il + 1;
    content->num_tables[lvl] = ntpoints;
    content->coarse_btables[lvl] = (ARKodeButcherTable*)malloc(ntpoints * sizeof(ARKodeButcherTable));
    if (content->coarse_btables[lvl] == NULL) return SUNBRAID_ALLOCFAIL;

    for (int ti = 0; ti < ntpoints; ti++)
    {
      content->coarse_btables[lvl][ti] = NULL;
      content->coarse_btables[lvl][ti] = ARKodeButcherTable_Alloc(ns, 0);
      if (content->coarse_btables[lvl][ti] == NULL) return SUNBRAID_ALLOCFAIL;
      /* Set a sensible default */
      _ARKBraidTheta_SetBtable(content->coarse_btables[lvl][ti], content, NULL);
    }
  }

  return SUNBRAID_SUCCESS;
}

int _ARKBraidTheta_FreeCGBtables(ARKBraidContent content)
{
  if (content->coarse_btables != NULL)
  {
    for (int level = 1; level < content->num_levels; level++)
    {
      for (int table = 0; table < content->num_tables[level]; table++)
      {
        if (content->coarse_btables[level][table] != NULL)
          ARKodeButcherTable_Free(content->coarse_btables[level][table]);
      }
      free(content->coarse_btables[level]);
      content->coarse_btables[level] = NULL;
    }
    free(content->num_tables+1);
    free(content->coarse_btables+1);
    content->num_tables     = NULL;
    content->coarse_btables = NULL;
  }

  return SUNBRAID_SUCCESS;
}

/* ---------------------------------
 * Solve dense nonlinear systems
 * --------------------------------- */

#define NLS_TOL RCONST(1.0e-10) /* nonlinear solver tolerance */

int ARKBraidTheta_NlsResidual(N_Vector thcor, N_Vector r, void* mem)
{
  ARKBraidNlsMem nlsmem;

  /* Check inputs */
  if (mem == NULL) return SUNBRAID_ILLINPUT;
  nlsmem = (ARKBraidNlsMem)mem;

  /* Update state using given correction */
  N_VLinearSum(ONE, nlsmem->th0, ONE, thcor, nlsmem->thcur);

  /* Compute residual */
  nlsmem->res(NV_DATA_S(r), NV_DATA_S(nlsmem->thcur), NV_DATA_S(nlsmem->rhs));

  return SUNBRAID_SUCCESS;
}

/* Setup linear solve for Newton iteration */
int ARKBraidTheta_NlsLSetup(booleantype jbad, booleantype* jcur, void* mem)
{
  int flag; /* return flag */
  ARKBraidNlsMem nls_mem;

  if (mem == NULL) return SUNBRAID_ILLINPUT;
  nls_mem = (ARKBraidNlsMem)mem;

  /* Compute Jacobian */
  nls_mem->jac(SUNDenseMatrix_Data(nls_mem->J), NV_DATA_S(nls_mem->thcur));

  /* Update Jacobian status */
  *jcur = SUNTRUE;

  /* setup linear solver */
  flag = SUNLinSolSetup(nls_mem->LS, nls_mem->J);

  return SUNBRAID_SUCCESS;
}

int ARKBraidTheta_NlsLSolve(N_Vector b, void* mem)
{
  int flag; /* return flag */
  ARKBraidNlsMem nls_mem;

  if (mem == NULL) return SUNBRAID_ILLINPUT;
  nls_mem = (ARKBraidNlsMem)mem;

  /* solve linear system */
  flag = SUNLinSolSolve(nls_mem->LS, nls_mem->J, nls_mem->x, b, ZERO);
  N_VScale(ONE, nls_mem->x, b);

  return flag;
}

int ARKBraidTheta_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                              realtype tol, N_Vector ewt, void* mem)
{
  realtype delnrm;

  /* Compute the norm of the correction */
  delnrm = N_VWrmsNorm(del, ewt);

  if (delnrm <= tol) return SUN_NLS_SUCCESS;
  else return SUN_NLS_CONTINUE;
}

int ARKBraidTheta_NlsSetup(ARKodeMem ark_mem, ARKBraidNlsMem nls_mem,
                           SUNNonlinearSolver* NLS_ptr)
{
  int flag;               /* return flag */
  SUNNonlinearSolver NLS; /* output, nonlinear solver object */

  /* Check input */
  if (ark_mem == NULL || nls_mem == NULL) return SUNBRAID_ILLINPUT;

  NLS = SUNNonlinSol_Newton(nls_mem->th0, ark_mem->sunctx);
  if (NLS == NULL) return SUNBRAID_MEMFAIL;

  /* set residual/jacobian functions */
  flag = SUNNonlinSolSetSysFn(NLS, ARKBraidTheta_NlsResidual);
  if (flag != SUN_NLS_SUCCESS) return flag;

  flag = SUNNonlinSolSetLSetupFn(NLS, ARKBraidTheta_NlsLSetup);
  if (flag != SUN_NLS_SUCCESS) return flag;

  flag = SUNNonlinSolSetLSolveFn(NLS, ARKBraidTheta_NlsLSolve);
  if (flag != SUN_NLS_SUCCESS) return flag;

  flag = SUNNonlinSolSetConvTestFn(NLS, ARKBraidTheta_NlsConvTest, NULL);
  if (flag != SUN_NLS_SUCCESS) return flag;

  /* Set the max iterations */
  flag = SUNNonlinSolSetMaxIters(NLS, 100);
  if (flag != SUN_NLS_SUCCESS) return flag;

  /* Set the nonlinear solver */
  *NLS_ptr = NLS;

  return SUNBRAID_SUCCESS;
}

/* Solve the nonlinear order conditions given rhs */
int ARKBraidTheta_NlsSolve(SUNNonlinearSolver NLS, ARKBraidNlsMem nls_mem,
                           realtype* rhs)
{
  int flag; /* return flag */

  /* Check input */
  if (NLS == NULL || nls_mem == NULL || rhs == NULL) return SUNBRAID_ILLINPUT;

  /* Set right hand side */
  for (int i = 0; i < nls_mem->nconds; i++) NV_Ith_S(nls_mem->rhs, i) = rhs[i];

  /* Solve nonlinear system */
  flag = SUNNonlinSolSolve(NLS, nls_mem->th0, nls_mem->thcor, nls_mem->weight,
                           NLS_TOL, SUNTRUE, nls_mem);
  if (flag != SUN_NLS_SUCCESS) printf("NLS failed with flag %d\n", flag);
  return flag;

  /* Copy solution to output */
  N_VLinearSum(ONE, nls_mem->th0, ONE, nls_mem->thcor, nls_mem->thcur);

  return SUNBRAID_SUCCESS;
}

int ARKBraidTheta_NlsMem_Create(ARKBraidContent content, ARKBraidNlsMem* nlsmem)
{
  int flag;           /* return flag */
  int nc;             /* number of order conditions */

  ARKBraidNlsMem mem; /* output, nonlinear solver memory object */

  /* Check input */
  if (content == NULL) return SUNBRAID_ILLINPUT;
  if (content->ark_mem == NULL || content->ark_mem->sunctx == NULL)
    return SUNBRAID_MEMFAIL;

  /* Access content */
  SUNContext sunctx = (SUNContext)content->ark_mem->sunctx;
  nc                = content->num_order_conditions;

  /* Allocate memory */
  mem = NULL;
  mem = (ARKBraidNlsMem)malloc(sizeof(struct _ARKBraidNlsMem));
  if (mem == NULL) return SUNBRAID_ALLOCFAIL;

  mem->nconds = nc;
  mem->order  = content->order_coarse;
  mem->rhs    = NULL;
  mem->th0    = NULL;
  mem->thcur  = NULL;
  mem->thcor  = NULL;
  mem->weight = NULL;
  mem->x      = NULL;
  mem->J      = NULL;
  mem->LS     = NULL;

  /* Create vectors and matrices */
  mem->th0 = N_VNew_Serial(nc, sunctx);
  if (mem->th0 == NULL) return SUNBRAID_ALLOCFAIL;

  mem->rhs = N_VNew_Serial(nc, sunctx);
  if (mem->rhs == NULL) return SUNBRAID_ALLOCFAIL;

  mem->thcur = N_VNew_Serial(nc, sunctx);
  if (mem->thcur == NULL) return SUNBRAID_ALLOCFAIL;

  mem->thcor = N_VNew_Serial(nc, sunctx);
  if (mem->thcor == NULL) return SUNBRAID_ALLOCFAIL;

  mem->weight = N_VNew_Serial(nc, sunctx);
  if (mem->weight == NULL) return SUNBRAID_ALLOCFAIL;

  mem->x = N_VNew_Serial(nc, sunctx);
  if (mem->x == NULL) return SUNBRAID_ALLOCFAIL;

  mem->J = SUNDenseMatrix(nc, nc, sunctx);
  if (mem->J == NULL) return SUNBRAID_ALLOCFAIL;

  mem->LS = SUNLinSol_Dense(mem->th0, mem->J, sunctx);
  if (mem->LS == NULL) return SUNBRAID_ALLOCFAIL;

  flag = SUNLinSolInitialize(mem->LS);
  if (flag != SUNLS_SUCCESS) return SUNBRAID_SUNFAIL;

  /* Initialize workspace variables */
  N_VConst(ZERO, mem->rhs);
  N_VConst(ZERO, mem->thcor);

  /* Set weight for norm */
  N_VConst(ONE, mem->weight);

  /* Attach residual/jacobian functions and set initial guess */
  switch (content->order_coarse)
  {
  case 2:
    mem->res = _theta_sdirk2_res;
    mem->jac = _theta_sdirk2_jac;
    _theta_sdirk2_guess(NV_DATA_S(mem->th0));
    _theta_sdirk2_guess(NV_DATA_S(mem->thcur));
    break;
  case 3:
    mem->res = _theta_sdirk3_res;
    mem->jac = _theta_sdirk3_jac;
    _theta_sdirk3_guess(NV_DATA_S(mem->th0));
    _theta_sdirk3_guess(NV_DATA_S(mem->thcur));
    break;
  default: return SUNBRAID_ILLINPUT;
  }

  /* Create linear solver */
  mem->LS = SUNLinSol_Dense(mem->th0, mem->J, sunctx);
  if (mem->LS == NULL) return SUNBRAID_ALLOCFAIL;

  flag = SUNLinSolInitialize(mem->LS);
  if (flag != SUNLS_SUCCESS) return SUNBRAID_SUNFAIL;

  /* Return nonlinear solver memory */
  *nlsmem = mem;

  return SUNBRAID_SUCCESS;
}

int ARKBraidTheta_NlsMem_Free(ARKBraidNlsMem nlsmem)
{
  if (nlsmem == NULL) return SUNBRAID_SUCCESS;

  N_VDestroy_Serial(nlsmem->rhs);
  N_VDestroy_Serial(nlsmem->th0);
  N_VDestroy_Serial(nlsmem->thcur);
  N_VDestroy_Serial(nlsmem->thcor);
  N_VDestroy_Serial(nlsmem->weight);
  N_VDestroy_Serial(nlsmem->x);
  SUNMatDestroy(nlsmem->J);
  SUNLinSolFree(nlsmem->LS);

  free(nlsmem);
  nlsmem = NULL;

  return SUNBRAID_SUCCESS;
}