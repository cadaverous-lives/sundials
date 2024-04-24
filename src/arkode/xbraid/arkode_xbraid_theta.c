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
 * This is an implementation file for the ARKODE + XBraid interface.
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

int _ARKBraid_rootprintf(braid_StepStatus step_status, braid_SyncStatus sync_status, const char *format, ...)
{
  // printf(const char * __restrict, ...)
  braid_Int level;
  braid_Int iu, il, root;
  va_list   ap;

  va_start(ap, format);

  if (step_status)
  {
    braid_StepStatusGetLevel(step_status, &level);
    braid_StepStatusGetTIUL(step_status, &iu, &il, level);
  }
  else if (sync_status)
  {
    braid_SyncStatusGetLevel(sync_status, &level);
    braid_SyncStatusGetTIUL(sync_status, &iu, &il, level);
  }
  root = (il == 0);

  if (!root) return 0;
  return vfprintf(stdout, format, ap);
}

int _ARKBraidTheta_GetNumOrderConditions(int fine_order, int coarse_order)
{
  int num_order_conditions = 0;

  // if (coarse_order <= fine_order) return num_order_conditions;
  if (coarse_order >= 2) num_order_conditions += 1;
  if (coarse_order >= 3) num_order_conditions += 2;
  if (coarse_order >= 4) num_order_conditions += 4;
  if (coarse_order >= 5) num_order_conditions += 9;

  return num_order_conditions;
}

/* Not currently used... */
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
  vdata->tprior   = ZERO;
  vdata->etascale = ONE;
  vdata->Phi      = NULL;

  /* Allocate for order conditions (Phi) */
  int num_conditions = _ARKBraidTheta_GetNumOrderConditions(content->order_fine,
                                                            content->order_coarse);
  if (num_conditions > 0)
  {
    vdata->Phi = (realtype*)calloc(num_conditions, sizeof(realtype));
    if (vdata->Phi == NULL) return SUNBRAID_ALLOCFAIL;
    for (size_t i = 0; i < num_conditions; i++) {
      vdata->Phi[i] = ZERO;
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
  // printf("  Packing buffer: tprior=%f, etascale=%f, Phi0=%f\n", vdata->tprior, vdata->etascale, vdata->Phi[0]);
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

  vdata = (ARKBraidThetaVecData)(*vdata_ptr);

  /* Copy data from buffer */
  realtype* buf   = (realtype*)buffer;
  vdata->tprior   = buf[0];
  vdata->etascale = buf[1];
  for (int i = 0; i < content->num_order_conditions; i++) {
    vdata->Phi[i] = buf[i + 2];
  }
  // printf("unPacking buffer: tprior=%f, etascale=%f, Phi0=%f\n", vdata->tprior, vdata->etascale, vdata->Phi[0]);

  return SUNBRAID_SUCCESS;
}

/* this is called during Sync to allocate memory for coarse grid butcher tables */
int ARKBraidTheta_InitHierarchy(braid_App app, braid_SyncStatus sstatus)
{
  int flag;           /* return flag          */
  int caller;         /* XBraid calling function */
  int level, nlevels; /* XBraid level, number of levels */
  int iter;           /* XBraid iteration */
  int cf;             /* coarsening factor */
  ARKBraidContent content; /* ARKBraid app content */

  /* Check input */
  if (app == NULL || sstatus == NULL) return SUNBRAID_ILLINPUT;

  /* Access app content */
  content = (ARKBraidContent)app->content;

  /* Get information from XBraid status */
  // flag = braid_SyncStatusGetLevel(sstatus, &level);
  // CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_SyncStatusGetCallingFunction(sstatus, &caller);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_SyncStatusGetNLevels(sstatus, &nlevels);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  flag = braid_SyncStatusGetIter(sstatus, &iter);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);

  /* Setup nonlinear solver */
  if (content->theta_mem == NULL)
  {
    flag = ARKBraidTheta_NlsMem_Create(content, &content->theta_mem);
    if (flag != SUNBRAID_SUCCESS) return flag;
    flag = ARKBraidTheta_NlsSetup(content->ark_mem, content->theta_mem, &content->NLS);
    if (flag != SUNBRAID_SUCCESS) return flag;
  }

  /* Allocate memory for coarse grid butcher tables */
  flag = _ARKBraidTheta_AllocCGBtables(content);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Turn on computation of coarse grid Butcher tables */
  content->flag_refine_downcycle = SUNTRUE;
  _ARKBraid_rootprintf(NULL, sstatus,"  ARKBraid: Computing theta methods\n");

  /* make sure this is still computed even when we skip the first downcycle */
  if (content->flag_skip_downcycle)
  {
    /* Compute theta methods assuming uniform time-grid */ 

    ARKBraidGridData   grid;
    ARKodeButcherTable B   = content->fine_btable;
    ARKBraidThetaOrdCondsMem theta_mem = content->theta_mem;
    sunrealtype       *phi   = theta_mem->phi1;
    sunrealtype       *phi_m = theta_mem->phi2;
    long m = 1; /* Keeps a running total of coarsening factor relative to fine-grid */

    // printf("  ARKBraid: WARNING: Skipping first downcycle, theta methods will not be computed until next refinement\n");

    /* Elementary weights for fine-grid Butcher table */
    _phi_order2(B->b, B->c, B->stages, &phi[0]);
    phi[0] -= RCONST(1.)/RCONST(2.);

    if (content->order_coarse >= 3)
    {
      _phi_order3a(B->b, B->c, B->c, B->stages, &phi[1]);
      _phi_order3b(B->b, B->A, B->c, B->stages, &phi[2]);
      phi[1] -= RCONST(1.)/RCONST(3.);
      phi[2] -= RCONST(1.)/RCONST(6.);
    }

    if (content->order_coarse >= 4)
    {
      _phi_order4a(B->b, B->c, B->c, B->c, B->stages, &phi[3]);
      _phi_order4b(B->b, B->c, B->A, B->c, B->stages, &phi[4]);
      _phi_order4c(B->b, B->A, B->c, B->c, B->stages, &phi[5]);
      _phi_order4d(B->b, B->A, B->A, B->c, B->stages, &phi[6]);
      phi[3] -= RCONST(1.)/RCONST(4.);
      phi[4] -= RCONST(1.)/RCONST(8.);
      phi[5] -= RCONST(1.)/RCONST(12.);
      phi[6] -= RCONST(1.)/RCONST(24.);
    }

    for (level = 1; level < nlevels; level++)
    {
      /* Compute running product of coarsening factor */
      flag = braid_SyncStatusGetCFactor(sstatus, &cf, level-1);
      CHECK_BRAID_RETURN(content->last_flag_braid, flag);

      /* Stop computing when m overflows */
      /* (As m -> oo, the order conditions reduce to the classical ones) */
      if (m > LONG_MAX/cf) break;

      m *= cf;
      grid = content->grids[level];

      /* Second order */
      phi_m[0] = RCONST(.5) + phi[0]/m;

      if (content->order_coarse >= 3)
      {
        /* Third order*/
        phi_m[1] = RCONST(1.)/RCONST(3.) + phi[0]/m + (phi[1] - phi[0])/m/m;
        phi_m[2] = RCONST(1.)/RCONST(6.) + phi[0]/m + (phi[2] - phi[0])/m/m;
      }

      if (content->order_coarse >= 4)
      {
        /* Fourth order */
        phi_m[3] = RCONST(1.)/RCONST(4.) + phi[0]/m + (RCONST(3.)/RCONST(2.))*(phi[1] - phi[0])/m/m
                 + (phi[3] + phi[0]/RCONST(2.) - RCONST(3.)*phi[1]/RCONST(2.))/m/m/m;
        phi_m[4] = RCONST(1.)/RCONST(8.) + (RCONST(5.)/RCONST(6.))*phi[0]/m
                 + (phi[0]*phi[0]/RCONST(2.) - phi[0] + phi[1]/RCONST(2.) + phi[2]/RCONST(2.))/m/m
                 + (phi[4] - phi[0]*phi[0]/RCONST(2.) + phi[0]/RCONST(6.) - phi[1]/RCONST(2.) - phi[2]/RCONST(2.))/m/m/m;
        phi_m[5] = RCONST(1.)/RCONST(12.) + (RCONST(2.)/RCONST(3.))*phi[0]/m
                 + (phi[2] + phi[1]/RCONST(2.) - RCONST(3.)*phi[0]/RCONST(2.))/m/m
                 + (phi[5] - phi[2] - phi[1]/RCONST(2.) + RCONST(5.)*phi[0]/RCONST(6.))/m/m/m;
        phi_m[6] = RCONST(1.)/RCONST(24.) + (RCONST(1.)/RCONST(2.))*phi[0]/m
                 + (phi[2] + phi[0]*phi[0]/RCONST(2.) - phi[0])/m/m
                 + (phi[6] - phi[2] - phi[0]*phi[0]/RCONST(2.) + phi[0]/RCONST(2.))/m/m/m;
      }

      /* Solve order conditions */
      ARKBraidTheta_NlsSolve(content->NLS, theta_mem, phi_m);

      /* Set Butcher tables */
      for (int i = 0; i < grid->num_steps_stored; i++)
      {
        _ARKBraidTheta_SetBtable(grid->coarse_btables[i], content, 
                                 NV_DATA_S(theta_mem->thcur));
      }

    }
  }

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
    // _ARKBraid_rootprintf(status, NULL, "  ARKBraid: Turning off elementary weight computation\n");
    content->flag_refine_downcycle = SUNFALSE;
    content->flag_skip_downcycle   = SUNFALSE;
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
    realtype eta2, eta3, eta4;
    realtype etapr2, etapr3, etapr4;
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
        vdata->Phi[i] = ZERO;

      vdata->tprior   = tstart;
      vdata->etascale = cfactor * (tstop - tstart);
      // printf("Beginning of F-interval: ti=%d, tstart=%f, tstop=%f\n", ti, tstart, tstop);
    }

    /* etascale is a prediction of the total length of the coarse interval
     * we normalize by this for numerical stability (scalars are O(1/m) rather
     * than O(h))
     */
    eta    = (tstop - tstart) / vdata->etascale;
    eta_pr = (tstart - vdata->tprior) / vdata->etascale;

    /* Get workspace arrays */
    realtype *phi_step, *phi_pr;
    phi_step = content->theta_mem->phi1;
    phi_pr   = content->theta_mem->phi2;

    /* Compute elementary weights */
 
    /* Temporarily copy prior weights */
    for (int i = 0; i < content->num_order_conditions; i++)
      phi_pr[i] = vdata->Phi[i];

    /* These recurrences explained in (upcoming) theta methods paper... */

    /* Second order */
    eta2 = eta*eta;

    _phi_order2(B->b, B->c, B->stages, &phi_step[0]);
    vdata->Phi[0] += eta2 * phi_step[0] + eta * eta_pr;

    /* Third order */
    if (content->order_coarse >= 3)
    {
      eta3 = eta2*eta;
      etapr2 = eta_pr*eta_pr;

      _phi_order3a(B->b, B->c, B->c, B->stages, &phi_step[1]);
      _phi_order3b(B->b, B->A, B->c, B->stages, &phi_step[2]);
      vdata->Phi[1] += eta3 * phi_step[1] + RCONST(2.) * eta_pr * eta2 * phi_step[0] + etapr2 * eta;
      vdata->Phi[2] += eta * phi_pr[0] + eta3 * phi_step[2] + eta_pr * eta2 * phi_step[0];
    }

    /* Fourth order */
    if (content->order_coarse >= 4)
    {
      eta4 = eta3*eta;
      etapr3 = etapr2*eta_pr;

      _phi_order4a(B->b, B->c, B->c, B->c, B->stages, &phi_step[3]);
      _phi_order4b(B->b, B->c, B->A, B->c, B->stages, &phi_step[4]);
      _phi_order4c(B->b, B->A, B->c, B->c, B->stages, &phi_step[5]);
      _phi_order4d(B->b, B->A, B->A, B->c, B->stages, &phi_step[6]);
      vdata->Phi[3] += eta4 * phi_step[3] + RCONST(3.) * eta_pr * eta3 * phi_step[1]
                    +  RCONST(3.) * etapr2 * eta2 * phi_step[0] + etapr3 * eta;
      vdata->Phi[4] += eta2 * phi_step[0] * phi_pr[0] + eta4 * phi_step[4] + eta_pr * eta3 * phi_step[1]
                    +  eta_pr * eta * phi_pr[0] + eta_pr * eta3 * phi_step[2] + etapr2 * eta2 * phi_step[0];
      vdata->Phi[5] += eta * phi_pr[1] + eta4 * phi_step[5] + RCONST(2.) * eta_pr * eta3 * phi_step[2] + etapr2 * eta2 * phi_step[0];
      vdata->Phi[6] += eta * phi_pr[2] + eta2 * phi_step[0] * phi_pr[0] + eta4 * phi_step[6] + eta_pr * eta3 * phi_step[2];
    }
    
    /* TODO: Fifth/Sixth order... */

    /* Last step of f-interval */
    if (_ARKBraid_IsCPoint(ti + 1, cfactor))
    {
      int iu, il; /* highest and lowest time indices on this processor */
      int tic;    /* time index on next coarsest level */

      /* Computing on level for level + 1 */
      flag = braid_StepStatusGetTIUL(status, &iu, &il, level + 1);
      CHECK_BRAID_RETURN(content->last_flag_braid, flag);

      tic = (ti+1) / cfactor - il;

      /* Normalize weights */
      eta_c = (tstop - vdata->tprior) / vdata->etascale;

      /* Second order */
      vdata->Phi[0] = vdata->Phi[0]/eta_c/eta_c;

      /* Third order */
      if (content->order_coarse >= 3)
      {
        vdata->Phi[1] = vdata->Phi[1]/eta_c/eta_c/eta_c; 
        vdata->Phi[2] = vdata->Phi[2]/eta_c/eta_c/eta_c; 
      }

      /* Fourth order */
      if (content->order_coarse >= 4)
      {
        vdata->Phi[3] = vdata->Phi[3]/eta_c/eta_c/eta_c/eta_c;
        vdata->Phi[4] = vdata->Phi[4]/eta_c/eta_c/eta_c/eta_c;
        vdata->Phi[5] = vdata->Phi[5]/eta_c/eta_c/eta_c/eta_c;
        vdata->Phi[6] = vdata->Phi[6]/eta_c/eta_c/eta_c/eta_c;
      }

      /* Solve order conditions (solution in theta_mem->thcur)*/
      ARKBraidTheta_NlsSolve(content->NLS, content->theta_mem, vdata->Phi);

      /* Set new Butcher table */
      _ARKBraidTheta_SetBtable(content->grids[level + 1]->coarse_btables[tic], content,
                               NV_DATA_S(content->theta_mem->thcur));
    }
  }

  /* Set the Butcher table for this step TODO: support for explicit/ImEx methods */
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
  int ns;            /* Number of stages                                 */
  booleantype guess; /* If we don't know theta, provide an initial guess */
  realtype* A;       /* Butcher table entries                            */
  ARKBraidThetaOrdCondsMem theta_mem;

  /* Check inputs */
  if (B == NULL || content == NULL || content->theta_mem == NULL) return SUNBRAID_ILLINPUT;

  theta_mem = content->theta_mem;
  ns = theta_mem->nstages;
  if (B->stages != ns) return SUNBRAID_ILLINPUT;

  guess = (theta == NULL);

  A = (realtype*)calloc(ns * ns, sizeof(realtype));
  if (A == NULL) return SUNBRAID_ALLOCFAIL;

  if (guess)
  {
    theta = (realtype*)calloc(content->num_order_conditions, sizeof(realtype));
    if (theta == NULL)
    {
      free(A);
      return SUNBRAID_ALLOCFAIL;
    }
  }

  /* Set Butcher table */
  if (guess) theta_mem->init(theta);
  theta_mem->btable_A(A, theta);
  theta_mem->btable_b(B->b, theta);
  theta_mem->btable_c(B->c, theta);

  B->q = content->order_fine;
  B->p = content->order_fine;

  for (int i = 0; i < ns; i++) {
    for (int j = 0; j < ns; j++) { 
      B->A[i][j] = A[i * ns + j]; 
    }
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
  int tir;    /* time index, relative to first time point owned by current proc */
  ARKBraidGridData grid;

  if (content == NULL || content->grids == NULL) return SUNBRAID_MEMFAIL;
  grid = content->grids[level];

  /* Get the coarse (relative) time index */
  flag = braid_StepStatusGetTIUL(status, &iu, &il, level);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  tir = ti - il + 1;

  // printf("GetBTable: il=%d, iu=%d, level=%d, ti=%d\n", il, iu, level, ti);
  if (level == 0 || !content->use_theta ||
      grid->coarse_btables == NULL || grid->coarse_btables[tir] == NULL)
  {
    /* Use the fine table */
    *B = content->fine_btable;
  }
  else
  {
    /* Use the coarse table */
    *B = grid->coarse_btables[tir];
  }
  return SUNBRAID_SUCCESS;
}

int _ARKBraidTheta_AllocCGBtables(ARKBraidContent  content)
{
  int flag;     /* Return flag */
  int ntpoints; /* Number of time-points on this level */
  int ns;       /* Number of theta method stages */
  int nlevels;  /* Number of levels in MGRIT hierarchy */
  ARKBraidGridData grid; /* the butcher tables are stored here */

  /* Check input */
  if (content == NULL || content->theta_mem == NULL || content->grids == NULL) return SUNBRAID_ILLINPUT;

  nlevels = content->num_levels_alloc;
  ns = content->theta_mem->nstages;

  /* Preallocate Butcher tables for each coarse level */

  for (int lvl = 1; lvl < nlevels; lvl++)
  {
    grid = content->grids[lvl];
    ntpoints = grid->num_steps_stored;

    grid->coarse_btables = NULL;
    grid->coarse_btables = (ARKodeButcherTable*)calloc(ntpoints, sizeof(ARKodeButcherTable));
    if (grid->coarse_btables == NULL) return SUNBRAID_ALLOCFAIL;

    for (int ti = 0; ti < ntpoints; ti++)
    {
      grid->coarse_btables[ti] = NULL;
      grid->coarse_btables[ti] = ARKodeButcherTable_Alloc(ns, 0);
      if (grid->coarse_btables[ti] == NULL) return SUNBRAID_ALLOCFAIL;
      /* Set a sensible default (max classical order) */
      _ARKBraidTheta_SetBtable(grid->coarse_btables[ti], content, NULL);
    }
  }

  return SUNBRAID_SUCCESS;
}

/* ---------------------------------
 * Solve dense nonlinear systems
 * --------------------------------- */

#define NLS_TOL RCONST(1.0e-10) /* nonlinear solver tolerance */

int ARKBraidTheta_NlsResidual(N_Vector thcor, N_Vector r, void* mem)
{
  ARKBraidThetaOrdCondsMem nlsmem;

  /* Check inputs */
  if (mem == NULL) return SUNBRAID_ILLINPUT;
  nlsmem = (ARKBraidThetaOrdCondsMem)mem;

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
  ARKBraidThetaOrdCondsMem nls_mem;

  if (mem == NULL) return SUNBRAID_ILLINPUT;
  nls_mem = (ARKBraidThetaOrdCondsMem)mem;

  /* Compute Jacobian */
  nls_mem->jac(SUNDenseMatrix_Data(nls_mem->J), NV_DATA_S(nls_mem->thcur));

  /* Update Jacobian status */
  *jcur = SUNTRUE;

  /* Setup linear solver */
  flag = SUNLinSolSetup(nls_mem->LS, nls_mem->J);

  return SUNBRAID_SUCCESS;
}

int ARKBraidTheta_NlsLSolve(N_Vector b, void* mem)
{
  int flag; /* return flag */
  ARKBraidThetaOrdCondsMem nls_mem;

  if (mem == NULL) return SUNBRAID_ILLINPUT;
  nls_mem = (ARKBraidThetaOrdCondsMem)mem;

  /* Solve linear system */
  flag = SUNLinSolSolve(nls_mem->LS, nls_mem->J, nls_mem->x, b, ZERO);
  N_VScale(ONE, nls_mem->x, b);

  return flag;
}

int ARKBraidTheta_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                              realtype tol, N_Vector ewt, void* data)
{
  SUNNonlinearSolverContent_Newton nls_content = (SUNNonlinearSolverContent_Newton)NLS->content;
  realtype delnrm;

  /* Compute the norm of the correction */
  delnrm = N_VWrmsNorm(del, ewt);

  if (delnrm <= tol) return SUN_NLS_SUCCESS;
  else if (nls_content->nconvfails < 100)
  {
    /* These flags force the newton solver to recompute the Jacobian every step */
    nls_content->jcur = SUNFALSE;
    return SUN_NLS_CONV_RECVR;
  }
  else return SUN_NLS_CONTINUE;
}

int ARKBraidTheta_NlsSetup(ARKodeMem ark_mem, ARKBraidThetaOrdCondsMem nls_mem,
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
int ARKBraidTheta_NlsSolve(SUNNonlinearSolver NLS, ARKBraidThetaOrdCondsMem nls_mem,
                           realtype* rhs)
{
  int flag = -1; /* return flag */
  int resets = 0;   /* Number of times solve has been restarted with different initial iterate */
  int reset_flag = 0;

  /* Check input */
  if (NLS == NULL || nls_mem == NULL || rhs == NULL) return SUNBRAID_ILLINPUT;

  /* Set right hand side */
  for (int i = 0; i < nls_mem->nconds; i++) NV_Ith_S(nls_mem->rhs, i) = rhs[i];

  /* Solve nonlinear system */
  while (flag != SUN_NLS_SUCCESS)
  {
    flag = SUNNonlinSolSolve(NLS, nls_mem->th0, nls_mem->thcor, nls_mem->weight,
                            NLS_TOL, SUNTRUE, nls_mem);
    if (flag != SUN_NLS_SUCCESS)
    {
      /* Reset with fall-back initial guess */
      if (nls_mem->alt_init)
      {
        reset_flag = nls_mem->alt_init(NV_DATA_S(nls_mem->th0), resets);
        reset_flag = nls_mem->alt_init(NV_DATA_S(nls_mem->thcur), resets);
        if (reset_flag == SUNFALSE) {
          // printf("Failed...\n");
          break;
        }
        // printf("resetting %d...\n", resets);
      }
      else
      {
        /* Reset once in case the data from previous solves messed up this one */
        nls_mem->init(NV_DATA_S(nls_mem->th0));
        nls_mem->init(NV_DATA_S(nls_mem->thcur));
        if (resets > 0) {
          // printf("Failed...\n");
          break;
        }
        // printf("resetting %d...\n", resets);
      }
      N_VConst(0., nls_mem->thcor);
      resets++;
    } 
  }

  if (flag != SUN_NLS_SUCCESS)
  {
    printf("  ARKBraid: WARNING: theta NLS failed with flag %d\n", flag);
    printf("            rhs was: [");
    for (int i = 0; i < nls_mem->nconds; i++)
    {
      printf("%f", rhs[i]);
      if (i < nls_mem->nconds - 1)
        printf(", ");
    }
    printf("]\n");
    return flag;
  }

  /* Copy solution to output */
  N_VLinearSum(ONE, nls_mem->th0, ONE, nls_mem->thcor, nls_mem->thcur);

  return SUNBRAID_SUCCESS;
}

int ARKBraidTheta_NlsMem_Create(ARKBraidContent content, ARKBraidThetaOrdCondsMem* nlsmem)
{
  int flag;           /* return flag */
  int nc;             /* number of order conditions */

  ARKBraidThetaOrdCondsMem mem; /* output, nonlinear solver memory object */

  /* Check input */
  if (content == NULL) return SUNBRAID_ILLINPUT;
  if (content->ark_mem == NULL || content->ark_mem->sunctx == NULL)
    return SUNBRAID_MEMFAIL;

  /* Access content */
  SUNContext sunctx = (SUNContext)content->ark_mem->sunctx;
  nc                = content->num_order_conditions;

  /* Allocate memory */
  mem = NULL;
  mem = (ARKBraidThetaOrdCondsMem)malloc(sizeof(struct _ARKBraidThetaOrdCondsMem));
  if (mem == NULL) return SUNBRAID_ALLOCFAIL;

  mem->order   = content->order_coarse;
  mem->nconds  = nc;
  mem->nstages = 0;
  mem->phi1    = NULL;
  mem->phi2    = NULL;
  mem->rhs     = NULL;
  mem->th0     = NULL;
  mem->thcur   = NULL;
  mem->thcor   = NULL;
  mem->weight  = NULL;
  mem->x       = NULL;
  mem->J       = NULL;
  mem->LS      = NULL;

  /* Create vectors and matrices */
  mem->phi1 = calloc(nc, sizeof(realtype));
  if (mem->phi1 == NULL) return SUNBRAID_ALLOCFAIL;

  mem->phi2 = calloc(nc, sizeof(realtype));
  if (mem->phi2 == NULL) return SUNBRAID_ALLOCFAIL;

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

  /* Initialize workspace variables */
  N_VConst_Serial(ZERO, mem->rhs);
  N_VConst_Serial(ZERO, mem->thcor);

  /* Set weight for norm */
  N_VConst_Serial(ONE, mem->weight);

  /* Attach residual/jacobian functions and set initial guess */
  mem->alt_init = NULL;
  switch (content->order_coarse)
  {
  case 2:
    mem->nstages = _theta_sdirk2_btable_ns;
    mem->res = _theta_sdirk2_res;
    mem->jac = _theta_sdirk2_jac;
    mem->btable_A = _theta_sdirk2_btable_A;
    mem->btable_b = _theta_sdirk2_btable_b;
    mem->btable_c = _theta_sdirk2_btable_c;
    mem->init = _theta_sdirk2_guess;
    _theta_sdirk2_guess(NV_DATA_S(mem->th0));
    _theta_sdirk2_guess(NV_DATA_S(mem->thcur));
    break;
  case 3:
    mem->nstages = _theta_sdirk3_btable_ns;
    mem->res = _theta_sdirk3_res;
    mem->jac = _theta_sdirk3_jac;
    mem->btable_A = _theta_sdirk3_btable_A;
    mem->btable_b = _theta_sdirk3_btable_b;
    mem->btable_c = _theta_sdirk3_btable_c;
    mem->init = _theta_sdirk3_guess;
    mem->alt_init = _theta_sdirk3_altguess;
    _theta_sdirk3_guess(NV_DATA_S(mem->th0));
    _theta_sdirk3_guess(NV_DATA_S(mem->thcur));
    break;
  case 4:
    mem->nstages = _theta_sdirk4_btable_ns;
    mem->res = _theta_sdirk4_res;
    mem->jac = _theta_sdirk4_jac;
    mem->btable_A = _theta_sdirk4_btable_A;
    mem->btable_b = _theta_sdirk4_btable_b;
    mem->btable_c = _theta_sdirk4_btable_c;
    mem->init = _theta_sdirk4_guess;
    _theta_sdirk4_guess(NV_DATA_S(mem->th0));
    _theta_sdirk4_guess(NV_DATA_S(mem->thcur));
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

int ARKBraidTheta_NlsMem_Free(ARKBraidThetaOrdCondsMem nlsmem)
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
  free(nlsmem->phi1);
  free(nlsmem->phi2);

  free(nlsmem);
  nlsmem = NULL;

  return SUNBRAID_SUCCESS;
}