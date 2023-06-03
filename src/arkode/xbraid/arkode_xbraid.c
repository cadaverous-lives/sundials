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
#include "sundials/sundials_math.h"
#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"
#include "arkode_xbraid_theta_methods.c"

#include "arkode/arkode.h"

#include "arkode_xbraid_impl.h"
#include "arkode_arkstep_impl.h"

#define ONE RCONST(1.0)


/* ------------------------
 * Private helper functions
 * ------------------------ */


static booleantype _ARKBraid_IsCPoint(int tindex, int cfactor)
{
  if (tindex % cfactor == 0) return SUNTRUE;
  return SUNFALSE;
}

static int _ARKBraid_GetNumOrderConditions(int fine_order, int coarse_order)
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


/* -------------------------------
 * Construct, initialize, and free
 * ------------------------------- */


/* Create XBraid app strucutre */
int ARKBraid_Create(void *arkode_mem, braid_App *app)
{
  int              flag;
  ARKBraidContent  content;
  ARKodeARKStepMem step_mem;

  /* Check input */
  if (arkode_mem == NULL) return SUNBRAID_ILLINPUT;

  /* Create XBraid interface object */
  flag = SUNBraidApp_NewEmpty(app);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Set operations */
  (*app)->ops->getvectmpl  = ARKBraid_GetVecTmpl;
  (*app)->ops->getbufsize  = ARKBraid_GetBufSize;
  (*app)->ops->bufpack     = ARKBraid_BufPack;
  (*app)->ops->bufunpack   = ARKBraid_BufUnpack;
  (*app)->ops->initvecdata = ARKBraid_InitVecData;
  (*app)->ops->freevecdata = ARKBraid_FreeVecData;

  /* Create ARKODE interface content */
  content = NULL;
  content = (ARKBraidContent) malloc(sizeof(*content));
  if (content == NULL)
  {
    (void) SUNBraidApp_FreeEmpty(app);
    return SUNBRAID_ALLOCFAIL;
  }

  /* Initialize content */

  /* Attach ARKODE memory */
  content->ark_mem = (ARKodeMem) arkode_mem;

  /* Interface functions */
  content->step   = ARKBraid_Step;
  content->init   = ARKBraid_Init;
  content->snorm  = SUNBraidVector_SpatialNorm;
  content->access = ARKBraid_Access;

  /* Saved return flags */
  content->last_flag_braid  = SUNBRAID_SUCCESS;
  content->last_flag_arkode = SUNBRAID_SUCCESS;

  /* Output time and solution (allocated in access if necessary) */
  content->tout = content->ark_mem->tn;
  content->yout = NULL;

  /* flags for theta method */
  content->flag_refine_downcycle = SUNFALSE;
  content->flag_skip_downcycle   = SUNTRUE;

  /* Fine and coarse grid method orders */

  if (content->ark_mem->step_mem == NULL) return SUNBRAID_ILLINPUT;
  step_mem = (ARKodeARKStepMem) content->ark_mem->step_mem;

  content->order_fine   = step_mem->q;
  // TODO: should this default to q + 1?
  content->order_coarse = step_mem->q;

  content->num_order_conditions = _ARKBraid_GetNumOrderConditions(
    content->order_fine, content->order_coarse
  );

  content->num_levels     = 0;
  content->num_tables     = NULL;
  content->coarse_btables = NULL;
  content->fine_btable    = NULL;

  /* Attach content */
  (*app)->content = content;

  return SUNBRAID_SUCCESS;
}

int ARKBraid_InitVecData(braid_App app, void **vdata_ptr)
{
  int flag;
  ARKBraidVecData vdata;

  /* Check input */
  if (app == NULL) return SUNBRAID_ILLINPUT;

  /* Access content */
  ARKBraidContent content = (ARKBraidContent) app->content;
  if (content == NULL) return SUNBRAID_ILLINPUT;

  /* Create vector data */
  vdata = NULL;
  vdata = (ARKBraidVecData)malloc(sizeof(struct _ARKBraidVecData));
  if (vdata == NULL) return SUNBRAID_ALLOCFAIL;

  /* Create vector data */
  vdata->tprior   = RCONST(0.);
  vdata->etascale = RCONST(1.);
  vdata->Phi      = NULL;

  /* Allocate for order conditions (Phi) */
  int num_conditions = _ARKBraid_GetNumOrderConditions(content->order_fine, content->order_coarse);
  if (num_conditions > 0)
  {
    vdata->Phi = (realtype *) malloc(num_conditions * sizeof(realtype));
    if (vdata->Phi == NULL) return SUNBRAID_ALLOCFAIL;
  }

  /* Attach vector data */
  *vdata_ptr = (void*)vdata;

  return SUNBRAID_SUCCESS;
}


/* Initialize XBraid, attach interface functions */
int ARKBraid_BraidInit(MPI_Comm comm_w, MPI_Comm comm_t, realtype tstart,
                       realtype tstop, sunindextype ntime, braid_App app,
                       braid_Core *core)
{
  braid_Int       braid_flag;
  ARKBraidContent content;

  /* Check inputs */
  if (comm_w == MPI_COMM_NULL || comm_t == MPI_COMM_NULL || ntime < 2 ||
      app == NULL)
    return SUNBRAID_ILLINPUT;

  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  /* Shortcut to content */
  content = (ARKBraidContent) app->content;

  /* Initialize XBraid */
  braid_flag = braid_Init(comm_w, comm_t, tstart, tstop, ntime, app,
                          content->step,
                          content->init,
                          SUNBraidVector_Clone,
                          SUNBraidVector_Free,
                          SUNBraidVector_Sum,
                          content->snorm,
                          content->access,
                          SUNBraidVector_BufSize,
                          SUNBraidVector_BufPack,
                          SUNBraidVector_BufUnpack,
                          core);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  /* Set sync function */
  braid_flag = braid_SetSync(*core, ARKBraid_Sync);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  return SUNBRAID_SUCCESS;
}

/* Set coarse grid Butcher table given theta parameters */
static int _ARKBraid_SetBtable(ARKodeButcherTable B, ARKBraidContent content, realtype *theta)
{
  /* TODO: get/store somewhere the actual number of stages */
  int         s;     /* Number of stages                                 */
  booleantype guess; /* If we don't know theta, provide an initial guess */
  realtype   *A;     /* Butcher table entries                            */

  /* Check inputs */
  if (B == NULL || content == NULL) return SUNBRAID_ILLINPUT;
  guess = (theta == NULL);

  s = B->stages;

  A = (realtype*) malloc(s*s*sizeof(realtype));
  if (A == NULL) return SUNBRAID_ALLOCFAIL;

  if (guess)
  {
    theta = (realtype*) malloc(content->num_order_conditions*sizeof(realtype));
    if (theta == NULL) { free(A); return SUNBRAID_ALLOCFAIL; }
  }

  if (content->order_coarse == 2)
  {
    if (guess) _theta_sdirk2_guess(theta);

    _theta_sdirk2_btable_A(A, theta);
    _theta_sdirk2_btable_b(B->b, theta);
    _theta_sdirk2_btable_c(B->c, theta);
    B->q = 1;
    B->p = 1;
  }
  else if (content->order_coarse == 3)
  {
    if (guess) _theta_sdirk3_guess(theta);

    _theta_sdirk3_btable_A(A, theta);
    _theta_sdirk3_btable_b(B->b, theta);
    _theta_sdirk3_btable_c(B->c, theta);
    B->q = 2;
    B->p = 2;
  }
  else
  {
    free(A);
    if (guess) free(theta);
    return SUNBRAID_ILLINPUT;
  }

  for (int i=0; i<s; i++)
  {
    for (int j=0; j<s; j++)
    {
      B->A[i][j] = A[i*s+j];
    }
  } 

  free(A);
  if (guess) free(theta);
  return SUNBRAID_SUCCESS;
}

/* returns the Butcher table to be used for index ti on the given level */
static int _ARKBraid_GetBTable(ARKBraidContent content, braid_StepStatus status, braid_Int level, braid_Int ti, ARKodeButcherTable* B)
{
  int flag;    /* return flag */
  int iu, il;  /* upper and lower indices for the current level */
  int tir;     /* time index, relative to first time point owned by current proc */
  
  /* Get the coarse (relative) time index */
  flag = braid_StepStatusGetTIUL(status, &iu, &il, level);
  CHECK_BRAID_RETURN(content->last_flag_braid, flag);
  tir = ti - il;

  if (level == 0 || content->order_coarse <= content->order_fine 
                 || content->coarse_btables == NULL 
                 || content->coarse_btables[level] == NULL 
                 || content->coarse_btables[level][tir] == NULL) 
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

static int _ARKBraid_AllocCGBtables(ARKBraidContent content, braid_SyncStatus sstatus)
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

  content->coarse_btables = (ARKodeButcherTable **) malloc((nlevels-1) * sizeof(ARKodeButcherTable*));
  if (content->coarse_btables == NULL) return SUNBRAID_ALLOCFAIL;
  content->coarse_btables -= 1; /* offset array by one, since we are storing nothing for lvl 0 */

  content->num_tables = (int *) malloc((nlevels-1) * sizeof(int));
  if (content->num_tables == NULL) return SUNBRAID_ALLOCFAIL;
  content->num_tables -= 1; /* offset array by one, since we are storing nothing for lvl 0 */

  /* TODO: This is obviously not going to generalize past order 3 */
  ns = content->order_coarse;

  /* Preallocate Butcher tables for each coarse level */

  for (int lvl = 1; lvl < nlevels; lvl++)
  {
    content->coarse_btables[lvl] = NULL;

    flag = braid_SyncStatusGetTIUL(sstatus, &iu, &il, lvl);
    CHECK_BRAID_RETURN(content->last_flag_braid, flag);

    ntpoints = iu - il + 1;
    content->num_tables[lvl] = ntpoints;
    content->coarse_btables[lvl] = (ARKodeButcherTable *) malloc(ntpoints * sizeof(ARKodeButcherTable));
    if (content->coarse_btables[lvl] == NULL) return SUNBRAID_ALLOCFAIL;

    for (int ti = 0; ti < ntpoints; ti++)
    {
      content->coarse_btables[lvl][ti] = NULL;
      content->coarse_btables[lvl][ti] = ARKodeButcherTable_Alloc(ns, 0);
      if (content->coarse_btables[lvl][ti] == NULL) return SUNBRAID_ALLOCFAIL;
      /* Set a sensible default */
      _ARKBraid_SetBtable(content->coarse_btables[lvl][ti], content, NULL);
    }
  }

  return SUNBRAID_SUCCESS;
}

static int _ARKBraid_FreeCGBtables(ARKBraidContent content)
{
  if (content->coarse_btables != NULL)
  {
    for (int level = 0; level < content->num_levels; level++)
    {
      for (int table = 0; table < content->num_tables[level]; table++)
      {
        if (content->coarse_btables[level][table] != NULL)
          ARKodeButcherTable_Free(content->coarse_btables[level][table]);
      }
      free(content->coarse_btables[level]);
      content->coarse_btables[level] = NULL;
    }
    free(content->num_tables);
    free(content->coarse_btables);
    content->num_tables     = NULL;
    content->coarse_btables = NULL;
  }

  return SUNBRAID_SUCCESS;
}

/* Deallocate XBraid app structure */
int ARKBraid_Free(braid_App *app)
{
  ARKBraidContent content;  /* ARKBraid app content  */

  if (*app == NULL) return SUNBRAID_SUCCESS;

  if ((*app)->content != NULL)
  {
    content = (ARKBraidContent) (*app)->content;

    if (content->yout != NULL)
    {
      arkFreeVec(content->ark_mem, &(content->yout));
      content->yout = NULL;
    }

    if (content->fine_btable != NULL)
    {
      ARKodeButcherTable_Free(content->fine_btable);
      content->fine_btable = NULL;
    }
    _ARKBraid_FreeCGBtables(content);

    free((*app)->content);
    (*app)->content = NULL;
  }
  return SUNBraidApp_FreeEmpty(app);
}

int ARKBraid_FreeVecData(braid_App app, void *vdata_ptr)
{
  ARKBraidVecData vdata;

  if (vdata_ptr == NULL) return SUNBRAID_SUCCESS;

  vdata = (ARKBraidVecData)vdata_ptr;

  if (vdata->Phi != NULL)
  {
    free(vdata->Phi);
    vdata->Phi = NULL;
  }
  free(vdata);
  vdata_ptr = NULL;

  return SUNBRAID_SUCCESS;
}

/* ----------------------
 * ARKBraid Set Functions
 * ---------------------- */

int ARKBraid_SetCoarseOrder(braid_App app, int order)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent) app->content;

  /* Check order */
  if (order < 0 || order > 3) return SUNBRAID_ILLINPUT;

  /* Restore default or set default */
  if (order == 0)
    content->order_coarse = content->order_fine;
  else
    content->order_coarse = order;

  /* Get new number of order conditions */
  content->num_order_conditions = _ARKBraid_GetNumOrderConditions(
    content->order_fine, content->order_coarse
  );

  return SUNBRAID_SUCCESS;
}

int ARKBraid_SetStepFn(braid_App app, braid_PtFcnStep step)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent) app->content;

  /* Restore default or set function pointer */
  if (step == NULL)
    content->step = ARKBraid_Step;
  else
    content->step = step;

  return SUNBRAID_SUCCESS;
}


int ARKBraid_SetInitFn(braid_App app, braid_PtFcnInit init)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent) app->content;

  /* Restore default or set function pointer */
  if (init == NULL)
    content->init = ARKBraid_Init;
  else
    content->init = init;

  return SUNBRAID_SUCCESS;
}


int ARKBraid_SetSpatialNormFn(braid_App app, braid_PtFcnSpatialNorm snorm)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent) app->content;

  /* Restore default or set function pointer */
  if (snorm == NULL)
    content->snorm = SUNBraidVector_SpatialNorm;
  else
    content->snorm = snorm;

  return SUNBRAID_SUCCESS;
}


int ARKBraid_SetAccessFn(braid_App app, braid_PtFcnAccess access)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent) app->content;

  /* Restore default or set function pointer */
  if (access == NULL)
    content->access = ARKBraid_Access;
  else
    content->access = access;

  return SUNBRAID_SUCCESS;
}


/* ----------------------
 * ARKBraid Get Functions
 * ---------------------- */


int ARKBraid_GetVecTmpl(braid_App app, N_Vector *tmpl)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;
  *tmpl = content->ark_mem->yn;
  return SUNBRAID_SUCCESS;
}


int ARKBraid_GetARKStepMem(braid_App app, void **arkode_mem)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;
  *arkode_mem = (void*) content->ark_mem;
  return SUNBRAID_SUCCESS;
}


int ARKBraid_GetUserData(braid_App app, void **user_data)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;
  *user_data = content->ark_mem->user_data;
  return SUNBRAID_SUCCESS;
}


int ARKBraid_GetLastBraidFlag(braid_App app, int *last_flag)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
  *last_flag = content->last_flag_braid;
  return SUNBRAID_SUCCESS;
}


int ARKBraid_GetLastARKStepFlag(braid_App app, int *last_flag)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
  *last_flag = content->last_flag_arkode;
  return SUNBRAID_SUCCESS;
}


int ARKBraid_GetSolution(braid_App app, realtype *tout, N_Vector yout)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
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
  braid_Int       braid_flag;  /* braid function return flag       */
  int             ark_flag;    /* arkode step return flag          */
  int             flag;        /* arkode function return flag      */
  int             level;       /* current level                    */
  int             rfac;        /* refinement factor                */
  int             cfactor;     /* coarsening factor                */
  int             ti, tic;     /* time indices for fine and coarse */
  int             iu, il;      /* upper and lower indices this proc owns */
  int             caller;      /* caller of ARKBraid_Step          */
  int             fixedstep;   /* flag for fixed step size         */
  realtype        tstart;      /* current time                     */
  realtype        tstop;       /* evolve to this time              */
  realtype        eta, eta_pr; /* normalized time step sizes       */
  realtype        eta_c;       /* predicted normalized coarse step */
  realtype        hacc;        /* accuracy based step size         */
  ARKBraidContent  content;     /* ARKBraid app content             */

  ARKodeButcherTable  B = NULL;  /* Butcher table for the step  */
  ARKBraidVecData vdata = NULL;  /* vector data for storing order conditions */

  /* Check input */
  if (app == NULL || status == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL || u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Access app content */
  content = (ARKBraidContent) app->content;

  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;

  /* Remember if we are using fixed time-stepping */
  fixedstep = (content->ark_mem->fixedstep);

  /* Get info from XBraid */
  braid_flag = braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  braid_flag = braid_StepStatusGetLevel(status, &level);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  braid_flag = braid_StepStatusGetCallingFunction(status, &caller);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);


  braid_flag = braid_StepStatusGetTIndex(status, &ti);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  /* Get Butcher table for this step */
  _ARKBraid_GetBTable(content, status, level, ti, &B);
  if (B == NULL) return SUNBRAID_MEMFAIL;
  
  /* Compute theta method order conditions */

  /* F-relax */
  booleantype frelax = (caller == braid_ASCaller_FRestrict 
                     || caller == braid_ASCaller_FASResidual);
  if (content->flag_refine_downcycle && frelax && u->vdata != NULL)
  {
    braid_flag = braid_StepStatusGetCFactor(status, &cfactor, level);
    CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);
    
    /* Access vector data */
    vdata = (ARKBraidVecData) u->vdata;

    /* First step of f-interval */
    if (_ARKBraid_IsCPoint(ti, cfactor))
    {
      /* Reset order conditions */
      for (int i=0; i < content->num_order_conditions; i++)
        vdata->Phi[i] = RCONST(0.0);

      vdata->tprior   = tstart;
      vdata->etascale = cfactor * (tstop - tstart);
    }

    /* etascale is a prediction of the total length of the coarse interval 
     * we normalize by this for numerical stability (scalars are O(1/m) rather than O(h)) 
     */
    eta = (tstop - tstart) / vdata->etascale;
    eta_pr = (tstart - vdata->tprior) / vdata->etascale;
    realtype *phi_step, *phi_pr;

    /* Allocate temporary arrays TODO: preallocate in app->content? */
    phi_step = (realtype*) malloc(content->num_order_conditions * sizeof(realtype));
    if (phi_step == NULL) return SUNBRAID_MEMFAIL;
    phi_pr = (realtype*) malloc(content->num_order_conditions * sizeof(realtype));
    if (phi_pr == NULL) return SUNBRAID_MEMFAIL;

    /* Compute elementary weights */

    /* Temporarily copy prior weights */
    for (int i=0; i < content->num_order_conditions; i++)
      phi_pr[i] = vdata->Phi[i];

    /* Second order */
    _phi_order2(B->b, B->c, B->stages, &phi_step[0]);
    vdata->Phi[0] += eta * eta * phi_step[0] + eta * eta_pr;

    /* Third order */
    if (content->order_coarse >= 3)
    {
      _phi_order3a(B->b, B->c, B->c, B->stages, &phi_step[1]);
      _phi_order3b(B->b, B->A, B->c, B->stages, &phi_step[2]);
      vdata->Phi[1] += SUNRpowerI(eta, 3) * phi_step[1] 
                     + RCONST(2.) * eta * eta * eta_pr * phi_step[0] + eta * eta_pr * eta_pr;
      vdata->Phi[2] += SUNRpowerI(eta, 3) * phi_step[2] 
                     + eta * phi_pr[0] + eta * eta * eta_pr * phi_step[0];
    }

    /* Last step of f-interval */
    if (_ARKBraid_IsCPoint(ti+1, cfactor))
    {
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

      /* Here is where we would solve order conditions */
      
      /* For now just do nothing */
      /* Set the Butcher table */

      /* Computing on level for level + 1 */
      // braid_flag = braid_StepStatusGetTIUL(status, &iu, &il, level+1);
      // CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);
      // tic = ti / cfactor - il;
      // _ARKBraid_SetBtable(B, content, theta);
      // content->coarse_btables[level+1][tic] = B;
    }

    free(phi_step);
    free(phi_pr);
  }

  /* Turn off computation of elementary weights once upcycle starts */
  if (content->flag_refine_downcycle && caller == braid_ASCaller_FInterp)
  {
    printf("Turning off elementary weight computation\n");
    content->flag_refine_downcycle = SUNFALSE;
  }

  /* Turn off error estimation on coarse grids */
  // if (level > 0)
  //   arkSetFixedStep(content->ark_mem, tstop - tstart);

  /* Set the Butcher table TODO: support for explicit methods */
  // if (app->explicit && level == 0)
  //   flag = ARKStepSetTables(ark_mem, B->q, B->p, NULL, B);
  // else
  flag = ARKStepSetTables(content->ark_mem, B->q, B->p, B, NULL);
  CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

  /* Finally propagate the solution */
  flag = ARKBraid_TakeStep((void*)(content->ark_mem), tstart, tstop, u->y,
                           &ark_flag);
  CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

  /* Restore fixedstep value */
  content->ark_mem->fixedstep = fixedstep;

  /* Refine grid (XBraid will ignore if refinement is disabled) */

  /* Compute refinement factor */
  if (level == 0)
  {
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
    braid_flag = braid_StepStatusSetRFactor(status, rfac);
    CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);
  }

  return SUNBRAID_SUCCESS;
}


/* Create and initialize vectors */
int ARKBraid_Init(braid_App app, realtype t, braid_Vector *u_ptr)
{
  int             flag;        /* return flag          */
  N_Vector        y;           /* output N_Vector      */
  ARKBraidVecData vdata;       /* output ARKBraid vector data (used for Theta method) */
  ARKBraidContent content;     /* ARKBraid app content */

  /* Check input */
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  /* Access app content */
  content = (ARKBraidContent) app->content;

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
int ARKBraid_Access(braid_App app, braid_Vector u,
                    braid_AccessStatus astatus)
{
  braid_Int       braid_flag;  /* braid return flag    */
  braid_Int       done;        /* braid finished flag  */
  braid_Int       ntpoints;    /* num pts on fine grid */
  braid_Int       idx;         /* time index for u     */
  braid_Real      time;        /* time value for u     */
  ARKBraidContent content;     /* ARKBraid app content */

  /* Check input */
  if (app == NULL || u == NULL || astatus == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL || u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Access app content */
  content = (ARKBraidContent) app->content;

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
        if (!arkAllocVec(content->ark_mem, content->ark_mem->yn,
                         &(content->yout)))
          return SUNBRAID_ALLOCFAIL;
      }

      /* Save solution for output to user */
      content->tout = time;
      N_VScale(ONE, u->y, content->yout);
    }
  }

  return SUNBRAID_SUCCESS;
}

/* Sync is used by theta method to allocate memory for coarse grid butcher tables */
int ARKBraid_Sync(braid_App app, braid_SyncStatus sstatus)
{
  int             flag;            /* return flag          */
  ARKBraidContent content;         /* ARKBraid app content */
  int             caller;          /* XBraid calling function */
  int             level, nlevels;  /* XBraid level, number of levels */
  int             iter;            /* XBraid iteration */
  int             iu, il;          /* XBraid index of lowest and highest time indices owned by this processor */
  int             ncpoints;        /* XBraid number of coarse grid points */
  ARKodeARKStepMem step_mem;       /* ARKStep memory */

  /* Check input */
  if (app == NULL || sstatus == NULL) return SUNBRAID_ILLINPUT;

  /* Access app content */
  content = (ARKBraidContent) app->content;

  /* Need to store the fine-grid butcher table */
  if (content->fine_btable == NULL)
  {
    /* Initialize arkode_mem so that there is a valid Butcher table */   
    flag = arkStep_SetButcherTables(content->ark_mem);
    CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

    /* Get the ARKStep memory */
    step_mem = (ARKodeARKStepMem) content->ark_mem->step_mem;

    /* Get the fine-grid butcher table */
    if (step_mem->Be != NULL) 
      content->fine_btable = ARKodeButcherTable_Copy(step_mem->Be);
    else if (step_mem->Bi != NULL)
      content->fine_btable = ARKodeButcherTable_Copy(step_mem->Bi);
    else
      return SUNBRAID_ILLINPUT;
  }


  /* Check that theta method needs computing */
  if (content->order_fine >= content->order_coarse) return SUNBRAID_SUCCESS;

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
  if (caller == braid_ASCaller_Drive_AfterInit || caller == braid_ASCaller_FRefine_AfterInitHier)
  {
    _ARKBraid_FreeCGBtables(content);

    /* Allocate memory for coarse grid butcher tables */
    flag = _ARKBraid_AllocCGBtables(content, sstatus);
    if (flag != SUNBRAID_SUCCESS) return flag;

    /* Turn on computation of coarse grid Butcher tables */
    if (nlevels > 1)
    {
      content->flag_refine_downcycle = SUNTRUE;
      printf("ARKBraid: Computing coarse grid Butcher tables\n");
    }
  }

  /* make sure this is still computed even when we skip the first downcycle */
  if (content->flag_skip_downcycle && caller == braid_ASCaller_Drive_TopCycle && iter == 0)
  {
    content->flag_refine_downcycle = SUNTRUE;
    printf("ARKBraid: Skipping first downcycle, but still computing coarse grid Butcher tables\n");
  }

  return SUNBRAID_SUCCESS;
}

/* return the buffer space needed, this will get added to the buffer space needed for the vector data in SUNBraidVector_BufSize */
int ARKBraid_GetBufSize(braid_App app, braid_Int *size_ptr)
{
  ARKBraidContent content; /* ARKBraid app content */

  /* Check input */
  if (app == NULL || size_ptr == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  /* Access app content */
  content = (ARKBraidContent) app->content;

  /* Set buffer size */
  *size_ptr = (content->num_order_conditions + 2) * sizeof(realtype);

  return SUNBRAID_SUCCESS;
}
 
/* pack/unpack data in/out of the buffer, SUNBraidVector_BufPack will offset the input pointer past the vector data */
int ARKBraid_BufPack(braid_App app, void *buffer, void *vdata_ptr)
{
  /* Check input */
  if (app == NULL || buffer == NULL || vdata_ptr == NULL) return SUNBRAID_ILLINPUT;

  /* Access vector data */
  ARKBraidVecData vdata = (ARKBraidVecData) vdata_ptr;

  /* Access app content */
  ARKBraidContent content = (ARKBraidContent) app->content;

  /* Copy data into buffer */
  realtype *buf = (realtype *) buffer;
  buf[0] = vdata->tprior;
  buf[1] = vdata->etascale;
  for (int i = 0; i < content->num_order_conditions; i++)
    buf[i+2] = vdata->Phi[i];

  return SUNBRAID_SUCCESS;
}

int ARKBraid_BufUnpack(braid_App app, void *buffer, void **vdata_ptr)
{
  int             flag;  /* return flag                  */
  ARKBraidVecData vdata; /* ARKBraid vector data, output */

  /* Check input */
  if (app == NULL || buffer == NULL || vdata_ptr == NULL) return SUNBRAID_ILLINPUT;

  /* Access app content */
  ARKBraidContent content = (ARKBraidContent) app->content;
  
  /* Allocate vdata */
  flag = ARKBraid_InitVecData(content, vdata_ptr);
  if (flag != SUNBRAID_SUCCESS) return flag;

  vdata = (ARKBraidVecData) *vdata_ptr;

  /* Copy data from buffer */
  realtype *buf = (realtype *) buffer;
  vdata->tprior = buf[0];
  vdata->etascale = buf[1];
  for (int i = 2; i < content->num_order_conditions; i++)
    vdata->Phi[i] = buf[i];

  /* Return vector data */
  *vdata_ptr = vdata;

  return SUNBRAID_SUCCESS;
}


/* -----------------
 * Utility Functions
 * ----------------- */


/* Force a single step with ARKEvolve */
int ARKBraid_TakeStep(void *arkode_mem, realtype tstart, realtype tstop,
                      N_Vector y, int *ark_flag)
{
  int      flag;      /* generic return flag      */
  int      tmp_flag;  /* evolve return flag       */
  realtype tret;      /* return time              */

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

  /* Take step, check flag below */
  tmp_flag = ARKStepEvolve(arkode_mem, tstop, y, &tret, ARK_ONE_STEP);

  /* Re-enable temporal error test check */
  flag = arkSetForcePass(arkode_mem, SUNFALSE);
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



/* ------------------------
 * Solve nonlinear systems
 * ------------------------ */


int ARKBraidNlsMem_Create(ARKBraidContent content, ARKBraidNlsMem *nlsmem)
{
  int flag; /* return flag */
  int nc;   /* number of order conditions */

  ARKBraidNlsMem mem;  /* output, nonlinear solver memory object */

  /* Check input */ 
  if (content == NULL) return SUNBRAID_ILLINPUT; 
  if (content->ark_mem == NULL || content->ark_mem->sunctx == NULL)
    return SUNBRAID_MEMFAIL;

  /* Access content */
  ARKodeMem ark_mem = (ARKodeMem) content->ark_mem;
  SUNContext sunctx = (SUNContext) ark_mem->sunctx;
  nc = content->num_order_conditions;

  /* Allocate memory */
  mem = NULL;
  mem = (ARKBraidNlsMem) malloc(sizeof(struct _ARKBraidNlsMem));
  if (mem == NULL) return SUNBRAID_ALLOCFAIL;

  mem->y0 = N_VNew_Serial(nc, sunctx);
  if (mem->y0 == NULL) return SUNBRAID_ALLOCFAIL;

  mem->ycur = N_VClone_Serial(mem->y0);
  if (mem->ycur == NULL) return SUNBRAID_ALLOCFAIL;

  mem->ycor = N_VClone_Serial(mem->y0);
  if (mem->ycor == NULL) return SUNBRAID_ALLOCFAIL;

  mem->w = N_VClone_Serial(mem->y0);
  if (mem->w == NULL) return SUNBRAID_ALLOCFAIL;

  mem->x = N_VClone_Serial(mem->y0);
  if (mem->x == NULL) return SUNBRAID_ALLOCFAIL;

  mem->A = SUNDenseMatrix(nc, nc, sunctx);
  if (mem->A == NULL) return SUNBRAID_ALLOCFAIL;

  mem->LS = SUNLinSol_Dense(mem->y0, mem->A, sunctx);
  if (mem->LS == NULL) return SUNBRAID_ALLOCFAIL;

  flag = SUNLinSolInitialize(mem->LS);
  if (flag != SUNLS_SUCCESS) return SUNBRAID_SUNFAIL;

  /* Return nonlinear solver memory */
  *nlsmem = mem;

  return SUNBRAID_SUCCESS;
}

int ARKBraidNlsMem_Free(ARKBraidNlsMem nlsmem)
{
  if (nlsmem == NULL) return SUNBRAID_SUCCESS;

  N_VDestroy_Serial(nlsmem->y0);
  N_VDestroy_Serial(nlsmem->ycur);
  N_VDestroy_Serial(nlsmem->ycor);
  N_VDestroy_Serial(nlsmem->w);
  N_VDestroy_Serial(nlsmem->x);
  SUNMatDestroy(nlsmem->A);
  SUNLinSolFree(nlsmem->LS);

  free(nlsmem);
  nlsmem = NULL;

  return SUNBRAID_SUCCESS;
}
