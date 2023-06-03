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
 * This is the implementation file for the SUNDIALS + XBraid interface.
 * -------------------------------------------------------------------------- */

#include "sundials/sundials_xbraid.h"
#include "sundials/sundials_math.h"

#define ONE  RCONST(1.0)
#define ZERO RCONST(0.0)


/* -------------------------
 * Create and free utilities
 * ------------------------- */


/* Create an empty SUNBraidApp instance */
int SUNBraidApp_NewEmpty(braid_App *app)
{
  SUNBraidOps ops;

  /* Create XBraid interface object */
  *app = NULL;
  *app = (braid_App) malloc(sizeof(struct _braid_App_struct));
  if (*app == NULL) return SUNBRAID_ALLOCFAIL;

  /* Create operations structure */
  ops = NULL;
  ops = (SUNBraidOps) malloc(sizeof(struct _SUNBraidOps));
  if (ops == NULL)
  {
    free(*app);
    *app = NULL;
    return SUNBRAID_ALLOCFAIL;
  }

  /* Initialize operations to NULL */
  ops->getvectmpl   = NULL;
  ops->initvecdata  = NULL;
  ops->freevecdata  = NULL;
  ops->clonevecdata = NULL;
  ops->sumvecdata   = NULL;
  ops->getbufsize   = NULL;
  ops->bufpack      = NULL;
  ops->bufunpack    = NULL;

  /* Attach operations and initialize content to NULL */
  (*app)->ops     = ops;
  (*app)->content = NULL;

  return SUNBRAID_SUCCESS;
}


/* Free and empty SUNBraidApp instance */
int SUNBraidApp_FreeEmpty(braid_App *app)
{
  if (*app == NULL) return SUNBRAID_SUCCESS;

  if ((*app)->ops) free((*app)->ops);
  (*app)->ops = NULL;

  free(*app);
  *app = NULL;

  return SUNBRAID_SUCCESS;
}


/* ----------------------
 * Generic app operations
 * ---------------------- */


/* Get a template vector from the integrator */
int SUNBraidApp_GetVecTmpl(braid_App app, N_Vector *y)
{
  if (app->ops->getvectmpl == NULL) return SUNBRAID_OPNULL;
  return app->ops->getvectmpl(app, y);
}

/* Get the size of the buffer to be allocated */
int SUNBraidApp_GetBufSize(braid_App app, braid_Int *size_ptr)
{
  if (app->ops->getbufsize == NULL) return SUNBRAID_OPNULL;
  return app->ops->getbufsize(app, size_ptr);
}

/* Pack the buffer */
int SUNBraidApp_BufPack(braid_App app, void* buffer, void *vdata_ptr)
{
  if (app->ops->bufpack == NULL) return SUNBRAID_OPNULL;
  return app->ops->bufpack(app, buffer, vdata_ptr);
}

/* Unpack the buffer */
int SUNBraidApp_BufUnpack(braid_App app, void* buffer, void **vdata_ptr)
{
  if (app->ops->bufunpack == NULL) return SUNBRAID_OPNULL;
  return app->ops->bufunpack(app, buffer, vdata_ptr);
}


/* Initialize the vector data */
int SUNBraidApp_InitVecData(braid_App app, void** vdata_ptr)
{
  if (app->ops->initvecdata == NULL) return SUNBRAID_OPNULL;
  return app->ops->initvecdata(app, vdata_ptr);
}

/* Free the vector data */
int SUNBraidApp_FreeVecData(braid_App app, void *vdata_ptr)
{
  if (app->ops->freevecdata == NULL) return SUNBRAID_OPNULL;
  return app->ops->freevecdata(app, vdata_ptr);
}

/* Copy vector data */
int SUNBraidApp_CloneVecData(braid_App app, void* data_ptr, void** data_clone_ptr)
{
  if (app->ops->clonevecdata == NULL) return SUNBRAID_OPNULL;
  return app->ops->clonevecdata(app, data_ptr, data_clone_ptr);
}

/* Sum vector data */
int SUNBraidApp_SumVecData(braid_App app, braid_Real a, void* data_x_ptr, braid_Real b, void* data_y_ptr)
{
  if (app->ops->sumvecdata == NULL) return SUNBRAID_OPNULL;
  return app->ops->sumvecdata(app, a, data_x_ptr, b, data_y_ptr);
}

/* -------------------------
 * SUNBraid Vector Functions
 * ------------------------- */


/* Create a new vector wrapper */
int SUNBraidVector_New(braid_App app, N_Vector y, SUNBraidVector *u)
{
  /* Check for valid N_Vector */
  if (y == NULL) return SUNBRAID_ILLINPUT;

  /* Create new vector wrapper */
  *u = NULL;
  *u = (SUNBraidVector) malloc(sizeof(struct _braid_Vector_struct));
  if (*u == NULL) return SUNBRAID_ALLOCFAIL;

  /* Attach N_Vector */
  (*u)->y     = y;
  (*u)->vdata = NULL;

  /* Initialize vector data */
  SUNBraidApp_InitVecData(app, &((*u)->vdata));

  return SUNBRAID_SUCCESS;
}


/* Attach auxiliary vector data */
int SUNBraidVector_SetVecData(void *vdata, SUNBraidVector *u)
{
  /* Check for valid wrapper */
  if (u == NULL) return SUNBRAID_ILLINPUT;

  /* Attach vector data */
  (*u)->vdata = vdata;

  return SUNBRAID_SUCCESS;
}


/* Get the wrapped NVector */
int SUNBraidVector_GetNVector(SUNBraidVector u, N_Vector *y)
{
  /* Check for valid wrapper */
  if (u == NULL) return SUNBRAID_ILLINPUT;
  if (u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Extract NVector */
  *y = u->y;

  return SUNBRAID_SUCCESS;
}

/* Get the wrapped vector data */
int SUNBraidVector_GetVecData(SUNBraidVector u, void **vdata_ptr)
{
  /* Check for valid wrapper */
  if (u == NULL) return SUNBRAID_ILLINPUT;

  /* Extract vector data */
  *vdata_ptr = u->vdata;

  return SUNBRAID_SUCCESS;
}

/* Init wrapper function for XBraid */
int SUNBraidVector_Init(braid_App app, realtype t, braid_Vector *u_ptr)
{
  int      flag;
  N_Vector vy;

  /* Check for valid app */
  if (app == NULL) return SUNBRAID_ILLINPUT;

  /* Get template vector */
  flag = SUNBraidApp_GetVecTmpl(app, &vy);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Create new vector wrapper */
  flag = SUNBraidVector_New(app, vy, u_ptr);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Initialize vector data */
  flag = SUNBraidApp_InitVecData(app, &((*u_ptr)->vdata));
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Initialize vector data */
  N_VConst(ZERO, vy);

  return SUNBRAID_SUCCESS;
}

/* Create clone of an existing vector */
int SUNBraidVector_Clone(braid_App app, braid_Vector u, braid_Vector *v_ptr)
{
  int      flag;
  N_Vector vy;

  /* Check for valid wrapper */
  if (u == NULL) return SUNBRAID_ILLINPUT;
  if (u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Clone input NVector */
  vy = N_VClone(u->y);
  if (vy == NULL) return SUNBRAID_ALLOCFAIL;

  /* Create new vector wrapper */
  flag = SUNBraidVector_New(app, vy, v_ptr);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Copy data from u to v */
  N_VScale(ONE, u->y, vy);

  /* Copy the vector data */
  flag = SUNBraidApp_CloneVecData(app, u->vdata, &((*v_ptr)->vdata));
  if (flag != SUNBRAID_SUCCESS) return flag;

  return SUNBRAID_SUCCESS;
}


/* Free vector */
int SUNBraidVector_Free(braid_App app, braid_Vector u)
{
  /* Check for valid input */
  if (u == NULL) return SUNBRAID_SUCCESS;

  /* Destroy N_Vector */
  if (u->y != NULL)
  {
    N_VDestroy(u->y);
    u->y = NULL;
  }

  /* Destroy vector data */
  if (u->vdata != NULL)
  {
    SUNBraidApp_FreeVecData(app, u->vdata);
    u->vdata = NULL;
  }

  /* Destroy SUNBraidVector wrapper */
  free(u);
  u = NULL;

  return SUNBRAID_SUCCESS;
}


/* Compute alpha x + beta y -> y */
int SUNBraidVector_Sum(braid_App app, braid_Real alpha, braid_Vector x,
                       braid_Real beta, braid_Vector y)
{
  /* Check for valid wrappers */
  if (x == NULL || y == NULL) return SUNBRAID_ILLINPUT;
  if (x->y == NULL || y->y == NULL) return SUNBRAID_MEMFAIL;

  /* Compute linear sum */
  N_VLinearSum(alpha, x->y, beta, y->y, y->y);

  /* Sum vector data */
  SUNBraidApp_SumVecData(app, alpha, x->vdata, beta, y->vdata);
  return SUNBRAID_SUCCESS;
}


/* Compute L2 norm */
int SUNBraidVector_SpatialNorm(braid_App app, braid_Vector u,
                               braid_Real *norm_ptr)
{
  /* Check for valid wrapper */
  if (u == NULL) return SUNBRAID_ILLINPUT;
  if (u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Compute L2 norm */
  *norm_ptr = SUNRsqrt(N_VDotProd(u->y, u->y));

  return SUNBRAID_SUCCESS;
}


/* Compute message buffer size */
int SUNBraidVector_BufSize(braid_App app, braid_Int *size_ptr,
                           braid_BufferStatus bstatus)
{
  int      flag;  /* return flag     */
  N_Vector ytmpl; /* template vector */

  /* Get template vector */
  flag = SUNBraidApp_GetVecTmpl(app, &ytmpl);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Get buffer size */
  flag = N_VBufSize(ytmpl, size_ptr);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Get integrator buffer size */
  braid_Int tmp;
  SUNBraidApp_GetBufSize(app, &tmp);

  /* Update buffer size */
  *size_ptr += tmp;

  return SUNBRAID_SUCCESS;
}


/* Pack message buffer */
int SUNBraidVector_BufPack(braid_App app, braid_Vector u, void *buffer,
                           braid_BufferStatus bstatus)
{
  int flag; /* return flag */

  /* Check for valid wrapper */
  if (u == NULL) return SUNBRAID_ILLINPUT;
  if (u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Fill buffer */
  flag = N_VBufPack(u->y, buffer);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Get buffer offset */
  braid_Int offset;
  N_VBufSize(u->y, &offset);

  /* Pack integrator buffer starting at offset */
  void* tmp_buffer = (char*) buffer + offset;
  SUNBraidApp_BufPack(app, tmp_buffer, u->vdata);

  return SUNBRAID_SUCCESS;
}


/* Unpack message buffer */
int SUNBraidVector_BufUnpack(braid_App app, void *buffer, braid_Vector *u_ptr,
                             braid_BufferStatus bstatus)
{
  int      flag;  /* return flag     */
  N_Vector ytmpl; /* template vector */
  N_Vector y;     /* new NVector     */

  /* Check for valid input */
  if (buffer == NULL) return SUNBRAID_ILLINPUT;

  /* Get template vector */
  flag = SUNBraidApp_GetVecTmpl(app, &ytmpl);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Create new NVector */
  y = N_VClone(ytmpl);
  if (y == NULL) return SUNBRAID_ALLOCFAIL;

  /* Create new XBraid vector */
  flag = SUNBraidVector_New(app, y, u_ptr);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Unpack buffer */
  flag = N_VBufUnpack(y, buffer);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Get buffer offset */
  braid_Int offset;
  N_VBufSize(y, &offset);

  /* Unpack integrator buffer starting at offset */
  void* tmp_buffer = (char*) buffer + offset;
  SUNBraidApp_BufUnpack(app, tmp_buffer, (*u_ptr)->vdata);

  return SUNBRAID_SUCCESS;
}
