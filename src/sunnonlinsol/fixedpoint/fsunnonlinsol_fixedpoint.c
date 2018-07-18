/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *-----------------------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *-----------------------------------------------------------------------------
 * This file contains the implementation of functions needed for initialization
 * of the SUNNonlinearSolver fixed point module operations in Fortran.
 *---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "fsunnonlinsol_fixedpoint.h"

/* Define global nonlinsol variables */

SUNNonlinearSolver F2C_CVODE_nonlinsol;
SUNNonlinearSolver F2C_IDA_nonlinsol;
SUNNonlinearSolver F2C_ARKODE_nonlinsol;

/* Declarations of external global variables */

extern N_Vector F2C_CVODE_vec;
extern N_Vector F2C_IDA_vec;
extern N_Vector F2C_ARKODE_vec;

/* Fortran callable interfaces */

void FSUNFIXEDPOINT_INIT(int *code, int *m, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    if (F2C_CVODE_nonlinsol)  SUNNonlinSolFree(F2C_CVODE_nonlinsol);
    F2C_CVODE_nonlinsol = NULL;
    F2C_CVODE_nonlinsol = SUNNonlinSol_FixedPoint(F2C_CVODE_vec, *ier);
    if (F2C_CVODE_nonlinsol == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    if (F2C_IDA_nonlinsol)  SUNNonlinSolFree(F2C_IDA_nonlinsol);
    F2C_IDA_nonlinsol = NULL;
    F2C_IDA_nonlinsol = SUNNonlinSol_FixedPoint(F2C_IDA_vec, *ier);
    if (F2C_IDA_nonlinsol == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    if (F2C_ARKODE_nonlinsol)  SUNNonlinSolFree(F2C_ARKODE_nonlinsol);
    F2C_ARKODE_nonlinsol = NULL;
    F2C_ARKODE_nonlinsol = SUNNonlinSol_FixedPoint(F2C_ARKODE_vec, *ier);
    if (F2C_ARKODE_nonlinsol == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


void FSUNFIXEDPOINT_SETMAXITERS(int *code, int *maxiters, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    if (!F2C_CVODE_nonlinsol) {
      *ier = -1;
      return;
    }
    *ier = SUNNonlinSolSetMaxIters(F2C_CVODE_nonlinsol, *maxiters);
    break;
  case FCMIX_IDA:
    if (!F2C_IDA_nonlinsol) {
      *ier = -1;
      return;
    }
    *ier = SUNNonlinSolSetMaxIters(F2C_IDA_nonlinsol, *maxiters);
    break;
  case FCMIX_ARKODE:
    if (!F2C_ARKODE_nonlinsol) {
      *ier = -1;
      return;
    }
    *ier = SUNNonlinSolSetMaxIters(F2C_ARKODE_nonlinsol, *maxiters);
    break;
  default:
    *ier = -1;
  }
}
