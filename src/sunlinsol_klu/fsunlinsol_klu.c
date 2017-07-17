/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
 * -----------------------------------------------------------------
 * This file (companion of fsunlinsol_klu.h) contains the
 * implementation needed for the Fortran initialization of klu
 * linear solver operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fsunlinsol_klu.h"

/* Define global linsol variables */

SUNLinearSolver F2C_CVODE_linsol;
SUNLinearSolver F2C_IDA_linsol;
SUNLinearSolver F2C_KINSOL_linsol;
SUNLinearSolver F2C_ARKODE_linsol;

/* Declarations of external global variables */

extern SUNMatrix F2C_CVODE_matrix;
extern SUNMatrix F2C_IDA_matrix;
extern SUNMatrix F2C_KINSOL_matrix;
extern SUNMatrix F2C_ARKODE_matrix;

extern N_Vector F2C_CVODE_vec;
extern N_Vector F2C_IDA_vec;
extern N_Vector F2C_KINSOL_vec;
extern N_Vector F2C_ARKODE_vec;

/* Fortran callable interfaces */

void FSUNKLU_INIT(int *code, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_linsol = NULL;
    F2C_CVODE_linsol = SUNKLULinearSolver(F2C_CVODE_vec,
                                          F2C_CVODE_matrix);
    if (F2C_CVODE_linsol == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_linsol = NULL;
    F2C_IDA_linsol = SUNKLULinearSolver(F2C_IDA_vec,
                                        F2C_IDA_matrix);
    if (F2C_IDA_linsol == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_linsol = NULL;
    F2C_KINSOL_linsol = SUNKLULinearSolver(F2C_KINSOL_vec,
                                           F2C_KINSOL_matrix);
    if (F2C_KINSOL_linsol == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    F2C_ARKODE_linsol = NULL;
    F2C_ARKODE_linsol = SUNKLULinearSolver(F2C_ARKODE_vec,
                                           F2C_ARKODE_matrix);
    if (F2C_ARKODE_linsol == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


void FSUNKLU_REINIT(int *code, long int *NNZ, int *reinit_type, int *ier)
{
  *ier = 0;

  sunindextype nnz = (sunindextype) *NNZ;
  switch(*code) {
  case FCMIX_CVODE:
    *ier = SUNKLUReInit(F2C_CVODE_linsol, F2C_CVODE_matrix,
                        nnz, *reinit_type);
    break;
  case FCMIX_IDA:
    *ier = SUNKLUReInit(F2C_IDA_linsol, F2C_IDA_matrix,
                        nnz, *reinit_type);
    break;
  case FCMIX_KINSOL:
    *ier = SUNKLUReInit(F2C_KINSOL_linsol, F2C_KINSOL_matrix,
                        nnz, *reinit_type);
    break;
  case FCMIX_ARKODE:
    *ier = SUNKLUReInit(F2C_ARKODE_linsol, F2C_ARKODE_matrix,
                        nnz, *reinit_type);
    break;
  default:
    *ier = -1;
  }
}


void FSUNKLU_SETORDERING(int *code, int *ordering_choice, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    *ier = SUNKLUSetOrdering(F2C_CVODE_linsol, *ordering_choice);
    break;
  case FCMIX_IDA:
    *ier = SUNKLUSetOrdering(F2C_IDA_linsol, *ordering_choice);
    break;
  case FCMIX_KINSOL:
    *ier = SUNKLUSetOrdering(F2C_KINSOL_linsol, *ordering_choice);
    break;
  case FCMIX_ARKODE:
    *ier = SUNKLUSetOrdering(F2C_ARKODE_linsol, *ordering_choice);
    break;
  default:
    *ier = -1;
  }
}
