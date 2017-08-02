/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 *---------------------------------------------------------------
 * This module contains the routines necessary to interface with 
 * the ARKBBDPRE module and user-supplied Fortran routines.
 * The routines here call the generically named routines and 
 * providea standard interface to the C code of the ARKBBDPRE 
 * package.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "farkbbd.h"
#include <arkode/arkode_bbdpre.h>

/*=============================================================*/

/* Prototypes of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_GLOCFN(sunindextype *NLOC, realtype *T, 
			  realtype *YLOC, realtype *GLOC,
			  sunindextype *IPAR, realtype *RPAR,
			  int *ier);
  extern void FARK_COMMFN(sunindextype *NLOC, realtype *T, 
			  realtype *Y, sunindextype *IPAR, 
			  realtype *RPAR, int *ier);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface to C routine ARKBBDPrecInit; see farkbbd.h 
   for further details. */
void FARK_BBDINIT(sunindextype *Nloc, sunindextype *mudq,
                  sunindextype *mldq, sunindextype *mu,
                  sunindextype *ml, realtype* dqrely, 
		  int *ier)
{
  /* Notes: FARKgloc is a pointer to the ARKLocalFn function, 
     and FARKcfn is a pointer to the ARKCommFn function */
  *ier = ARKBBDPrecInit(ARK_arkodemem, *Nloc, *mudq, *mldq, 
			*mu, *ml, *dqrely, FARKgloc, FARKcfn);
  return; 
}

/*=============================================================*/

/* Fortran interface to C routine ARKBBDPrecReInit; see farkbbd.h 
   for further details. */
void FARK_BBDREINIT(sunindextype *mudq, sunindextype *mldq,
                    realtype* dqrely, int *ier)
{
  *ier = ARKBBDPrecReInit(ARK_arkodemem, *mudq, *mldq, *dqrely);
  return;
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKGLOCFN; see 
   farkbbd.h for further details. */
int FARKgloc(sunindextype Nloc, realtype t, N_Vector yloc, 
	     N_Vector gloc, void *user_data)
{
  int ier;
  realtype *yloc_data, *gloc_data;
  FARKUserData ARK_userdata;

  yloc_data = N_VGetArrayPointer(yloc);
  gloc_data = N_VGetArrayPointer(gloc);
  ARK_userdata = (FARKUserData) user_data;

  FARK_GLOCFN(&Nloc, &t, yloc_data, gloc_data, 
             ARK_userdata->ipar, ARK_userdata->rpar, &ier);
  return(ier);
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKCOMMFN; see 
   farkbbd.h for further details. */
int FARKcfn(sunindextype Nloc, realtype t, N_Vector y, void *user_data)
{
  int ier;
  realtype *yloc;
  FARKUserData ARK_userdata;

  yloc = N_VGetArrayPointer(y);
  ARK_userdata = (FARKUserData) user_data;
  FARK_COMMFN(&Nloc, &t, yloc, ARK_userdata->ipar, 
	      ARK_userdata->rpar, &ier);
  return(ier);
}

/*=============================================================*/

/* Fortran interface to C routines ARKBBDPrecGetWorkSpace and 
   ARKBBDPrecGetNumGfnEvals; see farkbbd.h for further details */
void FARK_BBDOPT(long int *lenrwbbd, long int *leniwbbd, 
		 long int *ngebbd)
{
  ARKBBDPrecGetWorkSpace(ARK_arkodemem, lenrwbbd, leniwbbd);
  ARKBBDPrecGetNumGfnEvals(ARK_arkodemem, ngebbd);
  return;
}

/*===============================================================
   EOF
===============================================================*/
