/* theta_esdirk2 */

#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

void _theta_esdirk2_res(sunrealtype* phi, const sunrealtype* th, const sunrealtype* rhs)
{
  phi[0] = RCONST(-1.0) * rhs[0] + th[0];
}

void _theta_esdirk2_jac(sunrealtype* phi_J, const sunrealtype* th)
{
  phi_J[0] = RCONST(1.0);
}

const sunindextype _theta_esdirk2_btable_ns = 2;

void _theta_esdirk2_btable_A(sunrealtype* A, const sunrealtype* th)
{
  A[0] = RCONST(0.0);
  A[1] = RCONST(0.0);
  A[2] = RCONST(1.0) + RCONST(-1.0) * th[0];
  A[3] = th[0];
}

void _theta_esdirk2_btable_b(sunrealtype* b, const sunrealtype* th)
{
  b[0] = RCONST(1.0) + RCONST(-1.0) * th[0];
  b[1] = th[0];
}

void _theta_esdirk2_btable_c(sunrealtype* c, const sunrealtype* th)
{
  c[0] = RCONST(0.0);
  c[1] = RCONST(1.0);
}

void _theta_esdirk2_guess(sunrealtype* guess)
{
  guess[0] = RCONST(0.5);
}

/* theta_sdirk2 */

void _theta_sdirk2_res(sunrealtype* phi, const sunrealtype* th, const sunrealtype* rhs)
{
  phi[0] = RCONST(-1.0) * rhs[0] + RCONST(2.0) * th[0] + RCONST(-1.0) * SUNRpowerI(th[0], 2);
}

void _theta_sdirk2_jac(sunrealtype* phi_J, const sunrealtype* th)
{
  phi_J[0] = RCONST(2.0) + RCONST(-2.0) * th[0];
}

const sunindextype _theta_sdirk2_btable_ns = 2;

void _theta_sdirk2_btable_A(sunrealtype* A, const sunrealtype* th)
{
  A[0] = th[0];
  A[1] = RCONST(0.0);
  A[2] = RCONST(1.0) + RCONST(-1.0) * th[0];
  A[3] = th[0];
}

void _theta_sdirk2_btable_b(sunrealtype* b, const sunrealtype* th)
{
  b[0] = RCONST(1.0) + RCONST(-1.0) * th[0];
  b[1] = th[0];
}

void _theta_sdirk2_btable_c(sunrealtype* c, const sunrealtype* th)
{
  c[0] = th[0];
  c[1] = RCONST(1.0);
}

void _theta_sdirk2_guess(sunrealtype* guess)
{
  guess[0] = RCONST(0.2928932188134524);
}

/* theta_esdirk3 */

void _theta_esdirk3_res(sunrealtype* phi, const sunrealtype* th, const sunrealtype* rhs)
{
  phi[0] = RCONST(-1.0) * rhs[0] + th[0] + th[1] + RCONST(-1.0) * th[0] * th[1] + RCONST(-1.0) * th[1] * th[2];
  phi[1] = RCONST(-1.0) * rhs[1] + th[0] + SUNRpowerI(th[1], 2) + RCONST(-1.0) * th[0] * SUNRpowerI(th[1], 2) + RCONST(-1.0) * SUNRpowerI(th[1], 2) * th[2];
  phi[2] = RCONST(-1.0) * rhs[2] + SUNRpowerI(th[0], 2) + RCONST(2.0) * th[0] * th[1] + RCONST(-2.0) * SUNRpowerI(th[0], 2) * th[1] + RCONST(-2.0) * th[0] * th[1] * th[2];
}

void _theta_esdirk3_jac(sunrealtype* phi_J, const sunrealtype* th)
{
  phi_J[0] = RCONST(1.0) + RCONST(-1.0) * th[1];
  phi_J[1] = RCONST(1.0) + RCONST(-1.0) * SUNRpowerI(th[1], 2);
  phi_J[2] = RCONST(2.0) * th[0] + RCONST(2.0) * th[1] + RCONST(-4.0) * th[0] * th[1] + RCONST(-2.0) * th[1] * th[2];
  phi_J[3] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2];
  phi_J[4] = RCONST(2.0) * th[1] + RCONST(-2.0) * th[0] * th[1] + RCONST(-2.0) * th[1] * th[2];
  phi_J[5] = RCONST(2.0) * th[0] + RCONST(-2.0) * SUNRpowerI(th[0], 2) + RCONST(-2.0) * th[0] * th[2];
  phi_J[6] = RCONST(-1.0) * th[1];
  phi_J[7] = RCONST(-1.0) * SUNRpowerI(th[1], 2);
  phi_J[8] = RCONST(-2.0) * th[0] * th[1];
}

const sunindextype _theta_esdirk3_btable_ns = 3;

void _theta_esdirk3_btable_A(sunrealtype* A, const sunrealtype* th)
{
  A[0] = RCONST(0.0);
  A[1] = RCONST(0.0);
  A[2] = RCONST(0.0);
  A[3] = RCONST(-1.0) * th[0] + th[1];
  A[4] = th[0];
  A[5] = RCONST(0.0);
  A[6] = th[2];
  A[7] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2];
  A[8] = th[0];
}

void _theta_esdirk3_btable_b(sunrealtype* b, const sunrealtype* th)
{
  b[0] = th[2];
  b[1] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2];
  b[2] = th[0];
}

void _theta_esdirk3_btable_c(sunrealtype* c, const sunrealtype* th)
{
  c[0] = RCONST(0.0);
  c[1] = th[1];
  c[2] = RCONST(1.0);
}

void _theta_esdirk3_guess(sunrealtype* guess)
{
  guess[0] = RCONST(0.21132486540518713);
  guess[1] = RCONST(0.4226497308103742);
  guess[2] = RCONST(0.1056624327025936);
}

/* theta_esdirk43 */

void _theta_esdirk43_res(sunrealtype* phi, const sunrealtype* th, const sunrealtype* rhs)
{
  phi[0] = RCONST(-0.5715450000000001) + RCONST(-1.0) * rhs[0] + RCONST(0.622) * th[0] + RCONST(0.768) * th[2];
  phi[1] = RCONST(-0.65499057) + RCONST(-1.0) * rhs[1] + RCONST(0.857116) * th[0] + RCONST(0.589824) * th[2];
  phi[2] = RCONST(-1.0) * rhs[2] + RCONST(-2.268594) * th[0] + RCONST(1.125504) * th[1] + RCONST(1.012) * SUNRpowerI(th[0], 2) + RCONST(-0.768) * th[0] * th[1] + RCONST(1.536) * th[0] * th[2];
}

void _theta_esdirk43_jac(sunrealtype* phi_J, const sunrealtype* th)
{
  phi_J[0] = RCONST(0.622);
  phi_J[1] = RCONST(0.857116);
  phi_J[2] = RCONST(-2.268594) + RCONST(2.024) * th[0] + RCONST(-0.768) * th[1] + RCONST(1.536) * th[2];
  phi_J[3] = RCONST(0.0);
  phi_J[4] = RCONST(0.0);
  phi_J[5] = RCONST(1.125504) + RCONST(-0.768) * th[0];
  phi_J[6] = RCONST(0.768);
  phi_J[7] = RCONST(0.589824);
  phi_J[8] = RCONST(1.536) * th[0];
}

const sunindextype _theta_esdirk43_btable_ns = 4;

void _theta_esdirk43_btable_A(sunrealtype* A, const sunrealtype* th)
{
  A[0] = RCONST(0.0);
  A[1] = RCONST(0.0);
  A[2] = RCONST(0.0);
  A[3] = RCONST(0.0);
  A[4] = RCONST(0.768) + RCONST(-1.0) * th[0];
  A[5] = th[0];
  A[6] = RCONST(0.0);
  A[7] = RCONST(0.0);
  A[8] = RCONST(0.378) + RCONST(-1.0) * th[1];
  A[9] = RCONST(-1.0) * th[0] + th[1];
  A[10] = th[0];
  A[11] = RCONST(0.0);
  A[12] = RCONST(1.0) + RCONST(-1.0) * th[2];
  A[13] = RCONST(-1.4655) + th[2];
  A[14] = RCONST(1.4655) + RCONST(-1.0) * th[0];
  A[15] = th[0];
}

void _theta_esdirk43_btable_b(sunrealtype* b, const sunrealtype* th)
{
  b[0] = RCONST(1.0) + RCONST(-1.0) * th[2];
  b[1] = RCONST(-1.4655) + th[2];
  b[2] = RCONST(1.4655) + RCONST(-1.0) * th[0];
  b[3] = th[0];
}

void _theta_esdirk43_btable_c(sunrealtype* c, const sunrealtype* th)
{
  c[0] = RCONST(0.0);
  c[1] = RCONST(0.768);
  c[2] = RCONST(0.378);
  c[3] = RCONST(1.0);
}

void _theta_esdirk43_guess(sunrealtype* guess)
{
  guess[0] = RCONST(0.435866521508459);
  guess[1] = RCONST(0.435866521508459);
  guess[2] = RCONST(0.2576482460664272);
}

int _theta_esdirk43_altguess(sunrealtype* guess, const sunindextype resets)
{
  int flag;
  switch (resets)
  {
  case 0:
    _theta_esdirk43_guess(guess);
    flag = SUNTRUE;
    break;

  case 1:
    guess[0] = RCONST(0.4362975330768013);
    guess[1] = RCONST(0.3357974919586735);
    guess[2] = RCONST(1.0418853312841532);
    flag = SUNTRUE;
    break;

  default:
    _theta_esdirk43_guess(guess);
    flag = SUNFALSE;
    break;
  }
  return flag;
}


/* theta_sdirk3 */

void _theta_sdirk3_res(sunrealtype* phi, const sunrealtype* th, const sunrealtype* rhs)
{
  phi[0] = RCONST(-1.0) * rhs[0] + RCONST(2.0) * th[0] + RCONST(-1.0) * SUNRpowerI(th[0], 2) + RCONST(-1.0) * th[0] * th[2] + th[1] * th[2];
  phi[1] = RCONST(-1.0) * rhs[1] + th[0] + SUNRpowerI(th[0], 2) + RCONST(-1.0) * SUNRpowerI(th[0], 3) + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[2] + SUNRpowerI(th[1], 2) * th[2];
  phi[2] = RCONST(-1.0) * rhs[2] + RCONST(3.0) * SUNRpowerI(th[0], 2) + RCONST(-2.0) * SUNRpowerI(th[0], 3) + RCONST(-3.0) * SUNRpowerI(th[0], 2) * th[2] + RCONST(3.0) * th[0] * th[1] * th[2];
}

void _theta_sdirk3_jac(sunrealtype* phi_J, const sunrealtype* th)
{
  phi_J[0] = RCONST(2.0) + RCONST(-2.0) * th[0] + RCONST(-1.0) * th[2];
  phi_J[1] = RCONST(1.0) + RCONST(2.0) * th[0] + RCONST(-3.0) * SUNRpowerI(th[0], 2) + RCONST(-2.0) * th[0] * th[2];
  phi_J[2] = RCONST(6.0) * th[0] + RCONST(-6.0) * SUNRpowerI(th[0], 2) + RCONST(-6.0) * th[0] * th[2] + RCONST(3.0) * th[1] * th[2];
  phi_J[3] = th[2];
  phi_J[4] = RCONST(2.0) * th[1] * th[2];
  phi_J[5] = RCONST(3.0) * th[0] * th[2];
  phi_J[6] = RCONST(-1.0) * th[0] + th[1];
  phi_J[7] = RCONST(-1.0) * SUNRpowerI(th[0], 2) + SUNRpowerI(th[1], 2);
  phi_J[8] = RCONST(-3.0) * SUNRpowerI(th[0], 2) + RCONST(3.0) * th[0] * th[1];
}

const sunindextype _theta_sdirk3_btable_ns = 3;

void _theta_sdirk3_btable_A(sunrealtype* A, const sunrealtype* th)
{
  A[0] = th[0];
  A[1] = RCONST(0.0);
  A[2] = RCONST(0.0);
  A[3] = RCONST(-1.0) * th[0] + th[1];
  A[4] = th[0];
  A[5] = RCONST(0.0);
  A[6] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2];
  A[7] = th[2];
  A[8] = th[0];
}

void _theta_sdirk3_btable_b(sunrealtype* b, const sunrealtype* th)
{
  b[0] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2];
  b[1] = th[2];
  b[2] = th[0];
}

void _theta_sdirk3_btable_c(sunrealtype* c, const sunrealtype* th)
{
  c[0] = th[0];
  c[1] = th[1];
  c[2] = RCONST(1.0);
}

void _theta_sdirk3_guess(sunrealtype* guess)
{
  guess[0] = RCONST(0.4358665215);
  guess[1] = RCONST(0.71793326075);
  guess[2] = RCONST(-0.6443631706532353);
}

int _theta_sdirk3_altguess(sunrealtype* guess, const sunindextype resets)
{
  int flag;
  switch (resets)
  {
  case 0:
    _theta_sdirk3_guess(guess);
    flag = SUNTRUE;
    break;

  case 1:
    guess[0] = RCONST(0.104356);
    guess[1] = RCONST(0.896898);
    guess[2] = RCONST(0.381277);
    flag = SUNTRUE;
    break;

  default:
    _theta_sdirk3_guess(guess);
    flag = SUNFALSE;
    break;
  }
  return flag;
}


/* theta_qesdirk4 */

void _theta_qesdirk4_res(sunrealtype* phi, const sunrealtype* th, const sunrealtype* rhs)
{
  phi[0] = RCONST(-1.0) * rhs[0] + th[0] + RCONST(0.5665418153404946) * th[5] + RCONST(0.955934110328081) * th[6] + th[0] * th[4];
  phi[1] = RCONST(-1.0) * rhs[1] + th[0] + RCONST(0.3209696285293031) * th[5] + RCONST(0.9138100232887396) * th[6] + SUNRpowerI(th[0], 2) * th[4];
  phi[2] = RCONST(-1.0) * rhs[2] + SUNRpowerI(th[0], 2) + RCONST(1.1330836306809893) * th[0] * th[5] + RCONST(1.911868220656162) * th[0] * th[6] + RCONST(0.5665418153404946) * th[3] * th[6] + RCONST(1.5) * SUNRpowerI(th[0], 2) * th[4] + th[0] * th[1] * th[5] + th[0] * th[2] * th[6];
  phi[3] = RCONST(-1.0) * rhs[3] + th[0] + RCONST(0.1818427160161556) * th[5] + RCONST(0.8735421716214042) * th[6] + SUNRpowerI(th[0], 3) * th[4];
  phi[4] = RCONST(-1.0) * rhs[4] + SUNRpowerI(th[0], 2) + RCONST(0.8875114438697977) * th[0] * th[5] + RCONST(1.8697441336168206) * th[0] * th[6] + RCONST(0.5415766462111716) * th[3] * th[6] + SUNRpowerI(th[0], 2) * th[4] + RCONST(0.5665418153404946) * th[0] * th[1] * th[5] + RCONST(0.955934110328081) * th[0] * th[2] * th[6] + RCONST(0.5) * SUNRpowerI(th[0], 3) * th[4];
  phi[5] = RCONST(-1.0) * rhs[5] + SUNRpowerI(th[0], 2) + RCONST(0.6419392570586062) * th[0] * th[5] + RCONST(1.8276200465774792) * th[0] * th[6] + RCONST(0.3209696285293031) * th[3] * th[6] + RCONST(1.5) * SUNRpowerI(th[0], 3) * th[4] + SUNRpowerI(th[0], 2) * th[1] * th[5] + SUNRpowerI(th[0], 2) * th[2] * th[6];
  phi[6] = RCONST(-1.0) * rhs[6] + SUNRpowerI(th[0], 3) + RCONST(1.699625446021484) * SUNRpowerI(th[0], 2) * th[5] + RCONST(2.867802330984243) * SUNRpowerI(th[0], 2) * th[6] + RCONST(1.699625446021484) * th[0] * th[3] * th[6] + RCONST(1.75) * SUNRpowerI(th[0], 3) * th[4] + RCONST(2.5) * SUNRpowerI(th[0], 2) * th[1] * th[5] + RCONST(2.5) * SUNRpowerI(th[0], 2) * th[2] * th[6] + th[0] * th[1] * th[3] * th[6];
}

void _theta_qesdirk4_jac(sunrealtype* phi_J, const sunrealtype* th)
{
  phi_J[0] = RCONST(1.0) + th[4];
  phi_J[1] = RCONST(1.0) + RCONST(2.0) * th[0] * th[4];
  phi_J[2] = RCONST(2.0) * th[0] + RCONST(1.1330836306809893) * th[5] + RCONST(1.911868220656162) * th[6] + RCONST(3.0) * th[0] * th[4] + th[1] * th[5] + th[2] * th[6];
  phi_J[3] = RCONST(1.0) + RCONST(3.0) * SUNRpowerI(th[0], 2) * th[4];
  phi_J[4] = RCONST(2.0) * th[0] + RCONST(0.8875114438697977) * th[5] + RCONST(1.8697441336168206) * th[6] + RCONST(2.0) * th[0] * th[4] + RCONST(0.5665418153404946) * th[1] * th[5] + RCONST(0.955934110328081) * th[2] * th[6] + RCONST(1.5) * SUNRpowerI(th[0], 2) * th[4];
  phi_J[5] = RCONST(2.0) * th[0] + RCONST(0.6419392570586062) * th[5] + RCONST(1.8276200465774792) * th[6] + RCONST(4.5) * SUNRpowerI(th[0], 2) * th[4] + RCONST(2.0) * th[0] * th[1] * th[5] + RCONST(2.0) * th[0] * th[2] * th[6];
  phi_J[6] = RCONST(3.0) * SUNRpowerI(th[0], 2) + RCONST(3.399250892042968) * th[0] * th[5] + RCONST(5.735604661968486) * th[0] * th[6] + RCONST(1.699625446021484) * th[3] * th[6] + RCONST(5.25) * SUNRpowerI(th[0], 2) * th[4] + RCONST(5.0) * th[0] * th[1] * th[5] + RCONST(5.0) * th[0] * th[2] * th[6] + th[1] * th[3] * th[6];
  phi_J[7] = RCONST(0.0);
  phi_J[8] = RCONST(0.0);
  phi_J[9] = th[0] * th[5];
  phi_J[10] = RCONST(0.0);
  phi_J[11] = RCONST(0.5665418153404946) * th[0] * th[5];
  phi_J[12] = SUNRpowerI(th[0], 2) * th[5];
  phi_J[13] = RCONST(2.5) * SUNRpowerI(th[0], 2) * th[5] + th[0] * th[3] * th[6];
  phi_J[14] = RCONST(0.0);
  phi_J[15] = RCONST(0.0);
  phi_J[16] = th[0] * th[6];
  phi_J[17] = RCONST(0.0);
  phi_J[18] = RCONST(0.955934110328081) * th[0] * th[6];
  phi_J[19] = SUNRpowerI(th[0], 2) * th[6];
  phi_J[20] = RCONST(2.5) * SUNRpowerI(th[0], 2) * th[6];
  phi_J[21] = RCONST(0.0);
  phi_J[22] = RCONST(0.0);
  phi_J[23] = RCONST(0.5665418153404946) * th[6];
  phi_J[24] = RCONST(0.0);
  phi_J[25] = RCONST(0.5415766462111716) * th[6];
  phi_J[26] = RCONST(0.3209696285293031) * th[6];
  phi_J[27] = RCONST(1.699625446021484) * th[0] * th[6] + th[0] * th[1] * th[6];
  phi_J[28] = th[0];
  phi_J[29] = SUNRpowerI(th[0], 2);
  phi_J[30] = RCONST(1.5) * SUNRpowerI(th[0], 2);
  phi_J[31] = SUNRpowerI(th[0], 3);
  phi_J[32] = SUNRpowerI(th[0], 2) + RCONST(0.5) * SUNRpowerI(th[0], 3);
  phi_J[33] = RCONST(1.5) * SUNRpowerI(th[0], 3);
  phi_J[34] = RCONST(1.75) * SUNRpowerI(th[0], 3);
  phi_J[35] = RCONST(0.5665418153404946);
  phi_J[36] = RCONST(0.3209696285293031);
  phi_J[37] = RCONST(1.1330836306809893) * th[0] + th[0] * th[1];
  phi_J[38] = RCONST(0.1818427160161556);
  phi_J[39] = RCONST(0.8875114438697977) * th[0] + RCONST(0.5665418153404946) * th[0] * th[1];
  phi_J[40] = RCONST(0.6419392570586062) * th[0] + SUNRpowerI(th[0], 2) * th[1];
  phi_J[41] = RCONST(1.699625446021484) * SUNRpowerI(th[0], 2) + RCONST(2.5) * SUNRpowerI(th[0], 2) * th[1];
  phi_J[42] = RCONST(0.955934110328081);
  phi_J[43] = RCONST(0.9138100232887396);
  phi_J[44] = RCONST(1.911868220656162) * th[0] + RCONST(0.5665418153404946) * th[3] + th[0] * th[2];
  phi_J[45] = RCONST(0.8735421716214042);
  phi_J[46] = RCONST(1.8697441336168206) * th[0] + RCONST(0.5415766462111716) * th[3] + RCONST(0.955934110328081) * th[0] * th[2];
  phi_J[47] = RCONST(1.8276200465774792) * th[0] + RCONST(0.3209696285293031) * th[3] + SUNRpowerI(th[0], 2) * th[2];
  phi_J[48] = RCONST(2.867802330984243) * SUNRpowerI(th[0], 2) + RCONST(1.699625446021484) * th[0] * th[3] + RCONST(2.5) * SUNRpowerI(th[0], 2) * th[2] + th[0] * th[1] * th[3];
}

const sunindextype _theta_qesdirk4_btable_ns = 5;

void _theta_qesdirk4_btable_A(sunrealtype* A, const sunrealtype* th)
{
  A[0] = RCONST(0.0);
  A[1] = RCONST(0.0);
  A[2] = RCONST(0.0);
  A[3] = RCONST(0.0);
  A[4] = RCONST(0.0);
  A[5] = RCONST(0.5) * th[0];
  A[6] = RCONST(0.5) * th[0];
  A[7] = RCONST(0.0);
  A[8] = RCONST(0.0);
  A[9] = RCONST(0.0);
  A[10] = RCONST(0.5665418153404946) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[1];
  A[11] = th[1];
  A[12] = th[0];
  A[13] = RCONST(0.0);
  A[14] = RCONST(0.0);
  A[15] = RCONST(0.955934110328081) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2] + RCONST(-1.0) * th[3];
  A[16] = th[2];
  A[17] = th[3];
  A[18] = th[0];
  A[19] = RCONST(0.0);
  A[20] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[4] + RCONST(-1.0) * th[5] + RCONST(-1.0) * th[6];
  A[21] = th[4];
  A[22] = th[5];
  A[23] = th[6];
  A[24] = th[0];
}

void _theta_qesdirk4_btable_b(sunrealtype* b, const sunrealtype* th)
{
  b[0] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[4] + RCONST(-1.0) * th[5] + RCONST(-1.0) * th[6];
  b[1] = th[4];
  b[2] = th[5];
  b[3] = th[6];
  b[4] = th[0];
}

void _theta_qesdirk4_btable_c(sunrealtype* c, const sunrealtype* th)
{
  c[0] = RCONST(0.0);
  c[1] = th[0];
  c[2] = RCONST(0.5665418153404946);
  c[3] = RCONST(0.955934110328081);
  c[4] = RCONST(1.0);
}

void _theta_qesdirk4_guess(sunrealtype* guess)
{
  guess[0] = RCONST(0.5065274202451187);
  guess[1] = RCONST(-0.1264307159555053);
  guess[2] = RCONST(146.3451690077889);
  guess[3] = RCONST(-305.0457869546525);
  guess[4] = RCONST(-0.31902898451398665);
  guess[5] = RCONST(0.6494281724470768);
  guess[6] = RCONST(0.001108203953233309);
}

/* theta_sdirk4 */

void _theta_sdirk4_res(sunrealtype* phi, const sunrealtype* th, const sunrealtype* rhs)
{
  phi[0] = RCONST(-1.0) * rhs[0] + RCONST(2.0) * th[0] + RCONST(0.757) * th[4] + RCONST(0.572) * th[5] + RCONST(0.234) * th[6] + RCONST(-1.0) * SUNRpowerI(th[0], 2) + RCONST(-1.0) * th[0] * th[4] + RCONST(-1.0) * th[0] * th[5] + RCONST(-1.0) * th[0] * th[6];
  phi[1] = RCONST(-1.0) * rhs[1] + th[0] + RCONST(0.573049) * th[4] + RCONST(0.3271839999999999) * th[5] + RCONST(0.054756000000000006) * th[6] + SUNRpowerI(th[0], 2) + RCONST(-1.0) * SUNRpowerI(th[0], 3) + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[4] + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[5] + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[6];
  phi[2] = RCONST(-1.0) * rhs[2] + RCONST(3.0) * SUNRpowerI(th[0], 2) + RCONST(2.271) * th[0] * th[4] + RCONST(1.7159999999999997) * th[0] * th[5] + RCONST(0.7020000000000001) * th[0] * th[6] + RCONST(0.757) * th[1] * th[5] + RCONST(0.757) * th[2] * th[6] + RCONST(0.572) * th[3] * th[6] + RCONST(-2.0) * SUNRpowerI(th[0], 3) + RCONST(-3.0) * SUNRpowerI(th[0], 2) * th[4] + RCONST(-3.0) * SUNRpowerI(th[0], 2) * th[5] + RCONST(-3.0) * SUNRpowerI(th[0], 2) * th[6] + RCONST(-1.0) * th[0] * th[1] * th[5] + RCONST(-1.0) * th[0] * th[2] * th[6] + RCONST(-1.0) * th[0] * th[3] * th[6];
  phi[3] = RCONST(-1.0) * rhs[3] + th[0] + RCONST(0.43379809300000005) * th[4] + RCONST(0.18714924799999993) * th[5] + RCONST(0.012812904000000002) * th[6] + SUNRpowerI(th[0], 3) + RCONST(-1.0) * SUNRpowerI(th[0], 4) + RCONST(-1.0) * SUNRpowerI(th[0], 3) * th[4] + RCONST(-1.0) * SUNRpowerI(th[0], 3) * th[5] + RCONST(-1.0) * SUNRpowerI(th[0], 3) * th[6];
  phi[4] = RCONST(-1.0) * rhs[4] + RCONST(2.0) * SUNRpowerI(th[0], 2) + RCONST(1.903098) * th[0] * th[4] + RCONST(1.226368) * th[0] * th[5] + RCONST(0.34351200000000004) * th[0] * th[6] + RCONST(0.43300399999999994) * th[1] * th[5] + RCONST(0.17713800000000002) * th[2] * th[6] + RCONST(0.133848) * th[3] * th[6] + RCONST(-1.7570000000000001) * SUNRpowerI(th[0], 2) * th[4] + RCONST(-1.572) * SUNRpowerI(th[0], 2) * th[5] + RCONST(-1.234) * SUNRpowerI(th[0], 2) * th[6] + RCONST(-0.572) * th[0] * th[1] * th[5] + RCONST(-0.234) * th[0] * th[2] * th[6] + RCONST(-0.234) * th[0] * th[3] * th[6] + RCONST(-1.0) * SUNRpowerI(th[0], 4) + RCONST(-1.0) * SUNRpowerI(th[0], 3) * th[4] + RCONST(-1.0) * SUNRpowerI(th[0], 3) * th[5] + RCONST(-1.0) * SUNRpowerI(th[0], 3) * th[6];
  phi[5] = RCONST(-1.0) * rhs[5] + SUNRpowerI(th[0], 2) + RCONST(1.146098) * th[0] * th[4] + RCONST(0.6543679999999998) * th[0] * th[5] + RCONST(0.10951200000000001) * th[0] * th[6] + RCONST(0.573049) * th[1] * th[5] + RCONST(0.573049) * th[2] * th[6] + RCONST(0.3271839999999999) * th[3] * th[6] + RCONST(2.0) * SUNRpowerI(th[0], 3) + RCONST(0.757) * SUNRpowerI(th[0], 2) * th[4] + RCONST(0.572) * SUNRpowerI(th[0], 2) * th[5] + RCONST(0.234) * SUNRpowerI(th[0], 2) * th[6] + RCONST(-2.0) * SUNRpowerI(th[0], 4) + RCONST(-3.0) * SUNRpowerI(th[0], 3) * th[4] + RCONST(-3.0) * SUNRpowerI(th[0], 3) * th[5] + RCONST(-3.0) * SUNRpowerI(th[0], 3) * th[6] + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[1] * th[5] + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[2] * th[6] + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[3] * th[6];
  phi[6] = RCONST(-1.0) * rhs[6] + RCONST(4.0) * SUNRpowerI(th[0], 3) + RCONST(4.542) * SUNRpowerI(th[0], 2) * th[4] + RCONST(3.4319999999999995) * SUNRpowerI(th[0], 2) * th[5] + RCONST(1.4040000000000001) * SUNRpowerI(th[0], 2) * th[6] + RCONST(3.028) * th[0] * th[1] * th[5] + RCONST(3.028) * th[0] * th[2] * th[6] + RCONST(2.288) * th[0] * th[3] * th[6] + RCONST(0.757) * th[1] * th[3] * th[6] + RCONST(-3.0) * SUNRpowerI(th[0], 4) + RCONST(-6.0) * SUNRpowerI(th[0], 3) * th[4] + RCONST(-6.0) * SUNRpowerI(th[0], 3) * th[5] + RCONST(-6.0) * SUNRpowerI(th[0], 3) * th[6] + RCONST(-4.0) * SUNRpowerI(th[0], 2) * th[1] * th[5] + RCONST(-4.0) * SUNRpowerI(th[0], 2) * th[2] * th[6] + RCONST(-4.0) * SUNRpowerI(th[0], 2) * th[3] * th[6] + RCONST(-1.0) * th[0] * th[1] * th[3] * th[6];
}

void _theta_sdirk4_jac(sunrealtype* phi_J, const sunrealtype* th)
{
  phi_J[0] = RCONST(2.0) + RCONST(-2.0) * th[0] + RCONST(-1.0) * th[4] + RCONST(-1.0) * th[5] + RCONST(-1.0) * th[6];
  phi_J[1] = RCONST(1.0) + RCONST(2.0) * th[0] + RCONST(-3.0) * SUNRpowerI(th[0], 2) + RCONST(-2.0) * th[0] * th[4] + RCONST(-2.0) * th[0] * th[5] + RCONST(-2.0) * th[0] * th[6];
  phi_J[2] = RCONST(6.0) * th[0] + RCONST(2.271) * th[4] + RCONST(1.7159999999999997) * th[5] + RCONST(0.7020000000000001) * th[6] + RCONST(-6.0) * SUNRpowerI(th[0], 2) + RCONST(-6.0) * th[0] * th[4] + RCONST(-6.0) * th[0] * th[5] + RCONST(-6.0) * th[0] * th[6] + RCONST(-1.0) * th[1] * th[5] + RCONST(-1.0) * th[2] * th[6] + RCONST(-1.0) * th[3] * th[6];
  phi_J[3] = RCONST(1.0) + RCONST(3.0) * SUNRpowerI(th[0], 2) + RCONST(-4.0) * SUNRpowerI(th[0], 3) + RCONST(-3.0) * SUNRpowerI(th[0], 2) * th[4] + RCONST(-3.0) * SUNRpowerI(th[0], 2) * th[5] + RCONST(-3.0) * SUNRpowerI(th[0], 2) * th[6];
  phi_J[4] = RCONST(4.0) * th[0] + RCONST(1.903098) * th[4] + RCONST(1.226368) * th[5] + RCONST(0.34351200000000004) * th[6] + RCONST(-3.5140000000000002) * th[0] * th[4] + RCONST(-3.144) * th[0] * th[5] + RCONST(-2.468) * th[0] * th[6] + RCONST(-0.572) * th[1] * th[5] + RCONST(-0.234) * th[2] * th[6] + RCONST(-0.234) * th[3] * th[6] + RCONST(-4.0) * SUNRpowerI(th[0], 3) + RCONST(-3.0) * SUNRpowerI(th[0], 2) * th[4] + RCONST(-3.0) * SUNRpowerI(th[0], 2) * th[5] + RCONST(-3.0) * SUNRpowerI(th[0], 2) * th[6];
  phi_J[5] = RCONST(2.0) * th[0] + RCONST(1.146098) * th[4] + RCONST(0.6543679999999998) * th[5] + RCONST(0.10951200000000001) * th[6] + RCONST(6.0) * SUNRpowerI(th[0], 2) + RCONST(1.514) * th[0] * th[4] + RCONST(1.144) * th[0] * th[5] + RCONST(0.468) * th[0] * th[6] + RCONST(-8.0) * SUNRpowerI(th[0], 3) + RCONST(-9.0) * SUNRpowerI(th[0], 2) * th[4] + RCONST(-9.0) * SUNRpowerI(th[0], 2) * th[5] + RCONST(-9.0) * SUNRpowerI(th[0], 2) * th[6] + RCONST(-2.0) * th[0] * th[1] * th[5] + RCONST(-2.0) * th[0] * th[2] * th[6] + RCONST(-2.0) * th[0] * th[3] * th[6];
  phi_J[6] = RCONST(12.0) * SUNRpowerI(th[0], 2) + RCONST(9.084) * th[0] * th[4] + RCONST(6.863999999999999) * th[0] * th[5] + RCONST(2.8080000000000003) * th[0] * th[6] + RCONST(3.028) * th[1] * th[5] + RCONST(3.028) * th[2] * th[6] + RCONST(2.288) * th[3] * th[6] + RCONST(-12.0) * SUNRpowerI(th[0], 3) + RCONST(-18.0) * SUNRpowerI(th[0], 2) * th[4] + RCONST(-18.0) * SUNRpowerI(th[0], 2) * th[5] + RCONST(-18.0) * SUNRpowerI(th[0], 2) * th[6] + RCONST(-8.0) * th[0] * th[1] * th[5] + RCONST(-8.0) * th[0] * th[2] * th[6] + RCONST(-8.0) * th[0] * th[3] * th[6] + RCONST(-1.0) * th[1] * th[3] * th[6];
  phi_J[7] = RCONST(0.0);
  phi_J[8] = RCONST(0.0);
  phi_J[9] = RCONST(0.757) * th[5] + RCONST(-1.0) * th[0] * th[5];
  phi_J[10] = RCONST(0.0);
  phi_J[11] = RCONST(0.43300399999999994) * th[5] + RCONST(-0.572) * th[0] * th[5];
  phi_J[12] = RCONST(0.573049) * th[5] + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[5];
  phi_J[13] = RCONST(3.028) * th[0] * th[5] + RCONST(0.757) * th[3] * th[6] + RCONST(-4.0) * SUNRpowerI(th[0], 2) * th[5] + RCONST(-1.0) * th[0] * th[3] * th[6];
  phi_J[14] = RCONST(0.0);
  phi_J[15] = RCONST(0.0);
  phi_J[16] = RCONST(0.757) * th[6] + RCONST(-1.0) * th[0] * th[6];
  phi_J[17] = RCONST(0.0);
  phi_J[18] = RCONST(0.17713800000000002) * th[6] + RCONST(-0.234) * th[0] * th[6];
  phi_J[19] = RCONST(0.573049) * th[6] + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[6];
  phi_J[20] = RCONST(3.028) * th[0] * th[6] + RCONST(-4.0) * SUNRpowerI(th[0], 2) * th[6];
  phi_J[21] = RCONST(0.0);
  phi_J[22] = RCONST(0.0);
  phi_J[23] = RCONST(0.572) * th[6] + RCONST(-1.0) * th[0] * th[6];
  phi_J[24] = RCONST(0.0);
  phi_J[25] = RCONST(0.133848) * th[6] + RCONST(-0.234) * th[0] * th[6];
  phi_J[26] = RCONST(0.3271839999999999) * th[6] + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[6];
  phi_J[27] = RCONST(2.288) * th[0] * th[6] + RCONST(0.757) * th[1] * th[6] + RCONST(-4.0) * SUNRpowerI(th[0], 2) * th[6] + RCONST(-1.0) * th[0] * th[1] * th[6];
  phi_J[28] = RCONST(0.757) + RCONST(-1.0) * th[0];
  phi_J[29] = RCONST(0.573049) + RCONST(-1.0) * SUNRpowerI(th[0], 2);
  phi_J[30] = RCONST(2.271) * th[0] + RCONST(-3.0) * SUNRpowerI(th[0], 2);
  phi_J[31] = RCONST(0.43379809300000005) + RCONST(-1.0) * SUNRpowerI(th[0], 3);
  phi_J[32] = RCONST(1.903098) * th[0] + RCONST(-1.7570000000000001) * SUNRpowerI(th[0], 2) + RCONST(-1.0) * SUNRpowerI(th[0], 3);
  phi_J[33] = RCONST(1.146098) * th[0] + RCONST(0.757) * SUNRpowerI(th[0], 2) + RCONST(-3.0) * SUNRpowerI(th[0], 3);
  phi_J[34] = RCONST(4.542) * SUNRpowerI(th[0], 2) + RCONST(-6.0) * SUNRpowerI(th[0], 3);
  phi_J[35] = RCONST(0.572) + RCONST(-1.0) * th[0];
  phi_J[36] = RCONST(0.3271839999999999) + RCONST(-1.0) * SUNRpowerI(th[0], 2);
  phi_J[37] = RCONST(1.7159999999999997) * th[0] + RCONST(0.757) * th[1] + RCONST(-3.0) * SUNRpowerI(th[0], 2) + RCONST(-1.0) * th[0] * th[1];
  phi_J[38] = RCONST(0.18714924799999993) + RCONST(-1.0) * SUNRpowerI(th[0], 3);
  phi_J[39] = RCONST(1.226368) * th[0] + RCONST(0.43300399999999994) * th[1] + RCONST(-1.572) * SUNRpowerI(th[0], 2) + RCONST(-0.572) * th[0] * th[1] + RCONST(-1.0) * SUNRpowerI(th[0], 3);
  phi_J[40] = RCONST(0.6543679999999998) * th[0] + RCONST(0.573049) * th[1] + RCONST(0.572) * SUNRpowerI(th[0], 2) + RCONST(-3.0) * SUNRpowerI(th[0], 3) + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[1];
  phi_J[41] = RCONST(3.4319999999999995) * SUNRpowerI(th[0], 2) + RCONST(3.028) * th[0] * th[1] + RCONST(-6.0) * SUNRpowerI(th[0], 3) + RCONST(-4.0) * SUNRpowerI(th[0], 2) * th[1];
  phi_J[42] = RCONST(0.234) + RCONST(-1.0) * th[0];
  phi_J[43] = RCONST(0.054756000000000006) + RCONST(-1.0) * SUNRpowerI(th[0], 2);
  phi_J[44] = RCONST(0.7020000000000001) * th[0] + RCONST(0.757) * th[2] + RCONST(0.572) * th[3] + RCONST(-3.0) * SUNRpowerI(th[0], 2) + RCONST(-1.0) * th[0] * th[2] + RCONST(-1.0) * th[0] * th[3];
  phi_J[45] = RCONST(0.012812904000000002) + RCONST(-1.0) * SUNRpowerI(th[0], 3);
  phi_J[46] = RCONST(0.34351200000000004) * th[0] + RCONST(0.17713800000000002) * th[2] + RCONST(0.133848) * th[3] + RCONST(-1.234) * SUNRpowerI(th[0], 2) + RCONST(-0.234) * th[0] * th[2] + RCONST(-0.234) * th[0] * th[3] + RCONST(-1.0) * SUNRpowerI(th[0], 3);
  phi_J[47] = RCONST(0.10951200000000001) * th[0] + RCONST(0.573049) * th[2] + RCONST(0.3271839999999999) * th[3] + RCONST(0.234) * SUNRpowerI(th[0], 2) + RCONST(-3.0) * SUNRpowerI(th[0], 3) + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[2] + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[3];
  phi_J[48] = RCONST(1.4040000000000001) * SUNRpowerI(th[0], 2) + RCONST(3.028) * th[0] * th[2] + RCONST(2.288) * th[0] * th[3] + RCONST(0.757) * th[1] * th[3] + RCONST(-6.0) * SUNRpowerI(th[0], 3) + RCONST(-4.0) * SUNRpowerI(th[0], 2) * th[2] + RCONST(-4.0) * SUNRpowerI(th[0], 2) * th[3] + RCONST(-1.0) * th[0] * th[1] * th[3];
}

const sunindextype _theta_sdirk4_btable_ns = 5;

void _theta_sdirk4_btable_A(sunrealtype* A, const sunrealtype* th)
{
  A[0] = th[0];
  A[1] = RCONST(0.0);
  A[2] = RCONST(0.0);
  A[3] = RCONST(0.0);
  A[4] = RCONST(0.0);
  A[5] = RCONST(0.757) + RCONST(-1.0) * th[0];
  A[6] = th[0];
  A[7] = RCONST(0.0);
  A[8] = RCONST(0.0);
  A[9] = RCONST(0.0);
  A[10] = RCONST(0.572) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[1];
  A[11] = th[1];
  A[12] = th[0];
  A[13] = RCONST(0.0);
  A[14] = RCONST(0.0);
  A[15] = RCONST(0.234) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2] + RCONST(-1.0) * th[3];
  A[16] = th[2];
  A[17] = th[3];
  A[18] = th[0];
  A[19] = RCONST(0.0);
  A[20] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[4] + RCONST(-1.0) * th[5] + RCONST(-1.0) * th[6];
  A[21] = th[4];
  A[22] = th[5];
  A[23] = th[6];
  A[24] = th[0];
}

void _theta_sdirk4_btable_b(sunrealtype* b, const sunrealtype* th)
{
  b[0] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[4] + RCONST(-1.0) * th[5] + RCONST(-1.0) * th[6];
  b[1] = th[4];
  b[2] = th[5];
  b[3] = th[6];
  b[4] = th[0];
}

void _theta_sdirk4_btable_c(sunrealtype* c, const sunrealtype* th)
{
  c[0] = th[0];
  c[1] = RCONST(0.757);
  c[2] = RCONST(0.572);
  c[3] = RCONST(0.234);
  c[4] = RCONST(1.0);
}

void _theta_sdirk4_guess(sunrealtype* guess)
{
  guess[0] = RCONST(0.2780787351155658);
  guess[1] = RCONST(-0.062252697221817865);
  guess[2] = RCONST(0.060267864764432834);
  guess[3] = RCONST(-0.08371700471114331);
  guess[4] = RCONST(-0.7640793347677586);
  guess[5] = RCONST(1.8115844527853433);
  guess[6] = RCONST(3.2977134950712017);
}

int _theta_sdirk4_altguess(sunrealtype* guess, const sunindextype resets)
{
  int flag;
  switch (resets)
  {
  case 0:
    _theta_sdirk4_guess(guess);
    flag = SUNTRUE;
    break;

  case 1:
    guess[0] = RCONST(0.25);
    guess[1] = RCONST(-0.04);
    guess[2] = RCONST(-0.05036764705882353);
    guess[3] = RCONST(0.027573529411764705);
    guess[4] = RCONST(-1.0208333333333333);
    guess[5] = RCONST(7.8125);
    guess[6] = RCONST(-7.083333333333333);
    flag = SUNTRUE;
    break;

  default:
    _theta_sdirk4_guess(guess);
    flag = SUNFALSE;
    break;
  }
  return flag;
}


