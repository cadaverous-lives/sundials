% KINSOL, nonlinear system solver
%
% Functions
%   KINSetOptions - creates an options structure for KINSOL.
%   KINMalloc     - allocates and initializes memory for KINSOL.
%   KINSol        - solves the nonlinear problem.
%   KINGetStats   - returns statistics for the KINSOL solver
%   KINFree       - deallocates memory for the KINSOL solver.
%
% User-supplied function types   
%   KINsysFn         -  system function
%   KINDenseJacFn    -  dense Jacobian function
%   KINJacTimesVecFn -  Jacobian times vector function
%   KINPrecSetupFn   -  preconditioner setup function
%   KINPrecSolveFn   -  preconditioner solve function
%   KINGlocalFn      -  RHS approximation function (BBDPre)
%   KINGcomFn        -  communication function (BBDPre)
%
% Serial examples
%   mkinDiagon_kry - simple serial diagonal example
%   mkinTest_dns   - simple 2 dimensional example
%
% Parallel examples
%   mkinDiagon_kry_p - simple parallel diagonal example
% Use the mpirun function to run any of the parallel examples
%
% See also nvector, putils
