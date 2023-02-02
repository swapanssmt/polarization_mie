
#include <string>
#define _USE_MATH_DEFINES
#define MIETHEORY_MEX
#include <cmath>
#include <limits>
#include <inttypes.h>
#include <string>
#include <vector>

#include "mex.h"
#include "Array.hpp"
#include "ArrayMEX.hpp"
#include "MT.hpp"
#include "matrix.h"
//#include "arrays.hpp"
//#include "complex.hpp"
//#include "mie.hpp"
//#include "arrays.cpp"
//#include "complex.cpp"
//#include "mie.cpp"
#include "versionstring.h"
// Compiling (from MATLAB prompt):
//   mex MTmex.cpp
//
// To compile with OpenMP (multithread) support (from MATLAB prompt):
//   mex -DUSE_OMP MTmex.cpp CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
// Do not use OpenMP version if the MATLAB does not support the compiler used

void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
  mexPrintf("                 Mie-Theory\n");
  char infobuf[5012];
  version_string(infobuf);
  mexPrintf("%s",infobuf);
  
  if ((nrhs != 8) || (nlhs != 5))
  {
    mexPrintf("nrhs %i nlhs %i", nrhs, nlhs);
    mexErrMsgTxt("Syntax:\n [mus, s11, s12, s33, s43] = MTmex(radius, lambda, nre_med, nim_med, nre_p, nim_p, nangles, rho)\n");
  }
  mexPrintf("Initializing MT...\n");
  Array<double> radius, lambda, nre_med, nim_med, nre_p, nim_p, rho;
  Array<int_fast64_t> nangles;

  Convert_mxArray(prhs[0], radius);
  Convert_mxArray(prhs[1], lambda);
  Convert_mxArray(prhs[2], nre_med);
  Convert_mxArray(prhs[3], nim_med);
  Convert_mxArray(prhs[4], nre_p);
  Convert_mxArray(prhs[5], nim_p);   
  Convert_mxArray(prhs[6], nangles);
  Convert_mxArray(prhs[7], rho);
  // Set parameters to MT
  MT MT;
  MT.radius = radius[0];
  MT.lambda = lambda[0];
  MT.nre_med = nre_med[0];
  MT.nim_med = nim_med[0];
  MT.nre_p = nre_p[0];
  MT.nim_p = nim_p[0]; // [AL]
  MT.nangles = nangles[0]; // [AL]
  MT.rho = rho[0];
  // Initialize
  try {
    //MT.ErrorChecks();
    MT.Init();
  } catch(mcerror e) {
    std::string message = "Error in initializing MT: " + std::string(errorstring(e)) + "\n"; 
    mexErrMsgTxt(message.c_str());
    return;
  }
  mexPrintf("Computing... \n");
  MT.MieTheory();
  mexPrintf("...done\n");
  printf("\n"); fflush(stdout);
  // Copy solution from MT to output
  plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(plhs[0])= MT.mus_n;
  Array<double> s11_or, s11_oi, s12_or, s12_oi, s33_or, s33_oi, s43_or, s43_oi;
  //Array<double> dbsolr, dbsoli; // [AL]
  
  Convert_mxArray(&plhs[1], s11_or, s11_oi, MT.s11_nr.Nx, MT.s11_nr.Ny);
  Convert_mxArray(&plhs[2], s12_or, s12_oi, MT.s12_nr.Nx, MT.s12_nr.Ny);
  Convert_mxArray(&plhs[3], s33_or, s33_oi, MT.s33_nr.Nx, MT.s33_nr.Ny);
  Convert_mxArray(&plhs[4], s43_or, s43_oi, MT.s43_nr.Nx, MT.s43_nr.Ny);

  long ii;
  for(ii = 0; ii < MT.s11_nr.N; ii++){
    s11_or[ii] = MT.s11_nr[ii];
    s11_oi[ii] = MT.s11_ni[ii];
  }
  for(ii = 0; ii < MT.s12_nr.N; ii++){
    s12_or[ii] = MT.s12_nr[ii];
    s12_oi[ii] = MT.s12_ni[ii];
  }
  for(ii = 0; ii < MT.s33_nr.N; ii++){
    s33_or[ii] = MT.s33_nr[ii];
    s33_oi[ii] = MT.s33_ni[ii];
  }
  for(ii = 0; ii < MT.s43_nr.N; ii++){
    s43_or[ii] = MT.s43_nr[ii];
    s43_oi[ii] = MT.s43_ni[ii];
  }

  mexPrintf("Done\n");
}
