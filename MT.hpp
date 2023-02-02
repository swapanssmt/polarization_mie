#ifndef __MT_HPP__
#define __MT_HPP__
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <inttypes.h>
#include <vector>
#include "Array.hpp"
#include "arrays.hpp"
#include "complex.hpp"
#include "mie.hpp"
#include "arrays.cpp"
#include "complex.cpp"
#include "mie.cpp"
#include "Errors.hpp"
#ifndef INT_FAST64_MAX
#define INT_FAST64_MAX __INT_FAST64_MAX__
#endif

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

double eps = std::numeric_limits<double>::epsilon();

// Class for Mie Theory
class MT
{
public:
  // Constructor, destructor & Assingment operator
  MT();
  ~MT();
  MT &operator=(const MT &ref);

  // Initializes the MT after all the problem definition parameters have been given
  // Ie. constructs missing parameters
  void Init();
  // Perform MieTheory computation
  void MieTheory();
public:
  double radius, lambda, nre_med, nim_med, nre_p, nim_p, rho;
  // Number of angular discritization
  int_fast64_t nangles;
  // Calculatable parameters
  double mus_n;
  Array<double> s11_nr, s11_ni;     // elements of muller matrix (real & imaginary)
  Array<double> s12_nr, s12_ni;     // elements of muller matrix (real & imaginary)
  Array<double> s33_nr, s33_ni;     // elements of muller matrix (real & imaginary)
  Array<double> s43_nr, s43_ni;     // elements of muller matrix (real & imaginary)
};

// Constuctor, set some default values for Monte Carlo
MT::MT()
{
  radius=0.0;
  lambda=0.0;
  nre_med=0.0;
  nim_med=0.0;
  nre_p=0.0;
  nim_p=0.0;
  rho=0.0;
  nangles=1;
}

// Nothing need to be done, Arrays will kill themselves when it's time
MT::~MT()
{
}

// Assingment operator:
//  This will copy references to geometrical and parametric variables
//  Only new variables in the left hand side will be ER/EI, EBR/EBI, VER/VEI
MT &MT::operator=(const MT &ref)
{
  if (this != &ref)
  {
    radius = ref.radius;
    lambda = ref.lambda;
    nre_med = ref.nre_med;
    nim_med = ref.nim_med;
    nre_p = ref.nre_p;
    nim_p = ref.nim_p;
    nangles = ref.nangles;
    rho = ref.rho;
    s11_nr.resize(ref.s11_nr.N);
    s11_ni.resize(ref.s11_ni.N);
    s12_nr.resize(ref.s12_nr.N);
    s12_ni.resize(ref.s12_ni.N);
    s33_nr.resize(ref.s33_nr.N);
    s33_ni.resize(ref.s33_ni.N);
    s43_nr.resize(ref.s43_nr.N);
    s43_ni.resize(ref.s43_ni.N);
    long ii;

    for (ii = 0; ii < s11_nr.N; ii++)
      s11_nr[ii] = s11_ni[ii] = 0.0;
    for (ii = 0; ii < s12_nr.N; ii++)
      s12_nr[ii] = s12_ni[ii] = 0.0;
    for (ii = 0; ii < s33_nr.N; ii++)
      s33_nr[ii] = s33_ni[ii] = 0.0;
    for (ii = 0; ii < s12_nr.N; ii++)
      s43_nr[ii] = s43_ni[ii] = 0.0;
    
    mus_n = ref.mus_n;
  }
  return (*this);
}
// Initialize Monte Carlo after geometry & material parameters have been assigned
// Under MPI also communicates relevant parameters to other computers and initializes
// mersenne twister with consequetive seed numbers
void MT::Init()
{
  // Reserve memory for output variables
  int_fast64_t ii;
  s11_nr.resize((nangles+1));
  s11_ni.resize((nangles+1));
  for (ii = 0; ii <= nangles; ii++)
    s11_nr[ii] = s11_ni[ii] = 0.0;
  s12_nr.resize((nangles+1));
  s12_ni.resize((nangles+1));
  for (ii = 0; ii <= nangles; ii++)
    s12_nr[ii] = s12_ni[ii] = 0.0;
  s33_nr.resize((nangles+1));
  s33_ni.resize((nangles+1));
  for (ii = 0; ii <= nangles; ii++)
    s33_nr[ii] = s33_ni[ii] = 0.0;
  s43_nr.resize((nangles+1));
  s43_ni.resize((nangles+1));
  for (ii = 0; ii <= nangles; ii++)
    s43_nr[ii] = s43_ni[ii] = 0.0;
  return;
}
// Run Monte Carlo
void MT::MieTheory()
{
  // Single thread implementation
  double pi = 3.1415926535897932384;

	/* Mie theory stuff */
	double A;
	long i;
	struct complex m;
	struct complex *s1 = NULL;
	struct complex *s2 = NULL;

	double *mu = NULL;
	double x, qext, qsca, qback, g, vol;
	/* other variables */
	double mus;	   /* scattering coefficient [cm^-1] */
	double musp;   /* reduced scattering coefficient [cm^-1] */
	/**** allocate matrices and arrays *******/
	double *s11 = NULL;
	double *s12 = NULL;
	double *s33 = NULL;
	double *s43 = NULL;
	/**** end  allocate matrices and arrays *******/

	/* CHOOSE MIE SCATTERING parameters */
	/* Setup MIE SCATTERING parameters */

	mu = new_darray(nangles);
	s1 = new_carray(nangles);
	s2 = new_carray(nangles);
	s11 = new_darray(nangles);
	s12 = new_darray(nangles);
	s33 = new_darray(nangles);
	s43 = new_darray(nangles);

	m.re = nre_p / nre_med;
	m.im = 0.0;
	x = 2 * pi * radius / (lambda / nre_med);
	vol = 4.0 / 3 * pi * radius * radius * radius;
	A = pi * radius * radius;

	for (i = 0; i <= nangles; i++)
		mu[i] = cos(pi * i / nangles);
	s11 = new_darray(nangles);
	s12 = new_darray(nangles);
	s33 = new_darray(nangles);
	s43 = new_darray(nangles);
	s1 = new_carray(nangles);
	s2 = new_carray(nangles);

	Mie(x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g); /* <---- Call Mie program ----- */

	mus = qsca * A * rho * 1e4; /* Mus is in cm^-1 */
	musp = mus * (1 - g);		/* [cm^-1] */
	free_darray(mu);
	/*Scattering parameters s11 s12 s33 s43*/
	for (i = 0; i <= nangles; ++i)
	{
		s11[i] = 0.5 * cabbs(s2[i]) * cabbs(s2[i]) + 0.5 * cabbs(s1[i]) * cabbs(s1[i]);
		s12[i] = 0.5 * cabbs(s2[i]) * cabbs(s2[i]) - 0.5 * cabbs(s1[i]) * cabbs(s1[i]);
		s33[i] = (cmul(conj(s1[i]), s2[i])).re;
		s43[i] = (cmul(conj(s1[i]), s2[i])).im;
	}
  mus_n=mus;
  for (i = 0; i <= nangles; ++i)
	{
		s11_nr[i] = s11[i];
		s12_nr[i] = s12[i];
		s33_nr[i] = s33[i];
		s43_nr[i] = s43[i];

	}
}

#endif
