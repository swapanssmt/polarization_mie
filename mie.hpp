#ifndef __MIE_HPP__
#define __MIE_HPP__
#include <inttypes.h>

struct complex Lentz_Dn(struct complex z,long n);
void Dn_down(struct complex z,long nstop,struct complex*D);
void Dn_up(struct complex z,long nstop,struct complex*D);
void small_Mie(double x,struct complex m,double*mu,
long nangles,struct complex*s1,
struct complex*s2,double*qext,double*qsca,
double*qback,double*g);
void Mie(double x,struct complex m,double*mu,long nangles,struct complex*s1,
struct complex*s2,double*qext,double*qsca,double*qback,double*g);
void ez_Mie(double x,double n,double*qsca,double*g);
#endif