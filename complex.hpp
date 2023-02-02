#ifndef __COMPLEX_HPP__
#define __COMPLEX_HPP__
#include <inttypes.h>
struct complex{double re,im;};
struct complex cset(double a,double b);
struct complex cpolarset(double r,double theta);
double cabbs(struct complex z);
double carg(struct complex z);
struct complex csqr(struct complex z);
struct complex conj(struct complex z);
double cnorm(struct complex z);
struct complex csqrt(struct complex z);
struct complex cinv(struct complex w);
struct complex cadd(struct complex z,struct complex w);
struct complex csub(struct complex z,struct complex w);
struct complex cmul(struct complex z,struct complex w);
struct complex cdiv(struct complex z,struct complex w);
double crdiv(struct complex z,struct complex w);
double crmul(struct complex z,struct complex w);
struct complex csadd(double x,struct complex z);
struct complex csdiv(double x,struct complex w);
struct complex csmul(double x,struct complex z);
struct complex csin(struct complex z);
struct complex ccos(struct complex z);
struct complex ctan(struct complex z);
struct complex casin(struct complex z);
struct complex cacos(struct complex z);
struct complex catan(struct complex z);
struct complex csinh(struct complex z);
struct complex ccosh(struct complex z);
struct complex ctanh(struct complex z);
struct complex catanh(struct complex z);
struct complex casinh(struct complex z);
struct complex cexp(struct complex z);
struct complex clog(struct complex z);
struct complex clog10(struct complex z);
struct complex*new_carray(long size);
void free_carray(struct complex*a);
struct complex*copy_carray(struct complex*a,long size);
void set_carray(struct complex*a,long size,struct complex z);
#endif