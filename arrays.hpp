#ifndef __ARRAYS_HPP__
#define __ARRAYS_HPP__
#include <inttypes.h>
/*3:*/
#line 27 "arrays.w"

/*7:*/
#line 51 "arrays.w"

double*new_darray(long size)

/*:7*/
#line 28 "arrays.w"
;
/*9:*/
#line 70 "arrays.w"

void free_darray(double*a)

/*:9*/
#line 29 "arrays.w"
;
/*11:*/
#line 82 "arrays.w"

double*copy_darray(double*a,long size)

/*:11*/
#line 30 "arrays.w"
;
/*13:*/
#line 100 "arrays.w"

void set_darray(double*a,long size,double x)

/*:13*/
#line 31 "arrays.w"
;
/*15:*/
#line 116 "arrays.w"

void min_max_darray(double*a,long size,double*min,double*max)

/*:15*/
#line 32 "arrays.w"
;
/*18:*/
#line 144 "arrays.w"

void sort_darray(double*a,long size)

/*:18*/
#line 33 "arrays.w"
;
/*21:*/
#line 188 "arrays.w"

void print_darray(double*a,long size,long ilow,long ihigh)

/*:21*/
#line 34 "arrays.w"
;

/*:3*/
#endif