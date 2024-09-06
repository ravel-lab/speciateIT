#ifndef CSTATUTILITIES_H
#define CSTATUTILITIES_H

/*
Copyright (C) 2016 Pawel Gajer pgajer@gmail.com and Jacques Ravel jravel@som.umaryland.edu

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifdef __cplusplus
extern "C"
{
#endif

double min( double * x, int n );
double max( double * x, int n );

int bsearchDbl ( double *a, int n, double x );

int which_max( double * x, int n );
void which_max2( double * x, int n, int *i1, int *i2 );

int which_min( double * x, int n );

double mean( double *a, int n );
double medianInPlace( double *a, int n );

double ** quantileNormalize( double ** m, int nrows, int ncols);

double ran1( long * idum);

double approxfun(double v, double *x, double *y, int n);

double maxAbsDer( const double *sllr, int start, int end );

double TukeyBiweight( double * x, int len, double c);

#ifdef __cplusplus
}
#endif

#endif
