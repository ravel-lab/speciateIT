/*
Copyright (C) 2016 Pawel Gajer (pgajer@gmail.com) and Jacques Ravel jravel@som.umaryland.edu

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

#include <iostream>
#include "StatUtilities.hh"
#include "CStatUtilities.h"
#include "CppStatUtilities.hh"
#include "CppUtilities.hh"
//#include "emd.h"

using namespace std;

//----------------------------------------------------------- cov ----
// covariance of two vectors of length n
double cov( double *x, double *y, int n )
{
  double z = 0;

  double mx = mean(x,n);
  double my = mean(y,n);

  for ( int i = 0; i < n; ++i )
    z += (x[i] - mx) * (y[i] - my);

  return z/(n-1);
}

//----------------------------------------------------------- rcov ----
// robust covariance of two vectors of length n
double rcov( double *x, double *y, int n )
{
  double z = 0;

  double mx = median(x,n);
  double my = median(y,n);

  for ( int i = 0; i < n; ++i )
    z += (x[i] - mx) * (y[i] - my);

  return z/(n-1);
}


//--------------------------------------------------------- scaledL1dist ----
// scaled L1 distance between two vectors
double scaledL1dist( double *v1, double *v2, int n)
{
  double delta = 0;

  for ( int i = 0; i < n; ++i )
    delta += fabs(v1[i]-v2[i]);

  return 1000*delta/16;
}


//--------------------------------------------------------- doubleCol ----
double * doubleCol( double **x, int nRows, int colIdx )
{
  double *col = (double*)malloc(nRows * sizeof(double));

  for ( int i = 0; i < nRows; ++i )
    col[i] = x[i][colIdx];

  return col;
}
