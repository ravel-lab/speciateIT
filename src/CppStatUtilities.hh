#ifndef CPPSTATUTILITIES_HH
#define CPPSTATUTILITIES_HH

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

#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <algorithm>    // std::sort
#include "CStatUtilities.h"

using namespace std;

const double MAD_CONST = 1.4826;

//! holds signal parameters of a signle signal peak
struct PeakPar_t
{
    int maxInd;      //!< index of the local maximum value
    double maxVal;   //!< local maximum value
    int startInd;    //!< index of the position of the beginning of the signal peak
    double startVal; //!< signal intensity at startInd
    int endInd;      //!< index of the position of the end of the signal peak
    double endVal;   //!< signal intensity at endInd
};

void signalPeakPar( const double * x, int xSize, PeakPar_t * & peakPar, int & nPeaks );


//! holds min, max, median, mean, first and third quadrants of a dataset, and MAD
struct summaryStats_t
{
    double min;
    double max;
    double med;
    double mean;
    double mad;
    double q1;
    double q3;

    summaryStats_t();
    summaryStats_t & operator= (const summaryStats_t & other);
    void print() const;
    void print( FILE * f ) const;
};

void summaryStats( double * a, int dataSize, summaryStats_t & summary, int sample = 10000, double constant = MAD_CONST);

void whichMax( double * data, int dataLen, double & maxVal, int & maxInd );
void whichMin( double * data, int dataLen, double & minVal, int & minInd );

double mad( double *a, int start, int end, double constant = MAD_CONST);
double mad( double *a, int start, int end, double & medianVal, double constant = MAD_CONST);
double madInPlace( double *a, int n, double *res, double constant = MAD_CONST);
void madInPlace( double * a, int n, double * res, double & median, double & mad, double constant = MAD_CONST);

double median( double *a, double *w, int start, int end );
double median( double *a, int start, int end );
double median( double *a, int dataLen );
double median(vector<double> v);

double * percentiles( const double * data, int dataLen, const double * prob, int probLen, int perSize );

void hist( double *x, int xSize, int nBins, double **hx, double **hy );

void hist( double *x, int xSize,
	   double *y, int ySize,
	   int nBins, double **_hmid, double **_hx, double **_hy );


//------------------------------------------------------- median ----
//! minimum of vector v
template <typename T>
T medianInPlace( vector<T> &v )
{
  int n = v.size();

  if ( !n )
    fprintf(stderr,"Error in %s:medianInPlace() at line %d: argument size=0\n",
	    __FILE__,__LINE__);

  sort(v.begin(),v.end());

  if( n % 2 == 0 )
    return ( v[n/2-1] + v[n/2] ) / 2;
  else
    return v[n/2];
}

//---------------------------------------------------------- min ----
//! minimum of vector v
template <typename T>
T min( vector<T> &v )
{
  T minEl = v[0];
  int n = v.size();

  for ( int i = 0; i < n; ++i )
    if ( minEl > v[i] )
      minEl = v[i];

  return minEl;
}

template<typename T>
T min( T *x, int n )
{
  T minEl = x[0];
  int i;

  for ( i = 0; i < n; ++i )
    if ( minEl > x[i] )
      minEl = x[i];

  return minEl;
}


//---------------------------------------------------------- max ----
//! maximum of vector v
template <typename T>
T max( vector<T> &v )
{
  T maxEl = v[0];
  int n = v.size();

  for ( int i = 0; i < n; ++i )
    if ( maxEl < v[i] )
      maxEl = v[i];

  return maxEl;
}

template<typename T>
T max( T *x, int n )
{
  T maxEl = x[0];
  int i;

  for ( i = 0; i < n; ++i )
    if ( maxEl < x[i] )
      maxEl = x[i];

  return maxEl;
}

//-------------------------------------------------------- permuteInPlace ----
/// Fisher-Yates shuffle as implemented by Durstenfeld
template<typename T>
void permuteInPlace(T *x, int n)
{
  int k;
  T tmp;

  while (n > 1)
  {
    k = random() % n; // 0 <= k < n.
    n--;              // n is now the last pertinent index;
    tmp  = x[n];      // swap x[n] with x[k] (does nothing if k == n).
    x[n] = x[k];
    x[k] = tmp;
  }
}

//------------------------------------------------- Euclidean squared dist ----
//! maximum of vector v
template <typename T>
T eucDist2( T *x, T *y, int n )
{
  T d = 0;
  T z;

  for ( int i = 0; i < n; ++i )
  {
    z = x[i] - y[i];
    d += z * z;
  }

  return d;
}

#endif
