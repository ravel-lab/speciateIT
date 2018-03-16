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


#include <iostream>
#include "CppStatUtilities.hh"
#include "CUtilities.h"

using namespace std;

//==================================================== summaryStats_t ====
//---------------------------------------------------- summaryStats_t ----
summaryStats_t::summaryStats_t()
  : min(0), max(0), med(0), mean(0), mad(0), q1(0), q3(0)
{ }

//---------------------------------------------------------- opeartor= ----
summaryStats_t & summaryStats_t::operator= (const summaryStats_t & other)
{
  if ( this != &other )
  {
    min  = other.min;
    max  = other.max;
    med  = other.med;
    mean = other.mean;
    mad  = other.mad;
    q1   = other.q1;
    q3   = other.q3;
  }

  return *this;
}

//--------------------------------------------------------------- print ----
void summaryStats_t::print() const
{
  cerr << "min=" << min << "\tq1=" << q1
       << "\tmed=" << med << "\tmean=" << mean
       << "\tq3=" << q3 << "\tmax=" << max
       << "\tmad=" << mad << endl;
}

//--------------------------------------------------------------- print ----
void summaryStats_t::print( FILE * f ) const
{
  fprintf ( f, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", min, q1, med, mean, q3, max, mad );
}

//---------------------------------------------------------- summaryStats ----
void summaryStats( double * data, int dataSize, summaryStats_t &summary, int sample, double constant )
{
  int len = sample;
  int xLen = ( dataSize < len ) ? dataSize : len;
  double * perc;
  double prob[] = {0, 1, 0.5, 0.25, 0.75 };
  double mad;

  if ( xLen == dataSize )
  {
    perc = percentiles (data, xLen, prob, 5, 100);

    double * res = (double *) malloc ( xLen * sizeof ( double) );
    CHECK_PTR ( res );

    for( int i = 0; i < xLen; ++i )
      res[i] = fabs(data[i] - perc[2]);

    mad = constant * medianInPlace ( res, xLen );

    free ( res );
  }
  else
  {
    double * x = (double *)malloc (xLen * sizeof ( double ) );
    CHECK_PTR ( x );

    double * res = (double *) malloc ( xLen * sizeof ( double) );
    CHECK_PTR ( res );

    //        srand();
    for ( int k = 0; k < xLen; k++ )
      x[k] = (double)data[ rand() % xLen ];

    perc = percentiles (x, xLen, prob, 5, 100);

    for( int i = 0; i < xLen; ++i )
      res[i] = fabs(x[i] - perc[2]);

    mad = constant * medianInPlace ( res, xLen );

    free ( res );
    free ( x );
  }

  summary.min = perc[0];
  summary.max = perc[1];
  summary.med = perc[2];
  summary.mad = mad;
  summary.q1  = perc[3];
  summary.q3  = perc[4];

  // computing mean
  double mean = 0;
  for( int i = 0; i < dataSize; ++i )
    mean += data[i];
  mean /= dataSize;

  summary.mean = mean;

  free ( perc );
}

//---------------------------------------------------------- whichMax ----
//! returns index of maximum value of data array of length dataLen
void whichMax( double * data, int dataLen, double & maxVal, int & maxInd )
{
  maxVal = data[0];
  maxInd = 0;
  int i  = 0;

  for ( ; i < dataLen; ++i )
  {
    if ( data[i] > maxVal )
    {
      maxVal = data[i];
      maxInd = i;
    }
  }
}


//---------------------------------------------------------- whichMin ----
//! returns index of minimum value of data array of length dataLen
void whichMin( double * data, int dataLen, double & minVal, int & minInd )
{
  minVal = data[0];
  minInd = 0;
  int i  = 0;

  for ( ; i < dataLen; ++i )
  {
    if ( data[i] < minVal )
    {
      minVal = data[i];
      minInd = i;
    }
  }
}

//---------------------------------------------------------- mad ----
//! MAD of a double array from start to end
double mad( double *a, int start, int end, double constant )
{
  int len = end - start + 1;
  double * res = (double *) malloc ( len * sizeof ( double) );
  CHECK_PTR ( res );

  double med = median ( a, start, end );

  int i;
  end++;
  for( i = start; i < end; ++i )
    res[i-start] = fabs(a[i] - med);

  double mad = medianInPlace ( res, len );

  free ( res );

  return constant * mad;
}

//---------------------------------------------------------- mad ----
//! MAD of a double array from start to end
double mad( double *a, int start, int end, double & medianVal, double constant )
{
  int len = end - start + 1;
  double * res = (double *) malloc ( len * sizeof ( double) );
  CHECK_PTR ( res );

  medianVal = median ( a, start, end );

  int i;
  end++;
  for( i = start; i < end; ++i )
    res[i-start] = fabs(a[i] - medianVal);

  double mad = medianInPlace ( res, len );

  free ( res );

  return constant * mad;
}

//--------------------------------------------------------- weighted mad ----
//! MAD of a double array from start to end using weights; w is assumed to have the same length as a
double mad( double *a, double *w, int start, int end, double & medianVal, double constant )
{
  int len = end - start + 1;
  double * res = (double *)malloc( len * sizeof(double) );
  CHECK_PTR(res);

  medianVal = median( a, w, start, end );

  int k = 0;
  end++;
  for( int i = start; i < end; ++i )
    if ( w[i] != 0 )
      res[k++] = fabs(w[i]*a[i] - medianVal);

  double mad = medianInPlace( res, k );

  free(res);

  return constant * mad;
}

//---------------------------------------------------------- madInPlace ----
//! MAD of a double array, a, of size n; done in place (array will be sorted)
//! user supplies array res of the same size as a
double madInPlace( double *a, int n, double *res, double constant )
{
  double med = medianInPlace( a, n );
  int i;
  double mad;

  for( i = 0; i < n; ++i )
    res[i] = fabs(a[i] - med);

  mad = medianInPlace( res, n );

  return constant * mad;
}


//---------------------------------------------------------- madInPlace ----
//! MAD of a double array, a, of size n; done in place (array will be sorted)
//! user supplies array res of the same size as a
void madInPlace( double * a, int n, double * res, double & median, double & mad, double constant )
{
  median = medianInPlace( a, n );

  int i;
  for( i = 0; i < n; ++i )
    res[i] = fabs(a[i] - median);

  mad = constant * medianInPlace( res, n );
}


//---------------------------------------------------------- median ----
//! median of double array from start to end
double median( double *a, int start, int end )
{
  double median = 0;
  int len = end - start + 1;

  if ( len == 1 )
    return a[start];
  else if ( len == 2 )
    return (a[start] + a[end])/2.0;

  double * b = (double *) malloc ( len * sizeof ( double ) );
  CHECK_PTR ( b );

  memcpy ( b, a + start, len * sizeof(double) );

  qsort ( b, len, sizeof ( double ), compareDouble );

  if( len % 2 == 0 )
    median = ( b[len/2-1] + b[len/2] ) / 2;
  else
    median = b[len/2];

  free(b);

  return median;
}

//---------------------------------------------- naive weighted median ----
//! weighted median of double array, a, from start to end;
//! weight array, w, is assumed to have the same length as 'a'
double median( double *a, double *w, int start, int end )
{
  double median = 0;
  int len = end - start + 1;

  if ( len == 1 )
    return a[start];
  else if ( len == 2 )
    return (a[start] + a[end])/2.0;

  double * b = (double*)malloc( len * sizeof(double) );
  CHECK_PTR(b);

  memcpy( b, a + start, len * sizeof(double) );

  for( int i = 0; i < len; ++i )
    b[i] *= w[start+i];

  qsort( b, len, sizeof(double), compareDouble );

  if( len % 2 == 0 )
    median = ( b[len/2-1] + b[len/2] ) / 2;
  else
    median = b[len/2];

  free(b);

  return median;
}


//---------------------------------------------------------- median ----
//! median of data array of length dataLen
double median( double * data, int dataLen )
{
  return median( data, 0, dataLen-1 );
}

//---------------------------------------------------------- percentiles ----
//! computes percentiles of array data of length dataLen
//! the precentiles are given by an array prob (its values have to be between 0 and 1)
//! of length probLen;
//! percentiles are computed for random sample of size perLen/100 * dataLen, where perLen
//! has to be greater than 0 and less than or equal to 100
//! returns array or size probLen containing percentiles corresponding to elements of prob
double * percentiles( const double *data, int dataLen, const double *prob, int probLen, int perLen )
{
  if ( perLen <= 0 || perLen > 100 )
  {
    cerr << "ERROR: " << __FILE__ << "\tat line=" << __LINE__
	 << ": perLen needs to be >= 0 and <= 100." << endl;
    exit(EXIT_FAILURE);
  }

  for ( int i = 0; i < probLen; ++i )
    if ( prob[i] < 0 || prob[i] > 1 )
    {
      cerr << "ERROR: " << __FILE__ << "\tat line=" << __LINE__
	   << ": Every prob value has to be >= 0 and <= 1." << endl;
      cerr << "prob[" << i << "]=" << prob[i] << endl;
      exit(EXIT_FAILURE);
    }

  int sampleSize;

  if ( perLen < 100 )
    sampleSize = (int)( perLen/100.0 * dataLen );
  else
    sampleSize = dataLen;

  double *b = (double *)malloc( sampleSize * sizeof(double) );
  CHECK_PTR ( b );

  if ( perLen == 100 )
  {
    for ( int k = 0; k < sampleSize; k++ )
      b[ k ] = data[ k ];
  }
  else
  {
    //        srand();
    for ( int k = 0; k < sampleSize; k++ )
      b[ k ] = (double)data[ rand() % dataLen ];
  }

  qsort(b, sampleSize, sizeof(double), compareDouble );

  double *perc = (double *)malloc( probLen * sizeof(double) );
  CHECK_PTR ( perc );

  for ( int k=0; k < probLen; k++ )
  {
    double x = (sampleSize-1) * prob[k];
    double t = x - floor(x);

    if ( (int)floor(x) == sampleSize-1 )
      perc[ k ] =  b[ sampleSize-1 ];
    else
      perc[ k ] = (1-t) * b[ (int)x ] + t * b[ (int)(x+1) ];
  }

  free(b);

  return perc;
}

//----------------------------------------------------------------- hist ----
//! compute histogram for array x of size xSize using nBins
// hx - array of bin mid points
// hy - array of relative frequencies of elements of x in each bin
void hist( double *x, int xSize, int nBins, double **_hx, double **_hy )
{
  double xmin = min( x, xSize );
  double xmax = max( x, xSize );

  double *hx = (double *)malloc( nBins * sizeof(double) );
  double *hy = (double *)calloc( nBins,  sizeof(double) );

  double dx  = (xmax - xmin)/nBins;
  double dx2 = dx/2;

  //fprintf(stderr,"xmin=%f\txmax=%f\tdx=%f\tdx2=%f\n",xmin,xmax,dx,dx2);

  for ( int i = 0; i < nBins; ++i )
    hx[i] = xmin + dx2 + i*dx;

  // computing frequencies of x elements within each bin
  for ( int i = 0; i < xSize; ++i )
    ++hy[ (int)floor( (x[i] - xmin)/dx ) ];

  // computing relative frequencies
  double f = xSize * dx;
  for ( int i = 0; i < nBins; ++i )
    hy[i] /= f;

  *_hx = hx;
  *_hy = hy;
}

//----------------------------------------------------------------- hist ----
//! compute histogram for arrays x, y of sizes xSize and ySize respoctively
//  using nBins
//  hmid - array of bin mid points
//  hx - array of relative frequencies of elements of x in each bin
//  hy - array of relative frequencies of elements of y in each bin
void hist( double *x, int xSize,
	   double *y, int ySize,
	   int nBins, double **_hmid, double **_hx, double **_hy )
{
  double xmin = min( x, xSize );
  double xmax = max( x, xSize );

  double ymin = min( y, ySize );
  double ymax = max( y, ySize );

  double tmin = (xmin < ymin) ? xmin : ymin;
  double tmax = (xmax > ymax) ? xmax : ymax;

  double *hmid = (double *)malloc( nBins * sizeof(double) );
  double *hx   = (double *)calloc( nBins, sizeof(double) );
  double *hy   = (double *)calloc( nBins,  sizeof(double) );

  double dt  = (tmax - tmin)/nBins;
  double dt2 = dt/2;

  //fprintf(stderr,"xmin=%f\txmax=%f\tdt=%f\tdt2=%f\n",xmin,xmax,dt,dt2);

  for ( int i = 0; i < nBins; ++i )
    hmid[i] = tmin + dt2 + i*dt;

  // computing frequencies of x elements within each bin
  for ( int i = 0; i < xSize; ++i )
    ++hx[ (int)floor( (x[i] - tmin)/dt ) ];

  // computing frequencies of y elements within each bin
  for ( int i = 0; i < ySize; ++i )
    ++hy[ (int)floor( (y[i] - tmin)/dt ) ];


  // computing relative frequencies
  double fx = xSize * dt;
  double fy = ySize * dt;

  for ( int i = 0; i < nBins; ++i )
  {
    hx[i] /= fx;
    hy[i] /= fy;
  }

  *_hmid = hmid;
  *_hx = hx;
  *_hy = hy;
}

//--------------------------------------------------------- signalPeakPar ----
//! returns array of peakPar_t for signal array x of size xSize
void signalPeakPar( const double * x, int xSize, PeakPar_t * & peakPar, int & nPeaks )
{
  // arrays holding indices of left local min, local max and right local min indices of xyData peaks
  int * minLInd = (int *) calloc( xSize, sizeof(int) );
  CHECK_PTR(minLInd);

  int * maxInd  = (int *) calloc( xSize, sizeof(int) );
  CHECK_PTR(maxInd);

  int * minRInd = (int *) calloc( xSize, sizeof(int) );
  CHECK_PTR(minRInd);

  // advance to first i: x[i+1] != x[i]
  int i = 0;
  if ( xSize > 1 && x[1] == x[0] )
  {
    while ( i+1 < xSize && x[i+1] == x[i] )
      ++i;
  }

  int count = 0;  // peak counter; gets advanced by 1 only when minR is found
  minLInd[count] = i;
  double minVal = x[i];
  double maxVal = x[i];

  ++i;

  while ( i < xSize )
  {
    // looking for minima
    if ( x[i] < minVal )
    {
      minVal = x[i];

    }
    else if ( x[i] > minVal && x[i-1] == minVal )
    {
      if ( i > 2 && x[i-2] > minVal )
      {
	minRInd[count] = i-1;
	++count;
      }

      minLInd[count] = i-1;
      maxVal         = x[i-1];

    }
    else if ( x[i] == minVal && x[i-1] == minVal )
    {
      minRInd[count] = i-1;
      ++count;
      maxVal = x[i-1];

      while ( i < xSize && x[i] == minVal )
      {
	++i;
      }
      i--;
    }

    // looking for maxima
    if ( x[i] > maxVal )
    {
      maxVal = x[i];
    }
    else if ( x[i] < maxVal && x[i-1] == maxVal ) // local max found
    {
      maxInd[count] = i-1;
      minVal        = x[i];
    }
    ++i;
  }
  --i;

  // The end of the range case. Now, i = xSize-1.
  if ( (minLInd[count] || count==0) && x[i] >= maxVal )
  {
    maxInd[count]  = i;
    minRInd[count] = i;
  }
  else if ( (minLInd[count] || count==0) && x[i] <= minVal )
  {
    minRInd[count] = i;
  }

  nPeaks = count + 1;
  peakPar = (PeakPar_t *) malloc( nPeaks * sizeof(PeakPar_t) );

  for ( i = 0; i < nPeaks; ++i )
  {
    peakPar[i].startVal = x[minLInd[i]];
    peakPar[i].startInd = minLInd[i];

    peakPar[i].endVal   = x[minRInd[i]];
    peakPar[i].endInd   = minRInd[i];

    peakPar[i].maxVal   = x[maxInd[i]];
    peakPar[i].maxInd   = maxInd[i];
  }

  free( minLInd );
  free( minRInd );
  free( maxInd );
}
