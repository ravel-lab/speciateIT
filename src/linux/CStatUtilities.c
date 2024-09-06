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

#include "CStatUtilities.h"
#include "CUtilities.h"

//---------------------------------------------------------- bsearchClosest ----
// search array, a, of doubles for given x
// returns index of the closest element
int bsearchDbl ( double *a, int n, double x )
{
  int l = 0;
  int u = n-1;
  int i;
  while (l <= u)
  {
    i = (int)((l + u)/2);

    if (a[i] < x)
    {
      l = i+1;
    }
    else if (a[i] > x)
    {
      u = i-1;
    }
    else
    {
      return i; // found x
    }
  }

  if ( fabs( a[l] - x) > fabs( a[u] - x ) )
  {
    return u;
  }
  else
  {
    return l;
  }
}

//---------------------------------------------------------- min ----
//! minimum of data array of length dataLen
double min( double * x, int n )
{
    double minEl = x[0];
    int i;

    for ( i = 0; i != n; ++i )
        if ( minEl > x[i] )
            minEl = x[i];

    return minEl;
}

//---------------------------------------------------------- max ----
//! maximum of data array of length dataLen
double max( double * x, int n )
{
    double maxEl = x[0];
    int i;

    for ( i = 0; i != n; ++i )
        if ( maxEl < x[i] )
            maxEl = x[i];

    return maxEl;
}

//----------------------------------------------------- which_max ----
//! index of a maximal element of an array of length n
int which_max( double * x, int n )
{
    double maxEl = x[0];
    int maxIdx = 0;
    int i;

    for ( i = 0; i != n; ++i )
    {
      if ( maxEl < x[i] )
	{
	  maxEl = x[i];
	  maxIdx = i;
	}
    }

    return maxIdx;
}


//----------------------------------------------------- which_max2 ----
//! index of two largest elements of an array of length n
//! i1 is the index of the largest element and i2 the second largest
void which_max2( double * x, int n, int *i1, int *i2 )
{
  double y[n];
  int i;
  for ( i = 0; i < n; ++i )
    y[i] = x[i];

  *i1 = which_max(y, n);
  int j = which_min(y, n);
  y[*i1] = y[j] - 1;
  *i2 = which_max(y, n);

  #if 0
  //qsort( x, n, sizeof(double), compareDouble2 );
  //fprintf(stderr, "in which_max2() x after sorting: ");
  fprintf(stderr, "in which_max2() x: ");
  for ( i = 0; i < n; i++)
    fprintf(stderr, "%.2f ", x[i]);
  fprintf(stderr, "\n");

  fprintf(stderr, "i1=%d\ti2=%d\n",*i1, *i2);
  #endif
}




//----------------------------------------------------- which_min ----
//! index of a minimal element of an array of length n
int which_min( double * x, int n )
{
    double minEl = x[0];
    int minIdx = 0;
    int i;

    for ( i = 0; i != n; ++i )
    {
      if ( minEl > x[i] )
	{
	  minEl = x[i];
	  minIdx = i;
	}
    }

    return minIdx;
}

//---------------------------------------------------------- medianInPlace ----
//! median of double array of size n; done in place (array will be sorted)
double medianInPlace( double *a, int n )
{
    double median = 0;

    if ( n == 1 )
        return a[0];
    else if ( n == 0 )
    {
        fprintf(stderr,"ERROR: median cannot be computed for array of size 0!\n");
        exit(EXIT_FAILURE);
    }

    qsort( a, n, sizeof(double), compareDouble );

    if( n % 2 == 0 )
        median = ( a[n/2-1] + a[n/2] ) / 2.0;
    else
        median = a[n/2];

    return median;
}


//---------------------------------------------------------- mean ----
double mean( double *a, int n )
{
  double m = 0;
  int i;

  for ( i = 0; i < n; ++i )
    m += a[i];

  return m/n;
}

//---------------------------------------------------------- TukeyBiweight ----
// Tukey biweight location estimate
double TukeyBiweight( double * x, int len, double c)
{
    double median = medianInPlace(x, len);

    // compute MAD
    double * res = (double *)malloc(len * sizeof(double));
    int i;

    for( i = 0; i < len; i++)
        res[i] = fabs(x[i] - median);

    double MAD = medianInPlace(res, len);

    free(res);

    // compute the biweighted average

    double epsilon = 0.0001; // constant to prevent the denominator to be 0
    double weightedsum = 0;
    double sum = 0;

    for( i = 0; i != len; i++)
    {
      double y = (x[i] - median)/(c * MAD+epsilon);
      double w = 0;

      if ( fabs(y) <= 1 )
          w = (1-y*y)*(1-y*y);

      weightedsum += w*x[i];
      sum += w;
    }

    return weightedsum / sum;
}



struct indexedData_t
{
    double signal;
    int index;
};

//---------------------------------------------------------- compareIndexedData ----
int compareIndexedData( const void * e1, const void * e2 )
{
    double v1 = ( ( struct indexedData_t * ) e1 )->signal;
    double v2 = ( ( struct indexedData_t * ) e2 )->signal;

    return ( v1 < v2 ) ? -1
         : ( v1 > v2 ) ?  1
         :                0;
}

//---------------------------------------------------------- quantileNormalize ----
// quantile normalization
//
// m - matrix storing signal data in columns m[i][j] element of i-th column and j-th row
// ncols = number of columns of m
// nrows = number of rows of m
double ** quantileNormalize( double ** m, int nrows, int ncols )
{
    // copy m into a matrix, t, of indexedData_t, so that after sorting its columns we can
    // get back to the original order of the elements
    struct indexedData_t ** t = (struct indexedData_t **)malloc(ncols * sizeof(struct indexedData_t *));
    CHECK_PTR ( t );

    // output matrix
    double ** out = (double **)malloc(ncols * sizeof(double *));
    CHECK_PTR ( out );
    int col, row;

    for ( col = 0; col != ncols; col++ )
    {
        t[col] = (struct indexedData_t *)malloc(nrows * sizeof(struct indexedData_t));
        CHECK_PTR ( t[col] );
        out[col] = (double *)malloc(nrows * sizeof(double));
        CHECK_PTR ( out[col] );

        for ( row = 0; row != nrows; row++ )
        {
            t[col][row].signal = m[col][row];
            t[col][row].index = row;
        }

        // sort each column of t w/r signal field
        qsort ( t[col], nrows, sizeof ( struct indexedData_t ), compareIndexedData );

    }

    // average elements of each row
    // and repopulate t[.][row] with the average value
    for ( row = 0; row != nrows; row++ )
    {
        double avg = 0;

        for ( col = 0; col < ncols; col++ )
            avg += t[col][row].signal;

        avg /= ncols;

        for ( col = 0; col < ncols; col++ )
            t[col][row].signal = avg;
    }

    // populate out matrix
    for ( col = 0; col != ncols; col++ )
    {
        for ( row = 0; row < nrows; row++ )
        {
            out[col][t[col][row].index] = t[col][row].signal;
        }
    }

    // free memory
    for ( col = 0; col != ncols; col++ )
        free(t[col]);
    free(t);

    return out;
}

//----------------------------------------------------------------- approxfun ----
//! simplification of approx1 function from R (assuming linear approximation)
//! 'x' is assumed to be sorted
double approxfun(double v, double *x, double *y, int n)
{
    /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
    int i, j, ij;

    if(!n)
    {
        fprintf(stderr,"ERROR: In approxfun():  n has to be >0: %s at line %d\n",__FILE__,__LINE__);
        exit( EXIT_FAILURE );
    }

    if (!x || !y)
    {
      return 0;
    }

    i = 0;
    j = n - 1;

    /* handle out-of-domain points */

    if(v < x[i]) return y[i];
    if(v > x[j]) return y[j];

    /* find the correct interval by bisection */

    while(i < j - 1) { /* x[i] <= v <= x[j] */
	  ij = (i + j)/2; /* i+1 <= ij <= j-1 */
	  if(v < x[ij]) j = ij;
	  else i = ij;
	/* still i < j */
    }
    /* provably have i == j-1 */

    /* interpolation */

    if(v == x[j]) return y[j];
    if(v == x[i]) return y[i];
    /* impossible: if(x[j] == x[i]) return y[i]; */

    return y[i] + (y[j] - y[i]) * ((v - x[i])/(x[j] - x[i]));
}

//----------------------------------------------------------------- ran0 ----
//! "Minimal" random number generator lifted from "Numerical Recipes" by W.H. Press at el. (p. 278)

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float ran0( long * idum )
{
    long k;

    float ans;

    *idum ^= MASK;
    k = (*idum)/IQ;
    *idum = IA * (*idum-k*IQ)-IR*k;

    if ( *idum < 0 ) *idum += IM;
    ans = AM *(*idum);
    *idum ^= MASK;

    return ans;
}


//----------------------------------------------------------------- ran1 ----
//! Another "Minimal" random number generator lifted from "Numerical Recipes" by W.H. Press at el. (p. 280)

#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1( long * idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    double temp;

    if ( *idum <= 0 || !iy )
    {
        if ( -(*idum) < 1 )
            *idum = 1;
        else
            *idum = -(*idum);

        for ( j = NTAB+7; j >= 0; --j )
        {
            k = (*idum)/IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;

            if ( *idum < 0 )
                *idum += IM;

            if ( j < NTAB )
                iv[j] = *idum;
        }
        iy = iv[0];
    }

    k = (*idum)/IQ;
    *idum = IA * (*idum-k*IQ) - IR * k;

    if ( *idum < 0 )
        *idum += IM;

    j = iy / NDIV;
    iy = iv[j];
    iv[j] = *idum;

    if ( (temp = AM * iy) > RNMX )
        return RNMX;
    else
        return temp;
}


//-------------------------------------------------------------- maxAbsDer ----
//! computes median absolute derivative within selected range
double maxAbsDer( const double *sllr, int start, int end )
{
    int n = end - start;
    int i, j;
    double absDer;
    double maxAbsDer = -1;

    for ( j = 0, i = start; j < n; ++i, ++j )
    {
        if ( (absDer = fabs(sllr[i+1] - sllr[i])) > maxAbsDer )
             maxAbsDer = absDer;
    }

    return maxAbsDer;
}
