/*
Copyright (C) 2015 Pawel Gajer, Adam M Phillippy and Jacques Ravel jravel@som.umaryland.edu

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

////////////////////////////////////////////////////////////////////////////////
//! \file
//! \author Adam M Phillippy
//! \date 01/02/06
//!
//! \brief Utility functions for C++
//!
////////////////////////////////////////////////////////////////////////////////

#include "CppUtilities.hh"
#include <new>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <sstream>

//--------------------------------------------------------- join ----
string join( const char *tok, vector<int> &v )
{
  int n = v.size() - 1;

  ostringstream s;
  int i;
  for ( i = 0; i < n; ++i )
    s << v[i] << tok;
  s << v[i];

  return s.str();
}

//---------------------------------------------------------------- ReadLine ----
bool ReadLine (FILE * file, char * & buff, size_t & size)
{
  buff [size-1] = '$';
  if ( ! fgets (buff, size, file) )
    return false;

  while ( buff [size-1] == '\0' && buff [size-2] != '\n' )
    {
      long int old = size;
      size *= 2;
      buff = (char *) SafeRealloc (buff, sizeof (char) * size);

      buff [size-1] = '$';
      if ( ! fgets (buff+old-1, size-old+1, file) )
        return false;
    }

  return true;
}


//---------------------------------------------------------------- ReadLine ----
bool ReadLine (std::istream & file, char * & buff, size_t & size)
{
  buff [size-1] = '$';
  if ( ! file.read (buff, size) )
    return false;

  while ( buff [size-1] == '\0' && buff [size-2] != '\n' )
    {
      long int old = size;
      size *= 2;
      buff = (char *) SafeRealloc (buff, sizeof (char) * size);

      buff [size-1] = '$';
      if ( ! file.read (buff+old-1, size-old+1) )
        return false;
    }

  return true;
}


//-------------------------------------------------------------- SafeCalloc ----
void * SafeCalloc (size_t num, size_t size)
{
  void * Q = calloc (num, size);
  if ( Q == NULL )
    throw std::bad_alloc ( );
  return Q;
}


//-------------------------------------------------------------- SafeMalloc ----
void * SafeMalloc (size_t size)
{
  void * Q = malloc (size);
  if ( Q == NULL )
    throw std::bad_alloc ( );
  return Q;
}


//------------------------------------------------------------- SafeRealloc ----
void * SafeRealloc (void * P, size_t size)
{
  void * Q = realloc (P, size);
  if ( Q == NULL  &&  size != 0 )
    throw std::bad_alloc ( );
  return Q;
}


//-------------------------------------------------------------- SafeStrdup ----
char * SafeStrdup (const char * str)
{
  char * Q = strdup (str);
  if ( Q == NULL )
    throw std::bad_alloc ( );
  return Q;
}


//----------------------------------------------------------------- baseFileName ----
//! extracts base name from a file path
string baseFileName( const char *file, bool withSuffix )
{
    string fileStr(file);

    return baseFileName(fileStr, withSuffix);
}

//----------------------------------------------------------------- baseFileName ----
//! extracts base name from a file path
string baseFileName( const string & file, bool withSuffix )
{
    string base = file;
    string::size_type ind1 = file.find_last_of('/');


    if ( withSuffix && ind1 != string::npos )
    {
        base = file.substr(ind1+1);
    }
    else
    {
        string::size_type ind2 = file.find_last_of('.');


        if ( ind1 != string::npos && ind2 != string::npos )
        {
            base = file.substr(ind1+1,ind2-ind1-1);
        }
        else if ( ind1 == string::npos && ind2 != string::npos ) // no '/'
        {
            base = file.substr(0,ind2);
        }
        else if ( ind1 != string::npos && ind2 == string::npos )
        {
            base = file.substr(ind1+1);
        }

	//cerr << "file=" << file.c_str() << "\tind1=" << ind1 << "\tind2=" << ind2 << "\tbase=" << base << endl;
    }

    return base;
}

//-------------------------------------------------- pathWithoutSuffix ----
//! extracts files path without file's suffix
string pathWithoutSuffix( const string & file )
{
  string base = file;
  string::size_type ind = file.find_last_of('.');

  if ( ind != string::npos )
    base = file.substr(0,ind);

  return base;
}


//----------------------------------------------------------------- green2redRGB ----
//! interpolates between green and red color; color=red when signal<= min, color=green when signal>=max
//! color=yellow at signal=median(signal)
void green2redRGB( unsigned char& r, unsigned char& g, unsigned char& b,
                   double signal, double min, double max, double med )
{
    double t;

    if ( signal <= min )
    {
        r = 255;
        g = 0;
        b = 0;
    }
    else if ( signal >= max )
    {
        r = 0;
        g = 255;
        b = 0;
    }
    else if ( signal > min && signal <= med )
    {
        t = 1.0 / ( med - min ) * ( signal - min );
        t = t * t;
        r = 255;
        g = (unsigned char)( t * 255 );
        b = 0;
    }
    else if ( signal > med && signal <= max )
    {
        t = -1.0 / ( max - med ) * ( signal - med ) + 1.0;
        t = t * t;
        r = (unsigned char)( t * 255 );
        g = 255;
        b = 0;
    }


}

//----------------------------------------------------------------- green2redHSV ----
//! interpolates between green and red color; color=red when signal<= min, color=green when signal>=max
//! color=yellow at signal=median(signal)
void green2redHSV( unsigned char& h, unsigned char& s, unsigned char& v,
                   double signal, double min, double max, double med )
{
    double t;

    if ( signal <= min )
    {
        h = 0;
        s = 255;
        v = 255;
    }
    else if ( signal >= max )
    {
        h = 120;
        s = 255;
        v = 255;
    }
    else if ( signal > min && signal <= med )
    {
        t = 1.0/ ( med - min ) * ( signal - min );
        t = 0.5 * t;
        h = (unsigned char)( t * 120 );
        s = 255;
        v = 255;
    }
    else if ( signal > med && signal <= max )
    {
        t = 1.0 / ( max - med ) * ( signal - med );
        t = 0.5 * t + 0.5;
        h = (unsigned char)( t * 120 );
        s = 255;
        v = 255;
    }
}
