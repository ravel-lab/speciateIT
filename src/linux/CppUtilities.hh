/*
Copyright (C) 2016 Pawel Gajer, Adam M Phillippy and Jacques Ravel jravel@som.umaryland.edu

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

#ifndef CPPUTILITIES_HH
#define CPPUTILITIES_HH
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

//================================================================= Globals ====

/// join
string join( const char *tok, vector<int> &v );

/// Computer the log of a number in any base
double Log (double base, double val);

/// Read a line from a file into a buffer, resizing the buffer if necessary
bool ReadLine (FILE * file, char * & buff, size_t & size);

/// Read a line from a file into a buffer, resizing the buffer if necessary
bool ReadLine (std::istream & file, char * & buff, size_t & size);

/// Callocs an array, throws std::bad_alloc exception on failure
void * SafeCalloc (size_t num, size_t size);

/// Mallocs an array, throws std::bad_alloc exception on failure
void * SafeMalloc (size_t size);

/// Reallocs an array, throws std::bad_alloc exception on failure
void * SafeRealloc (void * P, size_t size);

/// Strdups a string, throws std::bad_alloc exception on failure
char * SafeStrdup (const char * str);

/// Extracts base name from a file path
string baseFileName( const string & file, bool withSuffix = false );
string baseFileName( const char *file, bool withSuffix );
string pathWithoutSuffix( const string & file );

void green2redRGB( unsigned char &r, unsigned char &g, unsigned char &b,
                   double signal, double min, double max, double med );

void green2redHSV( unsigned char &h, unsigned char &s, unsigned char &v,
                   double signal, double min, double max, double med );

//================================================================= Inlines ====

inline double Log (double base, double val)
{
  if ( base <= 0 || val <= 0 ) { return val; }
  else { return log (val) / log (base); }
}

#endif // CPPUTILITIES_HH
