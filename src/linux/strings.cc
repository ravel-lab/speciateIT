/*
Copyright (C) 2016 Pawel Gajer (pgajer@gmail.com), Adam Phillippy and Jacques Ravel jravel@som.umaryland.edu

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

#include "strings.hh"
#include <ctype.h>

//---------------------------------------------------------------- strclean ----
/// Removes funky characters and converts string to lowercase. "Funky"
/// characters include everything but alphanums, '-', and spaces.
char * strclean(char *s)
{
  char *c = s;
  char *p = s;
  for ( ; *c; ++c )
    if ( isalnum(*c) || *c == '-' || isspace(*c) )
      *(p++) = tolower(*c);
  *p = '\0';
  return s;
}

//------------------------------------------------------------------ strpre ----
/// Returns NULL if pre is not a prefix of str, otherwise returns the
/// first position in str after the prefix pre.
char * strpre(const char* str, const char* pre)
{
  while ( 1 )
    {
      if ( *pre == '\0' ) return (char*)str;
      if ( *pre != *str ) return NULL;
      ++str; ++pre;
    }
}


//-------------------------------------------------------------- strcasepre ----
/// Same as strpre but case insensitive.
char * strcasepre(const char* str, const char* pre)
{
  while ( 1 )
    {
      if ( *pre == '\0' ) return (char*)str;
      if ( tolower(*pre) != tolower(*str) ) return NULL;
      ++str; ++pre;
    }
}


//---------------------------------------------------------------- strchomp ----
/// Removes consecutive newlines from the end of a string s.
char* strchomp(char *s)
{
  char *last = s + strlen(s);
  while ( last != s && *(--last) == '\n' )
    {}
  *(last+1) = '\0';
  return s;
}


//------------------------------------------------------------------ strlop ----
/// Removes everything after the first ch found in s.
char* strlop(char *s, char ch)
{
  char *last = s;
  while ( *last && *last != ch )
    ++last;
  *last = '\0';
  return s;
}


//------------------------------------------------------------------ strlop ----
/// Removes everything after the first whitespace found in s.
char* strlopspace(char *s)
{
  char *last = s;
  while ( *last && !isspace(*last) )
    ++last;
  *last = '\0';
  return s;
}


//---------------------------------------------------------------- strsplit ----
/// Inserts terminating character at the first space in s. Returns
/// pointer to first non space character following new terminating
/// character. Can be placed in a while loop to jump through all words
/// of s. Similar to strtok.
char* strsplit(char *s)
{
  if ( *s == '\0' ) return NULL;
  while ( *s && !isspace(*s) ) ++s;
  if ( *s )
    {
      *(s++) = '\0';
      while ( *s && isspace(*s) ) ++s;
    }
  return s;
}


//---------------------------------------------------------------- strsplit ----
/// Inserts terminating character at the first occurance of ch in s. Returns
/// pointer to first non ch character following new terminating
/// character. Can be placed in a while loop to jump through all words
/// of s. Similar to strtok.
char* strsplit(char *s, char ch)
{
  if ( *s == '\0' ) return NULL;
  while ( *s && *s != ch ) ++ s;
  if ( *s )
    {
      *(s++) = '\0';
      while ( *s && *s == ch ) ++s;
    }
  return s;
}


//----------------------------------------------------------------- strtrim ----
/// Trims leading and trailing whitespace.
char* strtrim(char *s)
{
  char *last = s + strlen(s) - 1;
  while ( last >= s && isspace(*last) )
    --last;
  *(last+1) = '\0';
  while ( *s && isspace(*s) )
    ++s;
  return s;
}


//--------------------------------------------------------------- strtrimto ----
/// Trims to the first and last occurrence of ch
char* strtrimto(char *s, char ch)
{
  char *last = strrchr(s, ch);
  if ( !last ) last = s;
  *last = '\0';

  char *first = strchr(s, ch);
  if ( !first ) first = last-1;
  return first+1;
}
