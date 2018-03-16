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


#include "CUtilities.h"
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

#define MAXMESSAGE 1024
static char Message[MAXMESSAGE];

//-------------------------------------------------------------- error_msg ----
void error_msg( const char *file, int line, const char *msg )
{
  fprintf(stderr, "\n\n\tERROR: in %s line %d: %s\n\n", file, line, msg );
}


//--------------------------------------------------------------- FileCheck ----
/// Attempts to open/close file using mode. Validates the file by
/// actually performing an open/close operation. Therefore, if mode is
/// "w", file will be truncated.
void FileCheck(const char *file, const char *mode)
{
  FileClose(FileOpen(file, mode), file);
}


//--------------------------------------------------------------- FileClose ----
/// Closes file pointed to by fp. Aborts program with appropriate
/// message if operation fails.
void FileClose(FILE *fp, const char *file)
{
  if ( fclose(fp) != 0 )
    {
      if ( file == NULL )
        snprintf(Message, MAXMESSAGE, "Could not close file");
      else
        snprintf(Message, MAXMESSAGE, "Could not close %s", file);
      perror(Message);
      DIE();
    }
}


//---------------------------------------------------------------- FileOpen ----
/// Opens file using mode. Aborts program with appropriate message if
/// operation fails.
FILE* FileOpen(const char *file, const char *mode)
{
  FILE* retval = fopen(file,mode);
  if ( retval == NULL )
    {
      snprintf(Message, MAXMESSAGE, "Could not open %s", file);
      perror(Message);
      DIE();
    }
  return retval;
}


//-------------------------------------------------------------- ReadString ----
int ReadString(FILE *fp, char **buff, size_t *size)
{
  char c = 1;
  size_t cnt = 0;
  while ( c != '\0' && c != EOF )
    {
      if ( cnt >= *size )
        {
          *size *= 2;
          REALLOC(*buff, char*, *size * sizeof(char));
        }

      c = getc(fp);
      (*buff)[cnt++] = c;
    }

  (*buff)[--cnt] = '\0';

  if ( c == EOF )
    return 0;
  return 1;
}



//------------------------------------------------------------- WriteString ----
int WriteString (FILE *fp, const char *str)
{
  if ( !str ) str = "";
  return ( fputs(str,fp) != EOF && putc('\0',fp) != EOF );
}


int isDouble (const char * str)
{
  char * eptr;
  strtod (str, &eptr);
  return ( str != eptr && *eptr == '\0' );
}

char stringToChar (const char * str)
{
  for ( ; *str != '\0'; ++str )
    if ( !isspace (*str) )
      return *str;
  return ' ';
}


//------------------------------------------------------------------- chomp ----
/// Removes whitespace from the end of a string s.
char* chomp(char *s)
{
  char *last = s + strlen(s);
  while ( last != s && isspace(*(--last)) )
    {}
  *(last+1) = '\0';
  return s;
}

//----------------------------------------------------------------- fileSize ----
//! returns file size
// see stat_struct.h for members of stat struct (taken from bits/stat.h)
size_t fileSize ( const char * fileName )
{
    struct stat fileInfo;

    if ( stat( fileName, &fileInfo ) != 0 ) // Use stat( ) to get the info
    {
        fprintf(stderr,"Cannot read file info for %s: %s\n", fileName, strerror ( errno ) );
        exit( EXIT_FAILURE );
    }

    return fileInfo.st_size;
}

//----------------------------------------------------------------- readLine ----
// reads a line from a buffer of size bSize, starting at a specified offset
// of the buffer and writes the line into an array s of size sSize
// the new line character is not written into s
// getline returns the number of characters written into s
int readLine ( int offset, const char * buffer, int bSize, char * s, int sSize )
{
    int i,j = 0;

    if ( offset >= bSize )
    {
        fprintf ( stderr, "offset:%d >= bSize:%d\n",offset,bSize);
        return 0;
    }

    for ( j = 0, i = offset; i != bSize && j != sSize && buffer[i] != '\n' ; i++, j++ )
        s[j] = buffer[i];

    s[j] = '\0';

    return j+1;
}

//----------------------------------------------------------------- compareInt ----
int compareInt ( void const *e1, void const *e2 )
{
    int v1 = *( ( int * ) e1 );
    int v2 = *( ( int * ) e2 );

    return ( v1 < v2 ) ? -1 : ( v1 > v2 ) ? 1 : 0;
}

//----------------------------------------------------------------- compareDouble ----
int compareDouble ( void const *e1, void const *e2 )
{
    double v1 = *( ( double * ) e1 );
    double v2 = *( ( double * ) e2 );

    return ( 0 < v2 - v1) ? -1 : ( v1 - v2 > 0 ) ? 1 : 0;
}

//----------------------------------------------------------------- compareDouble ----
int compareDouble2 ( void const *e1, void const *e2 )
{
    double v1 = *( ( double * ) e1 );
    double v2 = *( ( double * ) e2 );

    return ( 0 < v2 - v1) ? 1 : ( v1 - v2 > 0 ) ? -1 : 0;
}



//----------------------------------------------------------------- newmem ----
//! returns char * to new block of memory; prints memory allocation error message
void *newmem ( int number, int size )
{
    void *memp = NULL;
    int len = number * size;

    errno = 0;

    if ( len > 0 )
    {
        memp = malloc ( len );

        if ( memp == NULL )
        {
            fprintf ( stderr, "Memory allocation error: %d requested: %s.\n",
                      len, strerror ( errno ) );
            exit ( EXIT_FAILURE );
        }
    }
    else if ( len < 0 )
    {
        fprintf ( stderr, "Negative memory (%d) requested!\n", len );
        exit ( EXIT_FAILURE );
    }

    return ( memp );
}

//----------------------------------------------------------------- newcmem ----
void *newcmem ( int number, int size )
{
    void *memp = NULL;
    int len = number*size;

    errno = 0;

    if ( len > 0 )
    {
        memp = calloc ( number, size);

        if ( memp == NULL )
        {
            fprintf ( stderr, "Memory allocation error: %d requested: %s.\n",
                      len, strerror ( errno ) );
            exit ( EXIT_FAILURE );
        }
    }
    else if ( len < 0 )
    {
        fprintf ( stderr, "Negative memory (%d) requested!\n", len );
        exit ( EXIT_FAILURE );
    }

    return ( memp );
}
