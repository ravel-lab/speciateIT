#ifndef CUTILITIES_H
#define CUTILITIES_H

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

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>

struct point3d_t
{
    double x;
    double y;
    double z;
};


void error_msg( const char *file, int line, const char *msg );
#define errorMsg(x) error_msg( __FILE__, __LINE__, (x) )


/// Exits if the pointer is NULL
#define CHECK_PTR(p)                                     \
  do{                                                    \
  if ( !(p) ) {                                          \
      fprintf(stderr,"ERROR: Out of memory in file %s at line %d\n",__FILE__,__LINE__); \
      exit (EXIT_FAILURE);                               \
    }                                                    \
  } while(0);                                            \

/// Print current line and exit
#define DIE()                                                           \
  do {                                                                  \
    fprintf(stderr,"ERROR in %s at line: %d\n",__FILE__,__LINE__);  \
    exit(EXIT_FAILURE);                                                 \
  } while(0)

/// Print a message and then DIE()
#define DIEM(m) do { fprintf(stderr,m"\n"); DIE(); } while(0)

/// DIE() if the pointer is NULL
#define DIENULL(p) do { if ( !(p) ) { DIEM("Out of memory"); } } while(0)

/// Calloc memory, DIE() if NULL
#define CALLOC(p,t,s) DIENULL( (p) = (t) calloc(s,1) )

/// Realloc memory, DIE() if NULL
#define REALLOC(p,t,s) DIENULL( (p) = (t) realloc(p, s) )

/// Malloc memory, DIE() if NULL
#define MALLOC(p,t,s) DIENULL( (p) = (t) malloc(s) )

/// Malloc memory for objects of type t, DIE() if NULL
#define NEW(p,t) MALLOC(p, t*, sizeof(t))

/// Strdup and DIE() if NULL
#define STRDUP(d,s) DIENULL( (d) = strdup(s) )

/// Check for a file by opening and closing a file
void FileClose(FILE *fp, const char *file);

/// Try to close a file, or DIE()
void FileCheck(const char *file, const char *mode);

/// Try to open a file, or DIE()
FILE* FileOpen(const char *file, const char *mode);

/// Read a string to a buffer, resize buffer if necessary, DIE() on error
int ReadString(FILE *fp, char **buff, size_t *size);

/// Write a string to a file and check for sucess
int WriteString(FILE *fp, const char *str);

/// Return 1 if string represents a double
int isDouble (const char * str);

/// Return ascii value of first non-whitespace character
char stringToChar (const char * str);

/// Clean whitespace off the end of a string
char* chomp(char *s);

/// Stat a file to return its filesize, exit on failure
size_t fileSize(const char * fileName);

/// Reads a line from a buffer
int readLine(int offset, const char * buffer, int bSize, char * s, int sSize);

/// Return -1 if *v1 < *v2, 0 if =, 1 if > (for quicksort)
int compareInt(void const *e1, void const *e2);

/// sorts in an increasing order
int compareDouble(void const *e1, void const *e2);
/// sorts in a decreasing order
int compareDouble2 ( void const *e1, void const *e2 );


#ifdef __cplusplus
}
#endif

#endif
