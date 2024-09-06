/*
 * C Library of Input/Output routines
 *
 * Copyright (C) 2016 Pawel Gajer pgajer@gmail.com and Jacques Ravel jravel@som.umaryland.edu
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 *
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

#ifndef IOUTILITIES_H
#define IOUTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef __cplusplus
extern "C"
{
#endif

int exists(const char *fname);

FILE *_fOpen ( const char *file, const char *format, const char *sourceFile, int line );

#define fOpen(x,y)   _fOpen((x), (y), __FILE__, __LINE__)

char * readTable( const char *inFile, double ***matrix, int *nrow, int *ncol,
		  char ***rowNames, char ***colNames );

int readCharTbl( const char *inFile, char ****tbl, int *nRows, int *nCols);

void printCharTbl(char ***tbl, int nRows, int nCols);
void writeCharTbl(char *inFile, const char ***tbl, int nRows, int nCols );

char* GetLine(FILE* inputfile);

#ifdef __cplusplus
}
#endif

#endif
