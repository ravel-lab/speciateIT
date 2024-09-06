#ifndef STRINGS_H
#define STRINGS_H

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

#include <string.h>

char * strchomp(char *s);
char * strclean(char *s);
char * strlop(char *s, char ch);
char * strlopspace(char *s);
char * strpre(const char* str, const char* pre);
char * strcasepre(const char* str, const char* pre);
char * strsplit(char* str);
char * strsplit(char* str, char ch);
char * strtrim(char *s);
char * strtrimto(char *s, char ch);

#endif // STRINGS_H
