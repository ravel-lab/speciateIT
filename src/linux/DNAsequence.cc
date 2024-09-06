/*
Copyright (C) 2016 Pawel Gajer pgajer@gmail.com, Adam M Phillippy and Jacques Ravel jravel@som.umaryland.edu

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

#include "DNAsequence.hh"
#include "CUtilities.h"

#include <iostream>
using namespace std;

//------------------------------------------------- getKmers ----
/// returns vector of all k-mers of size k
void getAllKmers( int k, vector<char *> &kmers )
{
  vector<char *> nucs;
  nucs.push_back(strdup("A"));
  nucs.push_back(strdup("C"));
  nucs.push_back(strdup("G"));
  nucs.push_back(strdup("T"));

  if ( k == 1 )
  {
    kmers = nucs;
    return;
  }

  vector<char *> km1mers;
  getAllKmers(k-1,km1mers);

  int n = km1mers.size();
  char str[100];
  char *kmer;

  for ( int i = 0; i < n; ++i )
  {
    for ( int j = 0; j < 4; ++j )
    {
      str[0] = '\0';
      strcat(str,nucs[j]);
      strcat(str,km1mers[i]);
      STRDUP(kmer,str);
      kmers.push_back(kmer);
    }
    free(km1mers[i]);
  }

  for ( int i = 0; i < 4; ++i )
    free(nucs[i]);
}
