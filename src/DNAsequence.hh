#ifndef DNASEQUENCE_HH
#define DNASEQUENCE_HH

/*
Copyright (C) 2016 Pawel Gajer (pgajer@gmail.com), Adam M Phillippy and Jacques Ravel jravel@som.umaryland.edu

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

// addopted from inVUE Adam's Sequence.hh

#include <vector>
using namespace std;

const int intACGTLookup[] =
  {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   -1,0,-1,1,-1,-1,-1,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   -1,0,-1,1,-1,-1,-1,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

// IUPAC nucleotide code
// R 	A or G
// Y 	C or T
// S 	G or C
// W 	A or T
// K 	G or T
// M 	A or C
// B 	C or G or T
// D 	A or G or T
// H 	A or C or T
// V 	A or C or G
// N 	any base


struct win_t
{
  int start;
  int end;
};

const char ComplementLookup[] =
  "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn"
  " nnnnnnnnn*nn-.nnnnnnnnnnnnnnnnn"
  "nTVGHNNCDNNMNKNNNNYSANBSNRNnnnn_"
  "ntvghnncdnnmnknnnnysanbsnrnnnnnn";

const bool ACGTLookup[] =
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
   0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0};

inline char Complement(char ch)
{
  return ComplementLookup[int(ch)];
}

inline bool isACGT(char ch)
{
  return ACGTLookup[int(ch)];
}

inline bool isACGT(char *s)
{
  for ( ; *s && ACGTLookup[int(*s)]; ++s )
  {}

  return !(bool)*s;
}

void getAllKmers( int k, vector<char *> &kmers );

#endif
