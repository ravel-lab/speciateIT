#ifndef MARKOVCHAINS2_HH
#define MARKOVCHAINS2_HH

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

#include <iostream>
#include <vector>
#include <map>
#include "DNAsequence.hh"

using namespace std;


// Flags for pseudo-count types
static const int zeroOffset0   = -1; ///< add 0 to all k-mer counts
static const int zeroOffset1   = 0;  ///< add 1 to all k-mer counts
static const int zeroOffset4mk = 1;  ///< add 1/4^k to k-mer counts
static const int recPdoCount   = 2;  ///< adding recursively defined pseudocount as described below

// Set the pseudocounts for a order k+1 model be alpha*probabilities from an
// order k model, recursively down to pseudocounts of alpha/num_letters for an order
// 0 model.

// Set alpha=1. Thus for order 0 model

// count[A] = count[A] + count[A]/total = count[A]*(1 + 1/total),
// count[T] = count[T]* (1 + 1/total),
// count[C] = count[C]* (1 + 1/total),
// count[G] = count[G]* (1 + 1/total),

// where

// total = count[A] + count[T] + count[C] + count[G].

// For order 1 model

// count[AA] = count[AA] + prob(A)*prob(A)
// count[AC] = count[AC] + prob(A)*prob(C)
// count[AG] = count[AG] + prob(A)*prob(G)
// count[AT] = count[AT] + prob(A)*prob(T)

// count[CA] = count[CA] + prob(C)*prob(A)
// count[CC] = count[CC] + prob(C)*prob(C)
// count[CG] = count[CG] + prob(C)*prob(G)
// count[CT] = count[CT] + prob(C)*prob(T)

// etc

// and for order 2 models

// count[ACA] = count[ACA] + prob(AC) * prob(A)
// count[ACC] = count[ACC] + prob(AC) * prob(C)
// count[ACG] = count[ACG] + prob(AC) * prob(G)
// count[ACT] = count[ACT] + prob(AC) * prob(T)

// where

// prob(AC) and prob(A) are probabilities of observing AC and A after addition of
// corresponding pseudo-counts.

// Clearly this is a generalization of 1/4^k zero-offset method, that reduces to
// 1/4^k zero-offset method for order 0 and 1 models when prob(b)= 1/4 for all bases
// b.


//============================================== MarkovChains2_t ====
/// Markov Chains probability model of order k
///
/// MC models are build for sequences in trgFiles
///
/// The class contains routines for generating and testing the models
/// To isolate training and test chromo fragments (called windows)
/// training and test window lists are generated for each chromosome
///
/// Parameters:
/// order - order of the current Markov Chains model
/// trgFiles - training fasta files; first we are going to check if models
/// corresponding to these fasta files already exist
///
class MarkovChains2_t
{
public:
  MarkovChains2_t( int order,
		   vector<char *> &trgFiles,
		   char *dir,
		   int maxNumAmbCodes=5,
		   int pseudoCountType=zeroOffset4mk );

  MarkovChains2_t( int order,
		   const char *seqID, // seq ID of a sequence to be excluded from model building (it will be used for leave-one-out validation)
		   char *&seq,        // sequence with seqID
		   int &seqLen,
		   vector<char *> &trgFiles,
		   char *dir,
		   int maxNumAmbCodes=5,
		   int pseudoCountType=zeroOffset4mk );

  ~MarkovChains2_t();

  double log10prob( const char *frag, int fragLen, int modelIdx );      // version for processing sequences without IUPAC ambiguous codes
  double log10probIUPAC( const char *frag, int fragLen, int modelIdx ); // version accepting IUPAC codes
  double log10probR( char *frag, int fragLen, int modelIdx );           // old (restart) version of log10prob()
  int log10probVect( const char *frag, int fragLen, int modelIdx, double *probs ); // computes conditional probabilities at each position of the sequence given the modelIdx-th model

  // normalized versions of the above routines where the output from the above functions is divided by the sequence length
  inline double normLog10prob( const char *frag, int fragLen, int modelIdx );
  inline double normLog10probIUPAC( const char *frag, int fragLen, int modelIdx );
  inline double normLog10probR( char *frag, int fragLen, int modelIdx );

  inline const vector<char *> & modelIds();
  inline void modelIds( vector<string> &v);
  inline const vector<vector<char *> > & wordStrgs() const;
  inline void printCounts();
  inline int order();

  void sample( const char *faFile, const char *txFile, int sampleSize, int seqLen=534 ); /// random samples from MC models
  void sample( char ***_seqTbl, int modelIdx, int sampleSize, int seqLen=534 );
  void sampleMF( char **seqTbl, int modelIdx, int sampleSize, int seqLen );
  void sample( char ***_seqTbl, map<string, string> &refSeqs, int modelIdx, int sampleSize, int seqLen );

  void printCounts( bool x ) { printCounts_m = x; }

private:
  inline int hashFn(const char *s, int sLen); /// ACGT string hash function
  void hashFnIUPAC(const char *s, int sLen, vector<int> &v);/// ACGT string hash function accepting IUPAC codes
  inline int hashFn(const char c);            /// ACGT string hash function for a single base
  inline int hashFnUL( int L );               /// upper limit for hashFn() on words of size L
  void setupIOfiles();
  bool readModelIds();
  void createModelIds();
  void createMooreMachine();
  void getKmers(int k);

  void wordCounts( int kIdx, const char *file, int modelIdx );
  void wordCountsR( int kIdx, const char *file, int modelIdx );
  void wordCountsR( int kIdx, const char *file, int modelIdx, const char *seqID);

  bool readLog10condProbTbl( int kIdx, const char *freqFile );

  void log10condProbs( int kIdx );
  void log10condProbs( int kIdx, const char *seqID );
  void log10condProbs( int kIdx, int modelIdx );

  void printLog10condProbTbl( const char *file, int kIdx, int idx);
  void printWordCountsTbl( const char *file, int kIdx, int idx );
  bool isIUPACambCode(char ch) const;
  int getIUPACambCodeHashVals(char c, vector<int> &v) const;
  void initIUPACambCodeHashVals();

  int order_m;                /// order of the current Markov Chains model
  vector<vector<char *> > wordStrgs_m;/// wordsStrgs_m[i] - vector of all words of size (i+1); to be deprecated (used only for testing)
  //-- Moore machine members: hashFn() see above
  int **tr_m;                 /// transition function; tr_m[hashFn(kmer)][nuc] = hashFn(suffix(kmer)nuc) for kmer of size order_m+1
  double **counts_m;          /// counts_m[i] - array of k-mer counts for the i-th OTU
  double **log10cProb_m;      /// log10cProb_m[i] array of log10 conditional probabilities of hashFn(kmer) for the i-th chromosome
  int *hashUL_m;              /// array of hashFnUL(k) for k=1, ... , order_m+1

  vector<string> cProbFile_m; /// vector of files containig conditional probability estimates, one for each k=1,...,(order_m+1)
  vector<string> countFile_m; /// vector of files containig k-mer counts data, one for each k=1,...,(order_m+1)

  vector<char *> trgFiles_m;  /// training fasta files
  char *dir_m;                /// input/output directory for MC model files
  vector<char *> modelIds_m;  /// list of model Ids (base names of training file names)
  int nAllWords_m;            /// number of all words of length order_m+1
  bool printCounts_m;         /// flag to print word counts
  int maxNumAmbCodes_m;       /// maximal acceptable number of ambiguity codes for a sequence; above this number log10probIUPAC() returns 1;
  int pseudoCountType_m;
  char *seq_m;                /// sequence with seqID to be selected for leave-one-out validation
  int seqLen_m;               /// length of seq_m
  vector<char *> nucs_m;      /// vector of all k-mers
  double *cProb_m;            /// variable holding posterior probabilities of a model

  // ambiguity vectors
  vector<int> Rcode_m;
  vector<int> Ycode_m;
  vector<int> Scode_m;
  vector<int> Wcode_m;
  vector<int> Kcode_m;
  vector<int> Mcode_m;
  vector<int> Bcode_m;
  vector<int> Dcode_m;
  vector<int> Hcode_m;
  vector<int> Vcode_m;
  vector<int> Ncode_m;
};

//-------------------- inlines -------------------------------
inline const vector<char *> & MarkovChains2_t::modelIds()
{ return modelIds_m; }

inline void MarkovChains2_t::modelIds( vector<string> &v )
{
  int n = modelIds_m.size();
  for ( int i = 0; i < n; i++ )
    v.push_back(string(modelIds_m[i]));
}

inline const vector<vector<char *> > & MarkovChains2_t::wordStrgs() const
{ return wordStrgs_m; }

inline int MarkovChains2_t::hashFnUL( int L )
{
  int h = 0;
  for ( int p = 1, i = 0; i < L; ++i, p *= 4 )
    h += p;

  return --h;
}

inline int MarkovChains2_t::hashFn(const char *s, int sLen)
{
  //int h = hashFnUL(sLen);
  int h = hashUL_m[sLen-1];

  for ( unsigned p = 1; *s && (intACGTLookup[int(*s)]>-1); ++s, p *= 4 )
    h += intACGTLookup[int(*s)]*p;

  if ( *s )
    h = -1;

  return h;
}

inline int MarkovChains2_t::hashFn(const char c)
{
  return intACGTLookup[int(c)];
}

inline void MarkovChains2_t::printCounts()
{
  printCounts_m = true;
}

inline double MarkovChains2_t::normLog10prob( const char *frag, int fragLen, int modelIdx )
{
  return log10prob( frag, fragLen, modelIdx ) / fragLen;
}

inline double MarkovChains2_t::normLog10probIUPAC( const char *frag, int fragLen, int modelIdx )
{
  return log10probIUPAC( frag, fragLen, modelIdx ) / fragLen;
}

inline double MarkovChains2_t::normLog10probR( char *frag, int fragLen, int modelIdx )
{
  return log10probR( frag, fragLen, modelIdx ) / fragLen;
}

inline int MarkovChains2_t::order()
{
  return order_m;
}


#endif
