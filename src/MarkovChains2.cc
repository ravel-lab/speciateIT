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

#include <string>
#include <algorithm>
#include "MarkovChains2.hh"
#include "CStatUtilities.h"
#include "CUtilities.h"
#include "IOCppUtilities.hh"
#include "IOCUtilities.h"
#include "CppStatUtilities.hh"
#include "CppUtilities.hh"
#include "strings.hh"
#include "DNAsequence.hh"

using namespace std;

//========================================================= MarkovChains2_t ====
/*
  March 10, 2009 changes:

   - added maxNumAmbCodes parameter that controls the number of acceptable ambiguity codes,
     before the sequence is dropped form log10prob() and printed (in buildMC.cc)
     to low quality seq's file.

  March 7, 2009 changes:

  - log10prob() has been moved to log10probR()

  - the new log10prob() assumes that the sequence does not have ambiguity codes
    if they are detected, the calculation of log10 prob. it passed to log10probIUPAC()

  - log10probIUPAC() averages over k-mers resulting from replacing ambiguous codes
    by the corresponding bases

  - by default, lower order MC probabilities for the initial fragment
    of the sequence are not computed (to avoid bias in prob difference due to
    mutations in the initial fragment of the sequence)

  - modified MC model building so that ambiguity codes are replaced by the corresonding
    bases

*/


//--------------------------------------------------------- MarkovChains2_t ----
/// MarkovChains2_t attempts to read relative k-mer frequencies
MarkovChains2_t::MarkovChains2_t(int order,
				 vector<char *> &trgFiles,
				 char *dir,
				 int maxNumAmbCodes,
				 int pseudoCountType)
  : order_m(order), dir_m(dir),  maxNumAmbCodes_m(maxNumAmbCodes), pseudoCountType_m(pseudoCountType)
{
  int maxWordLen = order_m+1;

  getAllKmers( 1, nucs_m );

  // wordStrgs_m is needed only for the header of cProbFile_m's
  wordStrgs_m.resize(maxWordLen);
  for ( int k = 1; k <= maxWordLen; ++k )
    getAllKmers(k, wordStrgs_m[k-1]);

  //MALLOC(cProb_m, double*, nAllWords_m * sizeof(double)); // pp_ref_sib_wr_ref_models(15223,0x7fff7f0b4300) malloc: *** mach_vm_map(size=18446744072073781248) failed (error code=3)
  #if 0
  cerr << "in MarkovChains2_t::MarkovChains2_t() order_m=" << order_m
       << "\tmaxWordLen=" << maxWordLen
       << endl;
  cerr << "wordStrgs_m[0]:" << endl;
  printVector(wordStrgs_m[0]);

  cerr << "wordStrgs_m[1]:" << endl;
  printVector(wordStrgs_m[1]);
  #endif

  char *file;
  int nTrgFiles = trgFiles.size();
  for ( int i = 0; i < nTrgFiles; ++i )
  {
    STRDUP(file, trgFiles[i]);
    trgFiles_m.push_back(file);
  }

  if ( dir_m )
  {
    setupIOfiles();

    if ( !readModelIds() && !nTrgFiles )
    {
      cerr << "Error in MarkovChains2_t::MarkovChains2_t(): cannot read model ids" << endl;
      exit(1);
    }
    else if ( !readModelIds() && nTrgFiles )
    {
      createModelIds();
    }

    createMooreMachine();

    MALLOC(cProb_m, double*, nAllWords_m * sizeof(double));

    initIUPACambCodeHashVals();
    bool ok = false;

    if ( !(ok=readLog10condProbTbl(0, cProbFile_m[0].c_str())) && nTrgFiles )
    {
      log10condProbs( 0 );
    }

    if ( 0 & !ok )
    {
      cerr << "Error in MarkovChains2_t::MarkovChains2_t(): cannot read " << cProbFile_m[0].c_str()
	   << "\nand cannot create Log10condProbTbl" << endl;
      exit(1);
    }

    for ( int i = 1; i < maxWordLen; ++i )
    {
      if ( !(ok=readLog10condProbTbl(i, cProbFile_m[i].c_str())) && nTrgFiles )
      {
	log10condProbs( i );
      }
      #if 0
      if ( !ok )
      {
	cerr << "Error in " << __FILE__ << " at line " << __LINE__ << ": cannot read and cannot create Log10condProbTbl" << endl;
	exit(1);
      }
      #endif
    }
  }
  else
  {
    createModelIds();
    createMooreMachine();
    MALLOC(cProb_m, double*, nAllWords_m * sizeof(double));
    initIUPACambCodeHashVals();

    for ( int i = 0; i < maxWordLen; ++i )
      log10condProbs( i );
  }
}


//--------------------------------------------------------- MarkovChains2_t ----
/// MarkovChains2_t version for leave-one-out validation setup (dir is NULL)
MarkovChains2_t::MarkovChains2_t(int order,
				 const char *seqID, // seq ID of a sequence to be excluded from model building (it will be used for leave-one-out validation)
				 char *&seq,        // sequence with seqID
				 int &seqLen,
				 vector<char *> &trgFiles,
				 char *dir,
				 int maxNumAmbCodes,
				 int pseudoCountType)
  : order_m(order), dir_m(dir),  maxNumAmbCodes_m(maxNumAmbCodes), pseudoCountType_m(pseudoCountType)
{
  int maxWordLen = order_m+1;

  getAllKmers( 1, nucs_m );

  // wordStrgs_m is needed only for the header of cProbFile_m's
  wordStrgs_m.resize(maxWordLen);
  for ( int k = 1; k <= maxWordLen; ++k )
    getAllKmers(k, wordStrgs_m[k-1]);

  char *file;
  int nTrgFiles = trgFiles.size();
  for ( int i = 0; i < nTrgFiles; ++i )
  {
    STRDUP(file, trgFiles[i]);
    trgFiles_m.push_back(file);
  }

  createModelIds();
  createMooreMachine();
  initIUPACambCodeHashVals();

  for ( int i = 0; i < maxWordLen; ++i )
    log10condProbs( i, seqID );

  seq = seq_m;
  seqLen = seqLen_m;
}


//-------------------------------------------------------- readModelIds ----
bool MarkovChains2_t::readModelIds()
{
  bool readModels = true;

  if ( !dir_m )
  {
    return 0;
  }
  else
  {
    string inFile(dir_m);
    inFile += "/modelIds.txt";
    FILE *in = fopen(inFile.c_str(), "r");

    if ( !in )
      return 0;

    fclose(in);

    readLines(inFile.c_str(), modelIds_m);
  }

  return readModels;
}

//-------------------------------------------------------- createModelIds ----
void MarkovChains2_t::createModelIds()
{
  int n = trgFiles_m.size();
  char *idStr;

  for ( int i = 0; i < n; ++i )
  {
    string id = baseFileName(trgFiles_m[i]);
    STRDUP(idStr,id.c_str());
    modelIds_m.push_back(idStr);
  }

  if ( dir_m )
  {
    string outFile(dir_m);
    outFile += "/modelIds.txt";
    FILE *out = fOpen(outFile.c_str(), "w");
    for ( int i = 0; i < n; ++i )
      fprintf(out,"%s\n",modelIds_m[i]);
    fclose(out);
  }

  #if 0
  cerr << "in createModelIds(): modelIds_m=";
  for ( int i = 0; i < n; ++i )
    cerr << "\t" << modelIds_m[i];
  cerr << endl;
  #endif
}

//-------------------------------------------------------- ~MarkovChains2_t ----
MarkovChains2_t::~MarkovChains2_t()
{
  int nModels = modelIds_m.size();
  unsigned n = 0;

  if ( log10cProb_m )
  {
    for ( int i = 0; i < nModels; ++i )
      free(log10cProb_m[i]);
    free(log10cProb_m);
  }

  if ( counts_m )
  {
    for ( int i = 0; i < nModels; ++i )
      free(counts_m[i]);
    free(counts_m);
  }

  if ( tr_m )
  {
    for ( int i = 0; i < nAllWords_m; ++i )
      free(tr_m[i]);
    free(tr_m);
  }

  if ( hashUL_m )
    free(hashUL_m);

  if ( (n=wordStrgs_m.size()) )
    for ( unsigned i = 0; i < n; ++i )
      for ( unsigned j = 0; j < wordStrgs_m[i].size(); ++j )
	free(wordStrgs_m[i][j]);

  if ( (n=modelIds_m.size()) )
    for ( unsigned i = 0; i < n; ++i )
      free(modelIds_m[i]);

  for ( int i = 0; i < 4; ++i )
    free(nucs_m[i]);

  free(cProb_m);
}

//------------------------------------------------------------------- setupIOfiles ----
/// setupIOfiles intialize basic parameter strings and file names
void MarkovChains2_t::setupIOfiles()
{
  if ( dir_m )
    printCounts_m = false;

  srand ( time(NULL) );

  char orderStr[10];
  sprintf(orderStr,"%d",order_m);

  int n = order_m+1;
  char countStr[5];

  string cmd("mkdir -p ");
  cmd += dir_m;
  system(cmd.c_str());

  for ( int j = 0; j < n; ++j )
  {
    sprintf(countStr,"%d",j);

    cProbFile_m.push_back( string(dir_m) + string("/MC") + string(countStr)
			   + string(".log10cProb") );

    countFile_m.push_back( string(dir_m) + string("/MC") + string(countStr)
			   + string(".count") );
  }

  #if 0
  cerr << "in setupIOfiles()" << endl;
  cerr << "cProbFile_m:\t";
  for ( int j = 0; j < n; ++j )
    cerr << cProbFile_m[j].c_str() << "\t";
  cerr << endl;

  cerr << "countFile_m:\t";
  for ( int j = 0; j < n; ++j )
    cerr << countFile_m[j].c_str() << "\t";
  cerr << endl;
  #endif
}

//----------------------------------------------------- createMooreMachine ----
void MarkovChains2_t::createMooreMachine( )
/// creating Moore machine as describe
/// in W.H. Majoros "Computatinal Gene Prediction" Section 7.2.3
{
  #define DEBUG_MOORE 0

  int nNucs = 4;
  int maxWordLen = order_m+1;
  int maxWordLen1 = maxWordLen+1;

  #if DEBUG_MOORE
  cerr << "order_m: " << order_m << endl;
  cerr << "maxWordLen1: " << maxWordLen1 << endl;
  #endif

  MALLOC(hashUL_m, int*, maxWordLen1 * sizeof(int));

  for ( int i = 0; i < maxWordLen1; ++i )
  {
    hashUL_m[i] = hashFnUL(i+1);
    #if DEBUG_MOORE
    cerr << "hashUL_m["  << i << "]=" << hashUL_m[i] << endl;
    #endif
  }
  #if DEBUG_MOORE
  cerr << endl;
  #endif

  nAllWords_m = hashFnUL(order_m+2);

  int nModels = modelIds_m.size();
  MALLOC(counts_m, double**, nModels * sizeof(double*));
  MALLOC(log10cProb_m, double**, nModels * sizeof(double*));

  for ( int i = 0; i < nModels; ++i )
  {
    CALLOC(counts_m[i], double*, nAllWords_m * sizeof(double));
    MALLOC(log10cProb_m[i], double*, nAllWords_m * sizeof(double));
  }

  //-- setting up tr_m
  MALLOC(tr_m, int**, nAllWords_m * sizeof(int*));

  for ( int i = 0; i < nAllWords_m; ++i )
    MALLOC(tr_m[i], int*, nNucs * sizeof(int));

  int l, p4, lL, lU;

  for( l = 1, p4 = 4; l < maxWordLen; ++l, p4 *= 4 )
  {
    lL = hashUL_m[l-1]; // lower limit for hash values of k-mers of length l
    lU = hashUL_m[l];   // upper limit for hash values of k-mers of length l

    #if DEBUG_MOORE
    cerr << endl << "> l=" << l << "\tlL=" << lL << "\tlU=" << lU << "\tp4=" << p4 << endl;
    #endif

    for ( int v = lL; v < lU; ++v )
    {
      for ( int i = 0; i < nNucs; ++i )
      {
	tr_m[v][i] = v + p4*(1 + i);
	#if DEBUG_MOORE
	cerr << "tr_m[" << v << "][" << i << "]=" << tr_m[v][i] << endl;
	#endif
      }
    }
  }

  lL = hashUL_m[l-1];
  lU = hashUL_m[l];
  p4 /= 4;

  #if DEBUG_MOORE
  cerr << "order_m=" << order_m << "\tmaxWordLen=" << maxWordLen << endl;
  cerr << "l=" << l << "\tlL=" << lL
       << "\tlU=" << lU
       << "\tp4=" << p4
       << endl;
  #endif

  for ( int v = lL; v < lU; ++v )
    for ( int i = 0; i < nNucs; ++i )
    {
      tr_m[v][i] = lL + (v - lL)/4 + p4*i;
      #if DEBUG_MOORE
      cerr << "v=" << v << "\ti=" << i
	   << "\t(v - lL)/4=" << (v - lL)/4
	   << "\tp4*i=" << p4*i
	   << "\ttr_m[" << v << "][" << i << "]=" << tr_m[v][i]
	   << endl;
      #endif
    }

  #if DEBUG_MOORE
  cerr << "tr_m:" << endl;
  for ( int i = 0; i < nAllWords_m; ++i )
    for ( int j = 0; j < 4; ++j )
      cerr << "tr_m[" << i << "][" << j << "]=" << tr_m[i][j] << endl;
  exit(1);
  #endif
}

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

//------------------------------------------------ isIUPACambCode -----------
/// check if the character is one of IUPAC ambiguity codes
bool MarkovChains2_t::isIUPACambCode(char c) const
{
  switch(c)
  {
    case 'R':
    case 'Y':
    case 'S':
    case 'W':
    case 'K':
    case 'M':
    case 'B':
    case 'D':
    case 'H':
    case 'V':
    case 'N':
      return true;
      break;
    default:
      cerr << "Error in MarkovChains2_t::isIUPACambCode(): Unrecognized IUPAC character: "
	   << (char)c << endl;
      return false;
  }
}

//------------------------------------------------ initIUPACambCodeHashVals -----------
/// initialize IUPAC ambiguity codes hash value vectors
void MarkovChains2_t::initIUPACambCodeHashVals()
{
  Rcode_m.push_back(intACGTLookup[int('A')]);
  Rcode_m.push_back(intACGTLookup[int('G')]);

  Ycode_m.push_back(intACGTLookup[int('C')]);
  Ycode_m.push_back(intACGTLookup[int('T')]);

  Scode_m.push_back(intACGTLookup[int('C')]);
  Scode_m.push_back(intACGTLookup[int('G')]);

  Wcode_m.push_back(intACGTLookup[int('A')]);
  Wcode_m.push_back(intACGTLookup[int('T')]);

  Kcode_m.push_back(intACGTLookup[int('T')]);
  Kcode_m.push_back(intACGTLookup[int('G')]);

  Mcode_m.push_back(intACGTLookup[int('A')]);
  Mcode_m.push_back(intACGTLookup[int('C')]);

  Bcode_m.push_back(intACGTLookup[int('C')]);
  Bcode_m.push_back(intACGTLookup[int('G')]);
  Bcode_m.push_back(intACGTLookup[int('T')]);

  Dcode_m.push_back(intACGTLookup[int('A')]);
  Dcode_m.push_back(intACGTLookup[int('G')]);
  Dcode_m.push_back(intACGTLookup[int('T')]);

  Hcode_m.push_back(intACGTLookup[int('A')]);
  Hcode_m.push_back(intACGTLookup[int('C')]);
  Hcode_m.push_back(intACGTLookup[int('T')]);

  Vcode_m.push_back(intACGTLookup[int('A')]);
  Vcode_m.push_back(intACGTLookup[int('C')]);
  Vcode_m.push_back(intACGTLookup[int('G')]);

  Ncode_m.push_back(intACGTLookup[int('A')]);
  Ncode_m.push_back(intACGTLookup[int('C')]);
  Ncode_m.push_back(intACGTLookup[int('G')]);
  Ncode_m.push_back(intACGTLookup[int('T')]);
}

//------------------------------------------------ getIUPACambCodeHashVals -----------
/// gets hash values of IUPAC ambiguity code bases
/// and returns the number of bases corresponding to the code
int MarkovChains2_t::getIUPACambCodeHashVals(char c, vector<int> &v) const
{
  v.clear();

  switch(c)
  {
    case 'R':
      v = Rcode_m;
      return 2;
      break;
    case 'Y':
      v = Ycode_m;
      return 2;
      break;
    case 'S':
      v = Scode_m;
      return 2;
      break;
    case 'W':
      v = Wcode_m;
      return 2;
      break;
    case 'K':
      v = Kcode_m;
      return 2;
      break;
    case 'M':
      v = Mcode_m;
      return 2;
      break;
    case 'B':
      v = Bcode_m;
      return 3;
      break;
    case 'D':
      v = Dcode_m;
      return 3;
      break;
    case 'H':
      v = Hcode_m;
      return 3;
      break;
    case 'V':
      v = Vcode_m;
      return 3;
      break;
    case 'N':
      v = Ncode_m;
      return 4;
      break;
    default:
      cerr << "Error in MarkovChains2_t::getIUPACambCodeHashVals(): Unrecognized IUPAC character: "
	   << (char)c << endl;
      v.clear();
      return 0;
  }
}

//------------------------------------------------- kmerState_t -----------
/// hold k-mer hashFn() index and cumulative log10 probability
struct kmerState_t
{
  int kmerIdx;
  double log10prob;
};


//------------------------------------------------- log10probIUPAC -----------
/// computes a Markov Chains estimate of log10 probability that frag
/// comes from i-th model, where i=modelIdx
/// if ambiguous IUPAC codes are encountered it averages over sequences
/// with the code being replaced by the corresponding nucleotides
double MarkovChains2_t::log10probIUPAC( const char *frag, int fragLen, int modelIdx )
{
  int k = 0, i;
  int codeSize;
  vector<int> codes;
  kmerState_t kmerState;
  vector<kmerState_t> thread; /// each ambiguity code changes thread.size() to m*thread.size()
  /// where m is the number of bases corresponding to the ambiguous code
  /// for example R = A or G, has m=2.

  int nAmbCodes = 0;

  //-- skip the first k bases to avoid lower MC probability calculations
  if ( (i=intACGTLookup[int(frag[k])]) > -1 )
  {
    kmerState.kmerIdx   = i+1;
    kmerState.log10prob = 0;
    thread.push_back(kmerState);
  }
  else if ( (codeSize=getIUPACambCodeHashVals(frag[k], codes)) )
  {
    nAmbCodes++;

    for ( i = 0; i < codeSize; ++i )
    {
      kmerState.kmerIdx = codes[i]+1;
      kmerState.log10prob = 0;
      thread.push_back(kmerState);
    }
  }
  else
  {
    return 1; // log10 of probability has to be <= 0, so returned value 1 means error
  }

  k++;
  int rank = order_m+1;

  while ( k < rank )
  {
    if ( nAmbCodes > maxNumAmbCodes_m )
    {
      //printf("\nfrag=%s\tnAmbCodes=%d\n",frag,nAmbCodes);
      return 1;
    }

    int n = thread.size();

    if ( (i=intACGTLookup[int(frag[k])]) > -1 )
    {
      for ( int j = 0; j < n; ++j )
	thread[j].kmerIdx = tr_m[thread[j].kmerIdx][i];
    }
    else if ( (codeSize=getIUPACambCodeHashVals(frag[k], codes)) )
    {
      nAmbCodes++;

      for ( int j = 0; j < n; ++j )
	thread[j].kmerIdx = tr_m[thread[j].kmerIdx][codes[0]];

      for ( i = 1; i < codeSize; ++i )
      {
	for ( int j = 0; j < n; ++j )
	{
	  kmerState.kmerIdx = tr_m[thread[j].kmerIdx][codes[i]];
	  kmerState.log10prob = 0;
	  thread.push_back(kmerState);
	}
      }
    }
    else
    {
      return 1; // log10 of probability has to be <= 0, so returned value 1 means error
    }

    k++;
  }

  //-- process the remaining k-mers
  while ( k < fragLen )
  {
    if ( nAmbCodes > maxNumAmbCodes_m )
    {
      //printf("\nfrag=%s\tnAmbCodes=%d\n",frag,nAmbCodes);
      return 1;
    }

    int n = thread.size();

    if ( (i=intACGTLookup[int(frag[k])]) > -1 )
    {
      for ( int j = 0; j < n; ++j )
      {
	thread[j].kmerIdx = tr_m[thread[j].kmerIdx][i];
	thread[j].log10prob += log10cProb_m[modelIdx][thread[j].kmerIdx];
      }
    }
    else if ( (codeSize=getIUPACambCodeHashVals(frag[k], codes)) )
    {
      nAmbCodes++;

      for ( int j = 0; j < n; ++j )
      {
	thread[j].kmerIdx = tr_m[thread[j].kmerIdx][codes[0]];
	thread[j].log10prob += log10cProb_m[modelIdx][thread[j].kmerIdx];
      }

      for ( i = 1; i < codeSize; ++i )
      {
	for ( int j = 0; j < n; ++j )
	{
	  kmerState.kmerIdx = tr_m[thread[j].kmerIdx][codes[i]];
	  kmerState.log10prob = thread[j].log10prob + log10cProb_m[modelIdx][kmerState.kmerIdx];
	  thread.push_back(kmerState);
	}
      }
    }
    else
    {
      return 1; // log10 of probability has to be <= 0, so returned value 1 means error
    }

    k++;
  }

  //-- compute the mean of all thread log10 values
  int n = thread.size();
  double meanLog10prob = 0;
  for ( i = 0; i < n; ++i )
    meanLog10prob += thread[i].log10prob;

  //cerr << "in log10probIUPAC() thread.size()=" << n << endl;

  return meanLog10prob / n;
}

// ---------------------------------------------------- log10prob -----------
/// computes a Markov Chains estimate of log10 probability that frag
/// comes from i-th model, where i=modelIdx
/// when IUPAC ambiguous code is encoutered, it passes log10 probability
/// calculation to log10probIUPAC()
/// by default, lower order MC probabilities for the initial fragment
/// of the sequence are not computed (to avoid bias in prob difference due to
/// mutations in the initial fragment of the sequence)
double MarkovChains2_t::log10prob( const char *frag, int fragLen, int modelIdx )
{
  double log10probVal = 0;
  int k = 0, v, i;

  //-- skip the first k bases to avoid lower MC probability calculations
  if ( (v=intACGTLookup[int(frag[k])]) > -1 )
  {
    v++;
    k++;
  }
  else
  {
    return log10probIUPAC(frag, fragLen, modelIdx);
  }

  int rank = order_m+1;

  while ( k < rank )
  {
    if ( (i=intACGTLookup[int(frag[k])]) > -1 )
    {
      v = tr_m[v][i];
    }
    else
    {
      return log10probIUPAC(frag, fragLen, modelIdx);
    }

    k++;
  }

  //-- process the remaining k-mers
  while ( k < fragLen )
  {
    if ( (i=intACGTLookup[int(frag[k])]) > -1 )
    {
      v = tr_m[v][i];
      log10probVal += log10cProb_m[modelIdx][v];
    }
    else
    {
      return log10probIUPAC(frag, fragLen, modelIdx);
    }
    k++;
  }

  return log10probVal;
}

// ---------------------------------------------------- log10probVect -----------

/// computes conditional probabilities at each position of the sequence given the
/// i-th model

/// when IUPAC ambiguous code is encoutered, we set the conditional probabilities
/// to 0

/// the conditional probabilities for the first orderMC initial positions are not
//  computed

// the user is responsible for allocation of memory for 'probs' vector

// returns the number of positions processed
//
int MarkovChains2_t::log10probVect( const char *frag, int fragLen, int modelIdx, double *probs )
{
  int v = 0, i = 0;

  //-- skip the first k bases to avoid lower MC probability calculations
  int rank = order_m+1;
  int k = rank;

  //-- process the remaining k-mers
  while ( k < fragLen )
  {
    if ( (i=intACGTLookup[int(frag[k])]) > -1 )
    {
      v = tr_m[v][i];
      probs[ k - rank ] = log10cProb_m[modelIdx][v];
    }
    else
    {
      probs[ k - rank ] = 0;
    }
    k++;
  }

  return k - rank;
}

// ---------------------------------------------------- log10probR -----------
/// computes a Markov Chains estimate of log10 probability that frag
/// comes from i-th model, where i=modelIdx
/// when IUPAC ambiguous code is encoutered all k-mers containing the base
/// are discarded and the probability calculation is restart after
/// the base as if it was the begininig of the sequence
double MarkovChains2_t::log10probR( char *frag, int fragLen, int modelIdx )
{
  #define LOG10CPROBR_DEBUG 0
  #if LOG10CPROBR_DEBUG
   cerr << "in MarkovChains2_t::log10prob()\tmodelIdx=" << modelIdx
       << "\tfragLen=" << fragLen << endl;
   cerr << "intACGTLookup[A]=" << intACGTLookup[int('A')] // 0
	<< "intACGTLookup[C]=" << intACGTLookup[int('C')] // 1
	<< "intACGTLookup[G]=" << intACGTLookup[int('G')] // 2
	<< "intACGTLookup[T]=" << intACGTLookup[int('T')] // 3
	<< endl;
  #endif

  double log10probVal = 0;
  char nuc[2];
  nuc[1] = '\0';
  int k = 0, v, i;

  while ( k < fragLen )
  {
    // find the first letter in frag which is from {A,C,G,T}
    while ( k < fragLen && intACGTLookup[int(frag[k])] == -1 )
      k++;

    nuc[0] = frag[k];
    v = hashFn(nuc,1);
    log10probVal += log10cProb_m[modelIdx][v];
    k++;

    #if LOG10CPROBR_DEBUG
    cerr << "k=" << k
	 << "\tnuc=" << nuc
	 << "\tv=hashFn(nuc)=" << v
	 << "\tlog10cProb[v]=" << log10cProb[v] << endl;
    #endif

    while ( k < fragLen && (i=intACGTLookup[int(frag[k])]) > -1 )
    {
      v = tr_m[v][i];
      log10probVal += log10cProb_m[modelIdx][v];

      #if LOG10CPROBR_DEBUG
      cerr << "k=" << k
	   << "\tfrag[k]=" << (char)frag[k]
	   << "\ti=idx[frag[k]]=" << i
	   << "\ttr_m[v][i]=" << v
	   << "\tlog10cProb_m[modelIdx][v]=" << log10cProb_m[modelIdx][v] << endl;
      #endif

      k++;
    }
  }

  return log10probVal;
}


//-------------------------------------------------- readLog10condProbTbl ----
/// Reads a log10 frequencies table
///
/// chromoID kmer1       kmer2       ... kmerN
/// chromo1  log10(f11)  log10(f12)      log10(f1N)
/// chromo2  log10(f21)  log10(f22)      log10(f2N)
/// ...
/// ...
/// chromoM  log10(fM1)  log10(fM2)      log10(fMN)
///
/// stored in inFile. Returns true on success, otherwise false.
/// kIdx = k-1
///
bool MarkovChains2_t::readLog10condProbTbl( int kIdx, const char *inFile )
{
  FILE *fp = fopen(inFile, "r");

  if ( !fp )
  {
    string cmd = string("rm -f ") + string(inFile);
    system(cmd.c_str());
    return false;
  }

  //-- Read the header
  size_t alloc = 1024*1024;
  char *data;
  MALLOC(data, char*, alloc * sizeof(char));

  if ( !ReadLine(fp, data, alloc) )
  {
    cerr << "Error: MarkovChains2_t::readLog10condProbTbl() at line " << __LINE__
	 << " could not read the data"  << endl;
    return false;
  }

  //-- Parse k-mers from the header
  char *brkt;
  char *word = strtok_r(data, "\t", &brkt);
  //char *kmer;

  for ( word = strtok_r(NULL, "\t", &brkt);// starting from second field
	word;                              // to get rid of chromoId label
	word = strtok_r(NULL, "\t", &brkt))
  {
    //STRDUP(kmer,chomp(word));
    //wordStrgs_m[kIdx].push_back(kmer);
  }

  int nModels = modelIds_m.size();
  int modelIdx = 0;
  int lL = hashUL_m[kIdx];  // lower bound for hash val's of words of size kIdx+1
  int lU = hashUL_m[kIdx+1];// upper bound for hash val's of words of size kIdx+2

  while ( ReadLine(fp, data, alloc) )
  {
    if ( modelIdx == nModels )
    {
      nModels += 100;
      REALLOC(log10cProb_m, double**, nModels * sizeof(double*));
    }

    //-- parse data; extract chromoId and log10 frequencies
    word = strtok_r(data, "\t", &brkt);

    int i = lL;
    for ( word = strtok_r(NULL, "\t", &brkt);
	  word;
	  word = strtok_r(NULL, "\t", &brkt), ++i)
    {
      if ( i == lU )
      {
	cerr << "Error: MarkovChains2_t::readLog10condProbTbl() at line " << __LINE__
	     << " i == nWords; i=" << i
	     << "\tnWords=" << (lU-lL)  << endl;
	return false;
      }

      log10cProb_m[modelIdx][i] = strtod(word, (char **)NULL);
    }
    modelIdx++;
  }

  //cerr << "Done parsing " << inFile << endl;

  fclose(fp);
  free(data);

  return true;
}


//------------------------------------------------- log10condProbs ----
void MarkovChains2_t::log10condProbs( int kIdx )
{
  int nModels = modelIds_m.size();

  for ( int i = 0; i < nModels; ++i )
  {
    //cerr << "Generating " << (kIdx+1) << "-mer FreqTbl for " << trgFiles_m[i] << endl;;

    //-- compute word counts and rel. frequencies
    wordCountsR( kIdx, trgFiles_m[i], i );

    if ( dir_m && printCounts_m )
      printWordCountsTbl( countFile_m[kIdx].c_str(), kIdx, i );

    log10condProbs( kIdx, i );

    if ( dir_m )
      printLog10condProbTbl( cProbFile_m[kIdx].c_str(), kIdx, i );

  }
}

//------------------------------------------------- log10condProbs ----
void MarkovChains2_t::log10condProbs( int kIdx, const char *seqID )
{
  int nModels = modelIds_m.size();

  for ( int i = 0; i < nModels; ++i )
  {
    //-- compute word counts and rel. frequencies
    wordCountsR( kIdx, trgFiles_m[i], i, seqID );

    //fprintf(stderr, "in log10condProbs() seq: %s\n", seq);

    if ( dir_m && printCounts_m )
      printWordCountsTbl( countFile_m[kIdx].c_str(), kIdx, i );

    log10condProbs( kIdx, i );

    if ( dir_m )
      printLog10condProbTbl( cProbFile_m[kIdx].c_str(), kIdx, i );

  }
}

//------------------------------------------------------ wordCounts ----
/// computes word counts for specified kIdx (wordLen=kIdx+1)
/// and places them in log10cProb_m array
/// pseudo-count is added to each word count
/// when ambiguity code is encountered, code's alternative k-mers are counted
void MarkovChains2_t::wordCounts( int kIdx,
				  const char *file,
				  int modelIdx )
{
  int wordLen = kIdx+1;
  int lL = hashUL_m[wordLen-1];
  int lU = hashUL_m[wordLen];

  char *id, *seq;
  int seqLen;
  vector<int> idxs;

  FILE *fp = fOpen(file, "r");

  while ( getNextFastaRecord( fp, id, seq, seqLen) )
  {
    int nWords = seqLen - wordLen + 1;
    for ( int j = 0; j < nWords; ++j )
    {
      char *word = seq + j;
      char c = word[wordLen];
      word[wordLen] = '\0';
      int hashIdx = hashFn(word, wordLen);

      if ( hashIdx > -1 )
      {
	counts_m[modelIdx][ hashIdx ]++;
      }
      else
      {
	hashFnIUPAC(word, wordLen, idxs);

	int n = idxs.size();
	for ( int j = 0; j < n; ++j )
	  counts_m[modelIdx][ idxs[j] ]++;
      }

      word[wordLen] = c;
    }

    free(id);
    free(seq);
  }

  if ( pseudoCountType_m == zeroOffset1 )
  {
    double pseudoCount = 1;
    for ( int i = lL; i < lU; ++i )
      counts_m[modelIdx][i] += pseudoCount;
  }
  else if ( pseudoCountType_m == zeroOffset0 )
  {
    #if 0
    cerr << "in pseudoCountType_m == zeroOffset0" << endl;
    for ( int i = lL; i < lU; ++i )
      cerr << counts_m[modelIdx][i] << endl;
    #endif
  }
  else if ( pseudoCountType_m == zeroOffset4mk )
  {
    double pseudoCount = 1.0/pow(4.0, wordLen) ;
    for ( int i = lL; i < lU; ++i )
      counts_m[modelIdx][i] += pseudoCount;
  }
  else if ( pseudoCountType_m == recPdoCount )
  {
    if ( wordLen==1 )
    {
      // for order 0 model

      // count[A] = count[A] + count[A]/total = count[A]*(1 + 1/total),
      // count[T] = count[T]* (1 + 1/total),
      // count[C] = count[C]* (1 + 1/total),
      // count[G] = count[G]* (1 + 1/total),

      double totalCount = 0;
      for ( int i = lL; i < lU; ++i )
	totalCount += counts_m[modelIdx][i];

      double pseudoCountFactor = 1.0 + 1.0/totalCount;
      for ( int i = lL; i < lU; ++i )
	counts_m[modelIdx][i] *= pseudoCountFactor;
    }
    else
    {
      // For order 1 model

      // count[AA] = count[AA] + prob(A)*prob(A)
      // count[AC] = count[AC] + prob(A)*prob(C)
      // count[AG] = count[AG] + prob(A)*prob(G)
      // count[AT] = count[AT] + prob(A)*prob(T)
      // etc

      // and for order 2 models

      // count[ACA] = count[ACA] + prob(AC) * prob(A)
      // count[ACC] = count[ACC] + prob(AC) * prob(C)
      // count[ACG] = count[ACG] + prob(AC) * prob(G)
      // count[ACT] = count[ACT] + prob(AC) * prob(T)
      // etc

      // where

      // prob(AC) and prob(A) are probabilities of observing AC and A after addition of
      // corresponding pseudo-counts.

      lL = hashUL_m[wordLen-2];
      lU = hashUL_m[wordLen-1];
      int nNucs  = 4;
      int w;

      double totalCount = 0;
      for ( int i = lL; i < lU; ++i )
	totalCount += counts_m[modelIdx][i];

      for ( int v = lL; v < lU; ++v )
      {
	for ( int i = 0; i < nNucs; ++i )
	{
	  w = tr_m[v][i];
	  //counts_m[modelIdx][w] += pow(10.0, log10cProb_m[modelIdx][v]) * pow(10.0, log10cProb_m[modelIdx][i]);
	  // log10cProb_m[modelIdx][v] is log10 of conditional probability not the probability of v
	  // here is a correct formula
	  counts_m[modelIdx][w] += counts_m[modelIdx][v]/totalCount * pow(10.0, log10cProb_m[modelIdx][i]);
	}
      }
    }
  }

  fclose(fp);
}

//------------------------------------------------------ wordCountsR ----
/// old version of wordCounts() that skips k-mers with ambiguity codes
/// computes word counts for specified kIdx (wordLen=kIdx+1)
/// and places them in log10cProb_m array
/// pseudo-count is added to each word count
void MarkovChains2_t::wordCountsR( int kIdx,
				   const char *file,
				   int modelIdx )
{
  #define GETWORDCOUNTSR_DEBUG 0

  int wordLen = kIdx+1;
  char *id, *seq;
  int seqLen;

  FILE *fp = fOpen(file, "r");

  while ( getNextFastaRecord( fp, id, seq, seqLen) )
  {
    int nWords = seqLen - wordLen + 1;

    #if GETWORDCOUNTSR_DEBUG
    cerr << "id=" << id
	 << "\nseq=" << seq
	 << "\nseqLen=" << seqLen
	 << endl;
    #endif

    for ( int j = 0; j < nWords; ++j )
    {
      char *word = seq + j;
      char c = word[wordLen];
      word[wordLen] = '\0';

      int hashIdx = hashFn(word, wordLen);

      #if GETWORDCOUNTSR_DEBUG
      cerr << "j=" << j
	   << "\thashIdx=" << hashIdx;
      #endif

      if ( hashIdx == -1 )
      {
	#if GETWORDCOUNTSR_DEBUG
	cerr << endl;
	#endif
	word[wordLen] = c;
	continue;
      }

      counts_m[modelIdx][ hashIdx ]++;

      #if GETWORDCOUNTSR_DEBUG
      cerr << "\tword=" << word
	   << "\tcounts_m[" << modelIdx << "][" << hashIdx << "]=" << counts_m[modelIdx][ hashIdx ]
	   << endl;
      #endif

      word[wordLen] = c;
    }

    free(id);
    free(seq);
  }

  int lL = hashUL_m[wordLen-1];
  int lU = hashUL_m[wordLen];

  if ( pseudoCountType_m == zeroOffset1 )
  {
    double pseudoCount = 1;
    for ( int i = lL; i < lU; ++i )
      counts_m[modelIdx][i] += pseudoCount;
  }
  else if ( pseudoCountType_m == zeroOffset0 )
  {
    #if 0
    cerr << "in MarkovChains2_t::wordCountsR pseudoCountType_m == zeroOffset0" << endl;
    for ( int i = lL; i < lU; ++i )
      cerr << counts_m[modelIdx][i] << endl;
    #endif
  }
  else if ( pseudoCountType_m == zeroOffset4mk )
  {
    double pseudoCount = 1.0/pow(4.0, wordLen) ;
    for ( int i = lL; i < lU; ++i )
      counts_m[modelIdx][i] += pseudoCount;
  }
  else if ( pseudoCountType_m == recPdoCount )
  {
    if ( wordLen==1 )
    {
      // for order 0 model

      // count[A] = count[A] + count[A]/total = count[A]*(1 + 1/total),
      // count[T] = count[T]* (1 + 1/total),
      // count[C] = count[C]* (1 + 1/total),
      // count[G] = count[G]* (1 + 1/total),

      double totalCount = 0;
      for ( int i = lL; i < lU; ++i )
	totalCount += counts_m[modelIdx][i];

      double pseudoCountFactor = 1.0 + 1.0/totalCount;
      for ( int i = lL; i < lU; ++i )
	counts_m[modelIdx][i] *= pseudoCountFactor;
    }
    else
    {
      // For order 1 model

      // count[AA] = count[AA] + prob(A)*prob(A)
      // count[AC] = count[AC] + prob(A)*prob(C)
      // count[AG] = count[AG] + prob(A)*prob(G)
      // count[AT] = count[AT] + prob(A)*prob(T)
      // etc

      // and for order 2 models

      // count[ACA] = count[ACA] + prob(AC) * prob(A)
      // count[ACC] = count[ACC] + prob(AC) * prob(C)
      // count[ACG] = count[ACG] + prob(AC) * prob(G)
      // count[ACT] = count[ACT] + prob(AC) * prob(T)
      // etc

      // where

      // prob(AC) and prob(A) are probabilities of observing AC and A after addition of
      // corresponding pseudo-counts.

      lL = hashUL_m[wordLen-2];
      lU = hashUL_m[wordLen-1];
      int nNucs  = 4;
      int w;

      double totalCount = 0;
      for ( int i = lL; i < lU; ++i )
	totalCount += counts_m[modelIdx][i];

      for ( int v = lL; v < lU; ++v )
      {
	for ( int i = 0; i < nNucs; ++i )
	{
	  w = tr_m[v][i];
	  //counts_m[modelIdx][w] += pow(10.0, log10cProb_m[modelIdx][v]) * pow(10.0, log10cProb_m[modelIdx][i]);
	  // log10cProb_m[modelIdx][v] is log10 of conditional probability not the probability of v
	  // here is a correct formula
	  counts_m[modelIdx][w] += counts_m[modelIdx][v]/totalCount * pow(10.0, log10cProb_m[modelIdx][i]);
	}
      }
    }
  }


  fclose(fp);
}


//------------------------------------------------------ wordCountsR ----
// version of wordCountsR that accepts seqID of a sequence to be excluded when
// computing word counts
void MarkovChains2_t::wordCountsR( int kIdx,
				   const char *file,
				   int modelIdx,
				   const char *seqID)
{
  #define GETWORDCOUNTSR_DEBUG 0

  static bool foundSeq = false;
  int wordLen = kIdx+1;
  char *id, *seq;
  int seqLen;

  FILE *fp = fOpen(file, "r");

  while ( getNextFastaRecord( fp, id, seq, seqLen) )
  {
    if ( !foundSeq && strcmp(id, seqID)==0 )
    {
      STRDUP(seq_m, seq);
      seqLen_m = seqLen;
      foundSeq = true;
      free(id);
      free(seq);
      //fprintf(stderr, "in wordCountsR() seq_m: %s\n", seq_m);
      //fprintf(stderr, "\nFound %s in wordCountsR()\tfile=%s\tmodelIdx=%d\n", seqID, file, modelIdx);
      continue;
    }

    int nWords = seqLen - wordLen + 1;

    #if GETWORDCOUNTSR_DEBUG
    cerr << "id=" << id
	 << "\nseq=" << seq
	 << "\nseqLen=" << seqLen
	 << endl;
    #endif

    for ( int j = 0; j < nWords; ++j )
    {
      char *word = seq + j;
      char c = word[wordLen];
      word[wordLen] = '\0';

      int hashIdx = hashFn(word, wordLen);

      #if GETWORDCOUNTSR_DEBUG
      cerr << "j=" << j
	   << "\thashIdx=" << hashIdx;
      #endif

      if ( hashIdx == -1 )
      {
	#if GETWORDCOUNTSR_DEBUG
	cerr << endl;
	#endif
	word[wordLen] = c;
	continue;
      }

      counts_m[modelIdx][ hashIdx ]++;

      #if GETWORDCOUNTSR_DEBUG
      cerr << "\tword=" << word
	   << "\tcounts_m[" << modelIdx << "][" << hashIdx << "]=" << counts_m[modelIdx][ hashIdx ]
	   << endl;
      #endif

      word[wordLen] = c;
    }

    free(id);
    free(seq);
  }

  int lL = hashUL_m[wordLen-1];
  int lU = hashUL_m[wordLen];

  if ( pseudoCountType_m == zeroOffset1 )
  {
    double pseudoCount = 1;
    for ( int i = lL; i < lU; ++i )
      counts_m[modelIdx][i] += pseudoCount;
  }
  else if ( pseudoCountType_m == zeroOffset0 )
  {
    #if 0
    cerr << "in MarkovChains2_t::wordCountsR pseudoCountType_m == zeroOffset0" << endl;
    for ( int i = lL; i < lU; ++i )
      cerr << counts_m[modelIdx][i] << endl;
    #endif
  }
  else if ( pseudoCountType_m == zeroOffset4mk )
  {
    double pseudoCount = 1.0/pow(4.0, wordLen) ;
    for ( int i = lL; i < lU; ++i )
      counts_m[modelIdx][i] += pseudoCount;
  }
  else if ( pseudoCountType_m == recPdoCount )
  {
    if ( wordLen==1 )
    {
      // for order 0 model

      // count[A] = count[A] + count[A]/total = count[A]*(1 + 1/total),
      // count[T] = count[T]* (1 + 1/total),
      // count[C] = count[C]* (1 + 1/total),
      // count[G] = count[G]* (1 + 1/total),

      double totalCount = 0;
      for ( int i = lL; i < lU; ++i )
	totalCount += counts_m[modelIdx][i];

      double pseudoCountFactor = 1.0 + 1.0/totalCount;
      for ( int i = lL; i < lU; ++i )
	counts_m[modelIdx][i] *= pseudoCountFactor;
    }
    else
    {
      // For order 1 model

      // count[AA] = count[AA] + prob(A)*prob(A)
      // count[AC] = count[AC] + prob(A)*prob(C)
      // count[AG] = count[AG] + prob(A)*prob(G)
      // count[AT] = count[AT] + prob(A)*prob(T)
      // etc

      // and for order 2 models

      // count[ACA] = count[ACA] + prob(AC) * prob(A)
      // count[ACC] = count[ACC] + prob(AC) * prob(C)
      // count[ACG] = count[ACG] + prob(AC) * prob(G)
      // count[ACT] = count[ACT] + prob(AC) * prob(T)
      // etc

      // where

      // prob(AC) and prob(A) are probabilities of observing AC and A after addition of
      // corresponding pseudo-counts.

      lL = hashUL_m[wordLen-2];
      lU = hashUL_m[wordLen-1];
      int nNucs  = 4;
      int w;

      double totalCount = 0;
      for ( int i = lL; i < lU; ++i )
	totalCount += counts_m[modelIdx][i];

      for ( int v = lL; v < lU; ++v )
      {
	for ( int i = 0; i < nNucs; ++i )
	{
	  w = tr_m[v][i];
	  //counts_m[modelIdx][w] += pow(10.0, log10cProb_m[modelIdx][v]) * pow(10.0, log10cProb_m[modelIdx][i]);
	  // log10cProb_m[modelIdx][v] is log10 of conditional probability not the probability of v
	  // here is a correct formula
	  counts_m[modelIdx][w] += counts_m[modelIdx][v]/totalCount * pow(10.0, log10cProb_m[modelIdx][i]);
	}
      }
    }
  }


  fclose(fp);
}

//------------------------------------------------- printWordCountsTbl ----
/// appends the row of modelId and tbl elements to file
/// if idx = 0, header of k-mer strings is added
void MarkovChains2_t::printWordCountsTbl( const char *file, int kIdx, int modelIdx )
{
  FILE *out;
  int wordLen = kIdx+1;
  int lL = hashUL_m[wordLen-1];// lower bound for hash val's of words of size l
  int lU = hashUL_m[wordLen];  // upper bound for hash val's of words of size l
  lU--;
  int i = 0;

  if ( modelIdx == 0 )
  {
    out = fOpen(file, "w");

    int nWords1 = wordStrgs_m[kIdx].size() - 1;

    fprintf(out,"modelId");
    for ( ; i < nWords1; ++i )
      fprintf(out,"\t%s",wordStrgs_m[kIdx][i]);
    fprintf(out,"\t%s\n",wordStrgs_m[kIdx][i]);
  }
  else
  {
    out = fOpen(file, "a");
  }

  fprintf(out,"%s",modelIds_m[modelIdx]);
  for ( i = lL; i < lU; ++i )
    fprintf(out,"\t%.1f", counts_m[modelIdx][i]);
  fprintf(out,"\t%.1f\n", counts_m[modelIdx][i]);

  fclose(out);
}


//------------------------------------------------ log10condProbs ----
/// computes log10 word frequencies
/// it is assumed that word counts of appropriate length are already computed
void MarkovChains2_t::log10condProbs( int kIdx, int modelIdx )
{
  int wordLen = kIdx+1;
  int lL = hashUL_m[wordLen-1];// lower bound for hash val's of words of size l
  int lU = hashUL_m[wordLen];  // upper bound for hash val's of words of size l
  double sum = 0;
  int nNucs  = 4;

  if ( kIdx == 0 )
  {
    for ( int i = lL; i < lU; ++i )
      sum += counts_m[modelIdx][i];

    double log10sum = log10(sum);
    for ( int i = lL; i < lU; ++i )
      log10cProb_m[modelIdx][i] = log10(counts_m[modelIdx][i]) - log10sum;
  }
  else
  {
    wordLen--;
    lL = hashUL_m[wordLen-1];
    lU = hashUL_m[wordLen];

    for ( int v = lL; v < lU; ++v )
    {
      sum = 0;
      for ( int i = 0; i < nNucs; ++i )
	sum += counts_m[modelIdx][tr_m[v][i]];

      sum = log10(sum);
      for ( int i = 0; i < nNucs; ++i )
      {
	int w = tr_m[v][i];
	log10cProb_m[modelIdx][w] = log10(counts_m[modelIdx][w]) - sum;
      }
    }
  }
}


//-------------------------------------------------- printLog10condProbTbl ----
/// appends the row of modelId and tbl elements to file
/// if idx = 0, header of k-mer strings is added
void MarkovChains2_t::printLog10condProbTbl( const char *file, int kIdx, int modelIdx )
{
  FILE *out;
  int wordLen = kIdx+1;
  int lL = hashUL_m[wordLen-1];// lower bound for hash val's of words of size l
  int lU = hashUL_m[wordLen];  // upper bound for hash val's of words of size l
  lU--;
  int i = 0;

  if ( modelIdx == 0 )
  {
    out = fOpen(file, "w");

    int nWords1 = wordStrgs_m[kIdx].size() - 1;

    fprintf(out,"modelId\t");
    for ( ; i < nWords1; ++i )
      fprintf(out,"%s\t",wordStrgs_m[kIdx][i]);
    fprintf(out,"%s\n",wordStrgs_m[kIdx][i]);
  }
  else
  {
    out = fOpen(file, "a");
  }

  fprintf(out,"%s\t",modelIds_m[modelIdx]);
  for ( i = lL; i < lU; ++i )
    fprintf(out,"%.15lf\t",log10cProb_m[modelIdx][i]);
  fprintf(out,"%.15lf\n",log10cProb_m[modelIdx][i]);

  fclose(out);
}

//--------------------------------------------------------------- getKmers ----
/// initializes words_m to contain all k-mers of size k
/// and creates hash table of word indices
void MarkovChains2_t::getKmers( int k )
{
  if ( (k < 1) || (k > order_m+1) )
  {
    cerr << "Error in " << __FILE__
	 << " MarkovChains2_t::getKmers() at line "
	 << __LINE__ << " arg out of range; k=" << k << endl;
    exit(1);
  }

  getAllKmers(k, wordStrgs_m[k-1]);
}

//------------------------------------------------------ hashFnIUPAC ----
/// ACGT string hash function accepting IUPAC codes
void MarkovChains2_t::hashFnIUPAC(const char *s, int sLen, vector<int> &idxs)
{
  idxs.clear();
  vector<int> codes;
  int codeSize;
  int h = hashUL_m[sLen-1];

  if ( intACGTLookup[int(*s)] > -1 )
  {
    h += intACGTLookup[int(*s)];
    idxs.push_back(h);
  }
  else if ( (codeSize=getIUPACambCodeHashVals(*s, codes)) )
  {
    for ( int i = 0; i < codeSize; ++i )
      idxs.push_back(h + codes[i]);
  }
  else
  {
    idxs.clear();
    cerr << "Error in MarkovChains2_t::hashFnIUPAC(): at line "
	 << __LINE__ << ": Unrecognized character in " << s << endl;
    return;
  }

  s++;
  int i;

  for ( ; *s; ++s )
  {
    int n = idxs.size();

    if ( (i=intACGTLookup[int(*s)]) > -1 )
    {
      for ( int j = 0; j < n; ++j )
	idxs[j] = tr_m[idxs[j]][i];
    }
    else if ( (codeSize=getIUPACambCodeHashVals(*s, codes)) )
    {
      for ( int j = 0; j < n; ++j )
	idxs[j] = tr_m[idxs[j]][codes[0]];

      for ( i = 1; i < codeSize; ++i )
      {
	for ( int j = 0; j < n; ++j )
	  idxs.push_back( tr_m[idxs[j]][codes[i]] );
      }
    }
    else
    {
      idxs.clear();
      cerr << "Error in MarkovChains2_t::hashFnIUPAC(): at line "
	   << __LINE__ << ": Unrecognized character in " << s << endl;
      break;
    }
  }
}

//------------------------------------------------------ sample ----
/// generating 'sampleSize' random samples from each MC model
void MarkovChains2_t::sample( const char *faFile, const char *txFile, int sampleSize, int seqLen )
{
  int nNucs    = 4;
  int nModels  = modelIds_m.size();
  double **cProb;

  char *dnaStr;
  MALLOC(dnaStr, char*, (seqLen+1) * sizeof(char));
  dnaStr[seqLen] = '\0';

  //-- setting up cProb table so we don't have to do pow(10.0, - ) during sampling
  MALLOC(cProb, double**, nModels * sizeof(double*));
  for ( int i = 0; i < nModels; ++i )
  {
    MALLOC(cProb[i], double*, nAllWords_m * sizeof(double));
    for ( int j = 0; j < nAllWords_m; ++j )
      cProb[i][j] = pow(10.0, log10cProb_m[i][j]);
  }

  int intMers[4];
  int intB; // hash value of chosed base
  int v;    // k-mer hash int value
  double u; // random number in [0,1]

  vector<char *> nucs;
  getAllKmers( 1, nucs );

  // initialize random seed so that consecutive calls of rand() do not generate similar numbers
  srand( rand() );
  srand( rand() );

  FILE *faOut = fOpen( faFile, "w");
  FILE *txOut = fOpen( txFile, "w");

  for ( int modelIdx = 0; modelIdx < nModels; modelIdx++ )
  {
    for ( int s = 0; s < sampleSize; s++ )
    {
      fprintf(faOut,">%s_%d\n",modelIds_m[modelIdx],s);
      fprintf(txOut,"%s_%d\t%s\n",modelIds_m[modelIdx],s,modelIds_m[modelIdx]);

      u = (double)rand() / RAND_MAX;

      // 1-mers
      if ( u < cProb[modelIdx][0] )
	intB = 0;
      else if ( u < cProb[modelIdx][0] + cProb[modelIdx][1] )
	intB = 1;
      else if ( u < cProb[modelIdx][0] + cProb[modelIdx][1] + cProb[modelIdx][2] )
	intB = 2;
      else
	intB = 3;

#if 0
      cerr << "cummProbs: " << cProb[modelIdx][0] << ", "
	   <<  (cProb[modelIdx][0] + cProb[modelIdx][1]) << ", "
	   <<  (cProb[modelIdx][0] + cProb[modelIdx][1] + cProb[modelIdx][2])
	   << endl;
      cerr << "\nu=" << u << "\tintB=" << intB << "\tnuc=" << nucs[intB] << endl;
#endif

      v = intB;
      dnaStr[0] = nucs[intB][0];

      // k-mers; k>1
      for( int i = 1; i < seqLen; ++i )
      {
	for ( int j = 0; j < nNucs; j++ )
	  intMers[j] = tr_m[v][j];

	u = (double)rand() / RAND_MAX;

	if ( u < cProb[modelIdx][intMers[0]] )
	  intB = 0;
	else if ( u < cProb[modelIdx][intMers[0]] + cProb[modelIdx][intMers[1]] )
	  intB = 1;
	else if ( u < cProb[modelIdx][intMers[0]] + cProb[modelIdx][intMers[1]] + cProb[modelIdx][intMers[2]] )
	  intB = 2;
	else
	  intB = 3;

	dnaStr[i] = nucs[intB][0];

#if 0
	cerr << "\ncummProbs: " << cProb[modelIdx][intMers[0]] << ", "
	     <<  (cProb[modelIdx][intMers[0]] + cProb[modelIdx][intMers[1]]) << ", "
	     <<  (cProb[modelIdx][intMers[0]] + cProb[modelIdx][intMers[1]] + cProb[modelIdx][intMers[2]])
	     << endl;
	cerr << "u=" << u << "\tintB=" << intB << "\tnuc=" << nucs[intB] << endl;
#endif
	v = intMers[intB];
      }

      fprintf(faOut,"%s\n", dnaStr);
    }
  }

  free(dnaStr);
  for ( int i = 0; i < nModels; ++i )
    free(cProb[i]);
  free(cProb);

  fclose(faOut);
  fclose(txOut);
}


//------------------------------------------------------ sample ----
/// generating 'sampleSize' random samples from an MC model with a given model index
void MarkovChains2_t::sample( char ***_seqTbl, int modelIdx, int sampleSize, int seqLen )
{
  char **seqTbl;
  MALLOC(seqTbl, char**, sampleSize * sizeof(char*));
  for ( int i = 0; i < sampleSize; ++i )
  {
    MALLOC(seqTbl[i], char*, (seqLen+1) * sizeof(char));
    seqTbl[i][seqLen] = '\0';
  }

  //-- setting up cProb table so we don't have to do pow(10.0, - ) during sampling
  double *cProb;
  MALLOC(cProb, double*, nAllWords_m * sizeof(double));
  for ( int j = 0; j < nAllWords_m; ++j )
      cProb[j] = pow(10.0, log10cProb_m[modelIdx][j]);

  #if 0
  for( int l = 1; l < maxWordLen; ++l )
  {
    lL = hashUL_m[l-1]; // lower limit for hash values of k-mers of length l
    lU = hashUL_m[l];   // upper limit for hash values of k-mers of length l
    for ( int v = lL; v < lU; ++v )
    {
      for ( int i = 0; i < nNucs; ++i )
      {
	tr_m[v][i] = v + p4*(1 + i);
      }
    }
  }
  #endif

  int nNucs    = 4;
  int intMers[4];
  int intB; // hash value of chosed base
  int v;    // k-mer hash int value
  double u; // random number in [0,1]

  vector<char *> nucs;
  getAllKmers( 1, nucs );

  // initialize random seed so that consecutive calls of rand() do not generate similar numbers
  srand( rand() );
  srand( rand() );

  for ( int s = 0; s < sampleSize; s++ )
  {
    u = (double)rand() / RAND_MAX;

    // 1-mers
    if ( u < cProb[0] )
      intB = 0;
    else if ( u < cProb[0] + cProb[1] )
      intB = 1;
    else if ( u < cProb[0] + cProb[1] + cProb[2] )
      intB = 2;
    else
      intB = 3;

#if 0
    cerr << "cummProbs: " << cProb[0] << ", "
	 <<  (cProb[0] + cProb[1]) << ", "
	 <<  (cProb[0] + cProb[1] + cProb[2])
	 << endl;
    cerr << "\nu=" << u << "\tintB=" << intB << "\tnuc=" << nucs[intB] << endl;
#endif

    v = intB;
    seqTbl[s][0] = nucs[intB][0];

    // k-mers; k>1
    for( int i = 1; i < seqLen; ++i )
    {
      for ( int j = 0; j < nNucs; j++ )
	intMers[j] = tr_m[v][j];

      srand( rand() );
      srand( rand() );
      u = (double)rand() / RAND_MAX;

      if ( u < cProb[intMers[0]] )
	intB = 0;
      else if ( u < cProb[intMers[0]] + cProb[intMers[1]] )
	intB = 1;
      else if ( u < cProb[intMers[0]] + cProb[intMers[1]] + cProb[intMers[2]] )
	intB = 2;
      else
	intB = 3;

      seqTbl[s][i] = nucs[intB][0];

#if 0
	cerr << "\ncummProbs: " << cProb[intMers[0]] << ", "
	     <<  (cProb[intMers[0]] + cProb[intMers[1]]) << ", "
	     <<  (cProb[intMers[0]] + cProb[intMers[1]] + cProb[intMers[2]])
	     << endl;
	cerr << "u=" << u << "\tintB=" << intB << "\tnuc=" << nucs[intB] << endl;
#endif
	v = intMers[intB];
    }
  }

  *_seqTbl = seqTbl;

  free(cProb);
}



//------------------------------------------------------ sampleMF ----
/// generating 'sampleSize' random samples from an MC model with a given model index
/// assuming the user allocated memory for seqTbl
void MarkovChains2_t::sampleMF( char **seqTbl, int modelIdx, int sampleSize, int seqLen )
{
  //-- setting up cProb table so we don't have to do pow(10.0, - ) during sampling
  // double *cProb;
  // MALLOC(cProb, double*, nAllWords_m * sizeof(double));
  for ( int j = 0; j < nAllWords_m; ++j )
      cProb_m[j] = pow(10.0, log10cProb_m[modelIdx][j]);

  #if 0
  for( int l = 1; l < maxWordLen; ++l )
  {
    lL = hashUL_m[l-1]; // lower limit for hash values of k-mers of length l
    lU = hashUL_m[l];   // upper limit for hash values of k-mers of length l
    for ( int v = lL; v < lU; ++v )
    {
      for ( int i = 0; i < nNucs; ++i )
      {
	tr_m[v][i] = v + p4*(1 + i);
      }
    }
  }
  #endif

  int nNucs    = 4;
  int intMers[4];
  int intB; // hash value of chosed base
  int v;    // k-mer hash int value
  double u; // random number in [0,1]

  // vector<char *> nucs;
  // getAllKmers( 1, nucs );

  // initialize random seed so that consecutive calls of rand() do not generate similar numbers
  srand( rand() );
  srand( rand() );

  for ( int s = 0; s < sampleSize; s++ )
  {
    u = (double)rand() / RAND_MAX;

    // 1-mers
    if ( u < cProb_m[0] )
      intB = 0;
    else if ( u < cProb_m[0] + cProb_m[1] )
      intB = 1;
    else if ( u < cProb_m[0] + cProb_m[1] + cProb_m[2] )
      intB = 2;
    else
      intB = 3;

#if 0
    cerr << "cummProbs: " << cProb_m[0] << ", "
	 <<  (cProb_m[0] + cProb_m[1]) << ", "
	 <<  (cProb_m[0] + cProb_m[1] + cProb_m[2])
	 << endl;
    cerr << "\nu=" << u << "\tintB=" << intB << "\tnuc=" << nucs_m[intB] << endl;
#endif

    v = intB;
    seqTbl[s][0] = nucs_m[intB][0];

    // k-mers; k>1
    for( int i = 1; i < seqLen; ++i )
    {
      for ( int j = 0; j < nNucs; j++ )
	intMers[j] = tr_m[v][j];

      srand( rand() );
      srand( rand() );
      u = (double)rand() / RAND_MAX;

      if ( u < cProb_m[intMers[0]] )
	intB = 0;
      else if ( u < cProb_m[intMers[0]] + cProb_m[intMers[1]] )
	intB = 1;
      else if ( u < cProb_m[intMers[0]] + cProb_m[intMers[1]] + cProb_m[intMers[2]] )
	intB = 2;
      else
	intB = 3;

      seqTbl[s][i] = nucs_m[intB][0];

#if 0
	cerr << "\ncummProbs: " << cProb_m[intMers[0]] << ", "
	     <<  (cProb_m[intMers[0]] + cProb_m[intMers[1]]) << ", "
	     <<  (cProb_m[intMers[0]] + cProb_m[intMers[1]] + cProb_m[intMers[2]])
	     << endl;
	cerr << "u=" << u << "\tintB=" << intB << "\tnuc=" << nucs_m[intB] << endl;
#endif
	v = intMers[intB];
    }
  }
}


//------------------------------------------------------ sample ----
/// generating 'sampleSize' random samples from an MC model and a set of reference sequences
/// random sequences are generated by
/// 1. selecting at random a ref sequence
/// 2. using high rank MC model to modify bases of the sequence
void MarkovChains2_t::sample( char ***_seqTbl, map<string, string> &refSeqs, int modelIdx, int sampleSize, int seqLen )
{
  char **seqTbl;
  MALLOC(seqTbl, char**, sampleSize * sizeof(char*));
  for ( int i = 0; i < sampleSize; ++i )
  {
    MALLOC(seqTbl[i], char*, (seqLen+1) * sizeof(char));
    seqTbl[i][seqLen] = '\0';
  }

  //-- setting up cProb table so we don't have to do pow(10.0, - ) during sampling
  double *cProb;
  MALLOC(cProb, double*, nAllWords_m * sizeof(double));
  for ( int j = 0; j < nAllWords_m; ++j )
      cProb[j] = pow(10.0, log10cProb_m[modelIdx][j]);

  #if 0
  for( int l = 1; l < maxWordLen; ++l )
  {
    lL = hashUL_m[l-1]; // lower limit for hash values of k-mers of length l
    lU = hashUL_m[l];   // upper limit for hash values of k-mers of length l
    for ( int v = lL; v < lU; ++v )
    {
      for ( int i = 0; i < nNucs; ++i )
      {
	tr_m[v][i] = v + p4*(1 + i);
      }
    }
  }
  #endif

  int nNucs    = 4;
  int intMers[4];
  int intB; // hash value of chosed base
  double u; // random number in [0,1]

  vector<char *> nucs;
  getAllKmers( 1, nucs );

  // initialize random seed so that consecutive calls of rand() do not generate similar numbers
  srand( rand() );
  srand( rand() );

  int nseqs = refSeqs.size();
  int nseqs1 = nseqs - 1;

  // extract IDs of reference sequences
  vector<string> ids;
  map<string, string>::iterator itr;
  //map<string, int> rCount;
  for ( itr=refSeqs.begin(); itr != refSeqs.end(); ++itr )
  {
    ids.push_back( itr->first );
    //rCount[ itr->first ] = 0;
  }

  int wordLen = order_m + 1;
  int nWords  = seqLen - wordLen + 1;

  for ( int s = 0; s < sampleSize; s++ )
  {
    // choose a ref sequence
    u = (double)rand() / RAND_MAX;
    int iseq = (int)(u * nseqs1);
    const char *refSeq = refSeqs[ ids[iseq] ].c_str();

    //rCount[ ids[iseq] ]++;

    for ( int i = 0; i < wordLen; ++i )
      seqTbl[s][i] = refSeq[i];

    for ( int j = 0; j < nWords; ++j )
    {
      char *word = seqTbl[s] + j;
      char c = word[wordLen];
      word[wordLen] = '\0';
      int v = hashFn(word, wordLen); // k-mer hash int value
      word[wordLen] = c;

      for ( int k = 0; k < nNucs; k++ )
	intMers[k] = tr_m[v][k];

      srand( rand() );
      srand( rand() );
      u = (double)rand() / RAND_MAX;

      if ( u < cProb[intMers[0]] )
	intB = 0;
      else if ( u < cProb[intMers[0]] + cProb[intMers[1]] )
	intB = 1;
      else if ( u < cProb[intMers[0]] + cProb[intMers[1]] + cProb[intMers[2]] )
	intB = 2;
      else
	intB = 3;

      seqTbl[s][ wordLen + j ] = nucs[intB][0];

#if 0
	cerr << "\ncummProbs: " << cProb[intMers[0]] << ", "
	     <<  (cProb[intMers[0]] + cProb[intMers[1]]) << ", "
	     <<  (cProb[intMers[0]] + cProb[intMers[1]] + cProb[intMers[2]])
	     << endl;
	cerr << "u=" << u << "\tintB=" << intB << "\tnuc=" << nucs[intB] << endl;
#endif
	//v = intMers[intB];
    }
  }

  #if 0
  map<string, int>::iterator itr2;
  cerr << "rCount\n";
  for ( itr2 = rCount.begin(); itr2 != rCount.end(); ++itr2 )
    cerr << itr2->first << "\t" << rCount[ itr2->first ] << endl;
  cerr << endl;
  exit(1);
  #endif

  *_seqTbl = seqTbl;
}
