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

#include "MarkovChains.hh"
#include "CStatUtilities.h"
#include "CUtilities.h"
#include "IOCppUtilities.hh"
#include "IOCUtilities.h"
#include "CppStatUtilities.hh"
#include "CppUtilities.hh"
#include "strings.hh"
#include "DNAsequence.hh"
#include "Newick.hh"

using namespace std;

//--------------------------------------------------------- MarkovChains_t ----
/// Attempts to read already created models from model files
MarkovChains_t::MarkovChains_t(int order,
                               char *dir,
                               int maxNumAmbCodes,
                               int pseudoCountType,
                               bool verbose)
  : order_m(order), dir_m(dir), maxNumAmbCodes_m(maxNumAmbCodes), pseudoCountType_m(pseudoCountType)
{
    #define DEBUG_MARKOVCHAINS1 0

    CHECK_DIR( dir_m );

    int maxWordLen = order_m+1;

    // ------------------------------------------------
    if ( verbose )
      fprintf(stderr,"\t\tCreating all k-mers ... ");

    getAllKmers( 1, nucs_m );

    if ( verbose )
      fprintf(stderr,"DONE\n");

    // ------------------------------------------------
    // wordStrgs_m is needed only for the header of cProbFile_m's
    wordStrgs_m.resize(maxWordLen);
    for ( int k = 1; k <= maxWordLen; ++k )
      getAllKmers(k, wordStrgs_m[k-1]);

    #if DEBUG_MARKOVCHAINS1
    fprintf(stderr, "\n\nIn MarkovChains_t::MarkovChains_t(order, dir, maxNumAmbCodes, pseudoCountType) order_m=%d  maxWordLen=%d\n",
            order_m, maxWordLen);
    //cerr << "wordStrgs_m[0]:" << endl;
    //printVector(wordStrgs_m[0]);
    //cerr << "wordStrgs_m[1]:" << endl;
    //printVector(wordStrgs_m[1]);
    #endif

    setupIOfiles();

    // ------------------------------------------------
    if ( verbose )
      fprintf(stderr,"\t\tCreating model IDs ... ");

    readModelIds();

    if ( verbose )
      fprintf(stderr,"DONE\n");


    // ------------------------------------------------
    if ( verbose )
      fprintf(stderr,"\t\tCreating Moore Machines ... ");

    createMooreMachine();

    if ( verbose )
      fprintf(stderr,"DONE\n");


    MALLOC(cProb_m, double*, nAllWords_m * sizeof(double));

    // ------------------------------------------------
    if ( verbose )
      fprintf(stderr,"\t\tInitializing IUPAC hash values ... ");

    initIUPACambCodeHashVals();

    if ( verbose )
      fprintf(stderr,"DONE\n");


    // ------------------------------------------------
    if ( verbose )
      fprintf(stderr,"\t\tReading MC models of oder 0 ... ");

    bool ok = readLog10condProbTbl(0, cProbFile_m[0].c_str());
    if ( !ok ) {
      fprintf(stderr, "ERROR: in file %s at line %d: cannot read order 0 MC models table",__FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }

    if ( verbose )
      fprintf(stderr,"DONE\n");


    for ( int i = 1; i < maxWordLen; ++i )
    {
      if ( verbose )
        fprintf(stderr,"\t\tReading MC models of oder %d ... ",i);

      ok = readLog10condProbTbl(i, cProbFile_m[i].c_str());
      if ( !ok ) {
        fprintf(stderr, "ERROR: in file %s at line %d: cannot read order %d MC models table",
                __FILE__, __LINE__, i);
        exit(EXIT_FAILURE);
      }

      if ( verbose )
        fprintf(stderr,"DONE\n");
    }
}

//--------------------------------------------------------- MarkovChains_t ----
/// Creates models from traing fasta files overwriting existing model files
MarkovChains_t::MarkovChains_t(int order,
                               vector<char *> &trgFiles,
                               char *dir,
                               int maxNumAmbCodes,
                               int pseudoCountType,
                               bool verbose)
  : order_m(order), dir_m(dir), maxNumAmbCodes_m(maxNumAmbCodes), pseudoCountType_m(pseudoCountType)
{
    #define DEBUG_MARKOVCHAINS2 0

    CHECK_DIR( dir_m );

    int maxWordLen = order_m+1;

    // ------------------------------------------------
    if ( verbose )
      fprintf(stderr,"\t\tCreating all k-mers ... ");

    getAllKmers( 1, nucs_m );

    if ( verbose )
      fprintf(stderr,"DONE\n");


    // ------------------------------------------------
    // wordStrgs_m is needed only for the header of cProbFile_m's
    wordStrgs_m.resize(maxWordLen);
    for ( int k = 1; k <= maxWordLen; ++k )
      getAllKmers(k, wordStrgs_m[k-1]);

    int nTrgFiles = trgFiles.size();
    char *file;
    for ( int i = 0; i < nTrgFiles; ++i )
    {
      STRDUP(file, trgFiles[i]);
      trgFiles_m.push_back(file);
    }

    #if DEBUG_MARKOVCHAINS2
    fprintf(stderr, "In MarkovChains_t::MarkovChains_t(order, trgFiles, ... ): order_m=%d  maxWordLen=%d  nTrgFiles=%d\n",
            order_m, maxWordLen, nTrgFiles);
    #endif

    setupIOfiles();

    // ------------------------------------------------
    if ( verbose )
      fprintf(stderr,"\t\tCreating model IDs ... ");

    createModelIds();

    if ( verbose )
      fprintf(stderr,"DONE\n");


    // ------------------------------------------------
    if ( verbose )
      fprintf(stderr,"\t\tCreating Moore Machines ... ");

    createMooreMachine();

    if ( verbose )
      fprintf(stderr,"DONE\n");


    MALLOC(cProb_m, double*, nAllWords_m * sizeof(double));


    // ------------------------------------------------
    if ( verbose )
      fprintf(stderr,"\t\tInitializing IUPAC hash values ... ");

    initIUPACambCodeHashVals();

    if ( verbose )
      fprintf(stderr,"DONE\n");


    // ------------------------------------------------
    for ( int i = 0; i < maxWordLen; ++i )
    {
      if ( verbose )
        fprintf(stderr,"\t\tCreating MC models of oder %d ... ",i);
      log10condProbs( i );
      if ( verbose )
        fprintf(stderr,"DONE\n");
    }
}


//--------------------------------------------------------- MarkovChains_t ----
/// MarkovChains_t version for leave-one-out validation setup (dir is NULL)
MarkovChains_t::MarkovChains_t(int order,
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
bool MarkovChains_t::readModelIds()
{
    CHECK_DIR( dir_m );

    string inFile(dir_m);
    inFile += "/modelIds.txt";
    CPP_CHECK_FILE( inFile );

    readLines(inFile.c_str(), modelIds_m);

    return true;
}

//-------------------------------------------------------- createModelIds ----
void MarkovChains_t::createModelIds()
{
    #define DEBUG_CREATEMODELIDS 0

    CHECK_DIR( dir_m );

    int n = (int)trgFiles_m.size();
    if ( n == 0 )
    {
      fprintf(stderr, "ERROR: in file %s at line %d: trgFiles_m.size()=0",__FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }

    #if DEBUG_CREATEMODELIDS
    fprintf(stderr, "In MarkovChains_t::createModelIds(): trgFiles_m.size()=%d", n);
    #endif

    char *idStr;
    for ( int i = 0; i < n; ++i )
    {
      string id = baseFileName(trgFiles_m[i]);
      STRDUP(idStr,id.c_str());
      modelIds_m.push_back(idStr);
    }

    string outFile(dir_m);
    outFile += "/modelIds.txt";
    FILE *out = fOpen(outFile.c_str(), "w");
    for ( int i = 0; i < n; ++i )
      fprintf(out,"%s\n",modelIds_m[i]);
    fclose(out);

    #if DEBUG_CREATEMODELIDS
    fprintf(stderr, "modelIds_m\n");
    for ( int i = 0; i < n; ++i )
      fprintf(stderr, "%s\n", modelIds_m[i]);
    fprintf(stderr, "\n");
    #endif
}

//-------------------------------------------------------- ~MarkovChains_t ----
MarkovChains_t::~MarkovChains_t()
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
void MarkovChains_t::setupIOfiles()
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
void MarkovChains_t::createMooreMachine( )
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
bool MarkovChains_t::isIUPACambCode(char c) const
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
        cerr << "Error in MarkovChains_t::isIUPACambCode(): Unrecognized IUPAC character: "
             << (char)c << endl;
        return false;
    }
}

//------------------------------------------------ initIUPACambCodeHashVals -----------
/// initialize IUPAC ambiguity codes hash value vectors
void MarkovChains_t::initIUPACambCodeHashVals()
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
int MarkovChains_t::getIUPACambCodeHashVals(char c, vector<int> &v) const
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
        cerr << "Error in MarkovChains_t::getIUPACambCodeHashVals(): Unrecognized IUPAC character: "
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
double MarkovChains_t::log10probIUPAC( const char *frag, int fragLen, int modelIdx )
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
double MarkovChains_t::log10prob( const char *frag, int fragLen, int modelIdx )
{
#define DEBUG_LOG10PROB 0

    double log10probVal = 0;
    int k = 0, v, i;

#if DEBUG_LOG10PROB
    fprintf(stderr, "log10prob(): Before init if\n");
#endif

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

#if DEBUG_LOG10PROB
    fprintf(stderr, "log10prob(): After init if v=%d k=%d rank=%d\n", v, k, rank);
#endif

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

#if DEBUG_LOG10PROB
    fprintf(stderr, "log10prob(): After init while v=%d k=%d log10probVal=%.10lf\n", v, k, log10probVal);
#endif


#if DEBUG_LOG10PROB
    modelIdx = 1375;
    for ( int i = 0; i < 100; i++ )
      fprintf(stderr, "log10cProb_m[modelIdx][%d]=%.10lf\n", i, log10cProb_m[modelIdx][i]);
#endif

    //-- process the remaining k-mers
    while ( k < fragLen )
    {
      if ( (i=intACGTLookup[int(frag[k])]) > -1 )
      {
        v = tr_m[v][i];
        log10probVal += log10cProb_m[modelIdx][v];

#if DEBUG_LOG10PROB
        fprintf(stderr, "k=%d v=%d modelIdx=%d log10probVal=%.10lf\n", k, v, modelIdx, log10probVal);
        exit(0);
#endif
      }
      else
      {
        return log10probIUPAC(frag, fragLen, modelIdx);
      }
      k++;
    }

#if DEBUG_LOG10PROB
    fprintf(stderr, "FINAL k=%d log10probVal=%.10lf\n", k, log10probVal);
    exit(0);
#endif

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
int MarkovChains_t::log10probVect( const char *frag, int fragLen, int modelIdx, double *probs )
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
double MarkovChains_t::log10probR( char *frag, int fragLen, int modelIdx )
{
#define LOG10CPROBR_DEBUG 0
#if LOG10CPROBR_DEBUG
    cerr << "in MarkovChains_t::log10prob()\tmodelIdx=" << modelIdx
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
bool MarkovChains_t::readLog10condProbTbl( int kIdx, const char *inFile )
{
    #define DEBUG_RLCPT 0

    #if DEBUG_RLCPT
    fprintf(stderr, "Entering readLog10condProbTbl(%d, %s)\n", kIdx, inFile);
    #endif

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
      fprintf(stderr, "ERROR: MarkovChains_t::readLog10condProbTbl() in file %s at line %d: Could not read the data\n", __FILE__, __LINE__);
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
          fprintf(stderr, "ERROR: MarkovChains_t::readLog10condProbTbl() in file %s at line %d: i == nWords; i=%d\tnWords=%d\n",
                  __FILE__, __LINE__, i, (lU-lL));
          return false;
        }

        log10cProb_m[modelIdx][i] = strtod(word, (char **)NULL);
      }
      modelIdx++;
    }

    #if DEBUG_RLCPT
    fprintf(stderr, "Done parsing %s\n", inFile);
    #endif

    fclose(fp);
    free(data);

    return true;
}


//------------------------------------------------- log10condProbs ----
void MarkovChains_t::log10condProbs( int kIdx )
{
    int nModels = modelIds_m.size();
    CHECK_DIR( dir_m );

    for ( int i = 0; i < nModels; ++i )
    {
      //cerr << "Generating " << (kIdx+1) << "-mer FreqTbl for " << trgFiles_m[i] << endl;;
      //-- compute word counts and rel. frequencies
      wordCountsR( kIdx, trgFiles_m[i], i );

      if ( printCounts_m )
        printWordCountsTbl( countFile_m[kIdx].c_str(), kIdx, i );

      log10condProbs( kIdx, i );
      printLog10condProbTbl( cProbFile_m[kIdx].c_str(), kIdx, i );
    }
}

//------------------------------------------------- log10condProbs ----
void MarkovChains_t::log10condProbs( int kIdx, const char *seqID )
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
void MarkovChains_t::wordCounts( int kIdx,
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
void MarkovChains_t::wordCountsR( int kIdx,
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
      cerr << "in MarkovChains_t::wordCountsR pseudoCountType_m == zeroOffset0" << endl;
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
void MarkovChains_t::wordCountsR( int kIdx,
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
      cerr << "in MarkovChains_t::wordCountsR pseudoCountType_m == zeroOffset0" << endl;
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
void MarkovChains_t::printWordCountsTbl( const char *file, int kIdx, int modelIdx )
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
void MarkovChains_t::log10condProbs( int kIdx, int modelIdx )
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
void MarkovChains_t::printLog10condProbTbl( const char *file, int kIdx, int modelIdx )
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
void MarkovChains_t::getKmers( int k )
{
    if ( (k < 1) || (k > order_m+1) )
    {
      cerr << "Error in " << __FILE__
           << " MarkovChains_t::getKmers() at line "
           << __LINE__ << " arg out of range; k=" << k << endl;
      exit(1);
    }

    getAllKmers(k, wordStrgs_m[k-1]);
}

//------------------------------------------------------ hashFnIUPAC ----
/// ACGT string hash function accepting IUPAC codes
void MarkovChains_t::hashFnIUPAC(const char *s, int sLen, vector<int> &idxs)
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
      cerr << "Error in MarkovChains_t::hashFnIUPAC(): at line "
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
        cerr << "Error in MarkovChains_t::hashFnIUPAC(): at line "
             << __LINE__ << ": Unrecognized character in " << s << endl;
        break;
      }
    }
}

//------------------------------------------------------ sample ----
/// generating 'sampleSize' random samples from each MC model
void MarkovChains_t::sample( const char *faFile, const char *txFile, int sampleSize, int seqLen )
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
void MarkovChains_t::sample( char ***_seqTbl, int modelIdx, int sampleSize, int seqLen )
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
/// Generates 'sampleSize' random samples from an MC model with a given model index
/// assuming the user allocated memory for seqTbl
void MarkovChains_t::sampleMF( char **seqTbl, int modelIdx, int sampleSize, int seqLen )
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

//------------------------------------------------------ sample_cprobs ----
/*!
    Generates 'sampleSize' random samples from an MC model with a given model index
    assuming the user allocated memory for seqTbl

    \param modelIdx   - A model index.
    \param sampleSize - The number of random sequences generated by the model.
    \param seqLen     - The length of random sequences generated by the model.
    \param cprobTbl   - The output conditional probabilities table (memory allocation is on the user).
*/
void MarkovChains_t::sample_cprobs(int modelIdx, int sampleSize, int seqLen, double **cprobTbl)
{
    //-- setting up cProb table so we don't have to do pow(10.0, - ) during sampling
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

    //char *seq;
    //MALLOC(seq, char*, seqLen * sizeof(char));

    int nNucs = 4;
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
      //seq[0] = nucs_m[intB][0];

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

        //seq[i] = nucs_m[intB][0];

        #if 0
        cerr << "\ncummProbs: " << cProb_m[intMers[0]] << ", "
             <<  (cProb_m[intMers[0]] + cProb_m[intMers[1]]) << ", "
             <<  (cProb_m[intMers[0]] + cProb_m[intMers[1]] + cProb_m[intMers[2]])
             << endl;
        cerr << "u=" << u << "\tintB=" << intB << "\tnuc=" << nucs_m[intB] << endl;
        #endif

        v = intMers[intB];

        cprobTbl[s][i] = log10cProb_m[modelIdx][v];

      } // END OF for( int i = 1; i < seqLen; ++i )

    } // END OF for ( int s = 0; s < sampleSize; s++ )

    //free(seq);
}


//------------------------------------------------------ sample ----
/// generating 'sampleSize' random samples from an MC model and a set of reference sequences
/// random sequences are generated by
/// 1. selecting at random a ref sequence
/// 2. using high rank MC model to modify bases of the sequence
void MarkovChains_t::sample( char ***_seqTbl, map<string, string> &refSeqs, int modelIdx, int sampleSize, int seqLen )
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

    free(cProb);
}

//------------------------------------------------------ sample_pp ----
/// Generates log10 pp's of random sequences from the given model and its siblings
/// and saves it to a file in the llpDir.
///
/// Random sequences are generated using corresponding models.
/// lpp's of all sequences are computed using the given node model.
///
/// node  - A node of the model tree.
/// ssize - The number of random sequences/llp's to be generated.
/// mcDir - A directory holding MC models, fasta files etc, within which rLpps output directory will be created.
/// seqLen - The sequence length of the random sequences
void MarkovChains_t::sample_pp( NewickNode_t *node, int ssize, char *mcDir, int seqLen )
{
    #define DEBUG_SAMPLE_PP 0

    string lppDir = string(mcDir) + string("/rLpps/");
    if ( !dir_exists( lppDir.c_str() ) ) { // checking for existence of lppDir
      string cmd("mkdir -p ");
      cmd += lppDir;
      system(cmd.c_str());
    }

    string lppFile = lppDir + node->label + string(".csv");
    FILE *lppFH = fOpen( lppFile.c_str(), "w");

    // format of the output file
    // <sp name>,<lpp1>,<lpp2>, ... ,<lppn>
    // <sib1 name>,<lpp1>,<lpp2>, ... ,<lppn>
    // <sib2 name>,<lpp1>,<lpp2>, ... ,<lppn>
    // ...

    // llpi's are lpp's of the corresponding random sequences using the 'node'
    // model


    // Setting random sequence length to the maximum of the reference sequnces
    // of the node model. For this lengths of the reference sequences are
    // examined to find their maximum.

    // Reading fasta file of the reference sequences of the current node (model)
    #if 0
    string faFile = string(mcDir) + string("/fasta_files/") + node->label + string(".fa");
    map<string, string> seqRecs; // fasta file sequence records
    seqRecs.clear();
    readFasta( faFile.c_str(), seqRecs);

    int maxSeqLen = 0;
    int s;
    map<string, string>::iterator itr;
    for ( itr = seqRecs.begin(); itr != seqRecs.end(); ++itr )
    {
      s = (int)itr->second.size();
      if ( s > maxSeqLen )
        maxSeqLen = s;
    }
    int seqLen = maxSeqLen;
    #endif

    // Random sequences of 'node'
    vector<string> modelStrIds;
    modelIds( modelStrIds );

    map<string, int> modelIdx;
    int n = modelStrIds.size();
    for ( int i = 0; i < n; i++ )
      modelIdx[ modelStrIds[i] ] = i;

    int ref_model_idx = modelIdx[node->label];
    #if DEBUG_SAMPLE_PP
    fprintf(stderr,"sample_pp(): ref ref_model_idx=%d\n", ref_model_idx);
    #endif

    char **seqTbl;
    MALLOC(seqTbl, char**, ssize * sizeof(char*));
    for ( int i = 0; i < ssize; ++i )
    {
      MALLOC(seqTbl[i], char*, (seqLen+1) * sizeof(char));
      seqTbl[i][seqLen] = '\0';
    }

    sampleMF(seqTbl, ref_model_idx, ssize, seqLen);

    // lpp's of the random sequences of 'node'
    double lpp;
    fprintf(lppFH,"%s", node->label.c_str());
    for ( int i = 0; i < ssize; ++i )
    {
      lpp = normLog10prob(seqTbl[i], seqLen, ref_model_idx);
      fprintf(lppFH,",%lf", lpp);
    }
    fprintf(lppFH,"\n");


    // Siblings of node
    NewickNode_t *pnode = node->parent_m;
    vector<NewickNode_t *> siblings;
    n = pnode->children_m.size();
    for (int i = 0; i < n; i++)
      if ( pnode->children_m[i] != node )
        siblings.push_back(pnode->children_m[i]);

    int nSiblings = (int)siblings.size();
    #if DEBUG_SAMPLE_PP
    fprintf(stderr,"sample_pp(): nSiblings=%d\n", nSiblings);
    #endif

    // For each sib get random seq's from the corresponding model and lpp's from
    // the 'node' model
    NewickNode_t *sibnode;
    for (int j = 0; j < nSiblings; j++)
    {
      sibnode = siblings[j];
      int model_idx = modelIdx[sibnode->label];
      sampleMF(seqTbl, model_idx, ssize, seqLen);

      fprintf(lppFH,"%s", sibnode->label.c_str());
      for ( int i = 0; i < ssize; ++i )
      {
        #if DEBUG_SAMPLE_PP
        fprintf(stderr,"\rj=%d i=%d", j, i);
        #endif
        lpp = normLog10prob(seqTbl[i], seqLen, ref_model_idx);
        fprintf(lppFH,",%lf", lpp);
      }
      fprintf(lppFH,"\n");
    }

    #if DEBUG_SAMPLE_PP
    fprintf(stderr,"sample_pp(): After ji loop\n");
    #endif

    // Cleanup
    for ( int s = 0; s < ssize; ++s )
      free(seqTbl[s]);

    #if DEBUG_SAMPLE_PP
    fprintf(stderr,"Freeing seqTbl\n");
    #endif
    free(seqTbl);

    fclose(lppFH);
}


//------------------------------------------------------ sample_cp ----
/// Generates log10 conditional probabilities of random sequences from the given
/// model's children and saves them to a file in the rLcps subdirectory of the mcDir.
///
/// Random sequences are generated using corresponding models. Conditional
/// probabilities are from the same model.
///
/// node   - A node of the model tree.
/// ssize  - The number of random sequences/llp's to be generated.
/// seqLen - The sequence length of the random sequences
/// mcDir  - A directory holding MC models, fasta files etc, within which rLpps output directory will be created.
/// outDir - An output directory.
void MarkovChains_t::sample_cp( NewickNode_t *node, int ssize, int seqLen, char *mcDir, string &outDir)
{
    #define DEBUG_SAMPLE_CP 0

    if ( node->idx < 0 ) // The node has to be an internal node for it to have children (idx < 0 means internal node)
    {
      // Bulding model index table
      vector<string> modelStrIds;
      modelIds( modelStrIds );

      map<string, int> modelIdx;
      int n = modelStrIds.size();
      for ( int i = 0; i < n; i++ )
        modelIdx[ modelStrIds[i] ] = i;

      // Setting up the output
      char **seqTbl;
      MALLOC(seqTbl, char**, ssize * sizeof(char*));
      for ( int i = 0; i < ssize; ++i )
      {
        MALLOC(seqTbl[i], char*, (seqLen+1) * sizeof(char));
        seqTbl[i][seqLen] = '\0';
      }

      // Setting up output directory
      #if 0
      string outDir = string(mcDir) + string("/rLcps/");
      if ( !dir_exists( outDir.c_str() ) ) { // checking for existence of outDir
        string cmd("mkdir -p ");
        cmd += outDir;
        system(cmd.c_str());
      }
      #endif

      int numChildren = node->children_m.size();
      for (int i = 0; i < numChildren; i++)
      {
        NewickNode_t *child_node = node->children_m[i];

        string lppFile = outDir + node->label + string("__") + child_node->label + string(".csv");
        FILE *lppFH = fOpen( lppFile.c_str(), "w");

        int child_model_idx = modelIdx[child_node->label];
#if DEBUG_SAMPLE_CP
        fprintf(stderr,"sample_cp(): child_model_idx=%d\n", child_model_idx);
#endif
        sampleMF(seqTbl, child_model_idx, ssize, seqLen);

        // lpp's of the random sequences of 'child_node'
        double lpp;
        fprintf(lppFH,"%s", child_node->label.c_str());
        for ( int i = 0; i < ssize; ++i )
        {
          lpp = normLog10prob(seqTbl[i], seqLen, child_model_idx);
          fprintf(lppFH,",%lf", lpp);
        }
        fprintf(lppFH,"\n");

        fclose(lppFH);
      } // END OF for (int i = 0; i < numChildren; i++)

      #if DEBUG_SAMPLE_CP
      fprintf(stderr,"Freeing seqTbl\n");
      #endif
      for ( int i = 0; i < ssize; ++i )
        free(seqTbl[i]);
      free(seqTbl);
    } // END OF if ( node->idx < 0 )
}
