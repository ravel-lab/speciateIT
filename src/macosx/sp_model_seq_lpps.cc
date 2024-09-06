/*
  Generating random samples of sequences from species MC models, estimating
  log10 posterior probabilities on these sequences and the sibling sequences of
  the given species and saving it in a csv file named by the given species name.

  Copyright (C) 2017 Pawel Gajer pgajer@gmail.com and Jacques Ravel jravel@som.umaryland.edu

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

#include <algorithm>
#include <cctype>
#include <getopt.h>
#include <iostream>
#include <limits.h> /* PATH_MAX */
#include <queue>
#include <string>
#include <vector>
#include <filesystem>

#include "CUtilities.h"
#include "IOCUtilities.h"
#include "IOCppUtilities.hh"
#include "CppUtilities.hh"
#include "CppStatUtilities.hh"
#include "MarkovChains.hh"
#include "StatUtilities.hh"
#include "Newick.hh"
#include "CStatUtilities.h"

using namespace std;

//----------------------------------------------------------- printUsage ----
void printUsage( const char *s )
{
  cout << endl

       << "USAGE " << endl
       << endl
       << "  Estimating log10 posterior probabilities of random sequences from species models. The random sequences are from the model and the model's siblings.\n"
       << endl
       << "  " << s << " -d <MC models directory> -s <size>\n"
       << endl
       << "\tOptions:\n"
       << "\t-d | --dir <dir>          - directory containing MC model files and reference sequences fasta files.\n"
       << "\t-s | --sample-size <size> - the number of random sequences per species. The default value: 1000.\n"
       << "\t-v | --verbose            - verbose mode.\n"

       << "\n\tExample: \n"

       << "\t" << s << " -v -d mcDir -s 1000" << endl << endl;
}


//----------------------------------------------- printHelp ----
void printHelp( const char *s )
{
    printUsage(s);

    cout << endl
	 << "  The code generates lpp CSV files in rlpps sub-directory of the mcDir\n"
	 << "  File format: <lpp of ref sp>,<lpp of sib1>,<lpp of sib2>, ... ,<lpp of sibn>\n"
	 << endl << endl;
}

//================================================= inPar_t ====
//! holds input parameters
class inPar_t
{
public:
  inPar_t();
  ~inPar_t();

  char *mcDir;              /// input directory for MC model files
  char *treeFile;           /// reference tree file
  vector<char *> trgFiles;  /// list of paths to fasta training files
  vector<int> kMerLens;     /// list of word lengths
  int maxNumAmbCodes;       /// maximal acceptable number of ambiguity codes for a sequence; above this number log10probIUPAC() returns 1;
  int randSampleSize;       /// number of random sequences of each model (seq length = mean ref seq). If 0, no random samples will be generated.
  int pseudoCountType;      /// pseudo-count type; see MarkovChains.hh for possible values
  bool verbose;

  void print();
};

//------------------------------------------------- constructor ----
inPar_t::inPar_t()
{
  mcDir           = NULL;
  treeFile        = NULL;
  maxNumAmbCodes  = 5;
  randSampleSize  = 1000;
  pseudoCountType = recPdoCount;
  verbose         = false;
}

//------------------------------------------------- constructor ----
inPar_t::~inPar_t()
{
  if ( mcDir )
    free(mcDir);

  if ( treeFile )
    free(treeFile);

  int n = trgFiles.size();
  for ( int i = 0; i < n; ++i )
    free(trgFiles[i]);
}

//------------------------------------------------------- print ----
void inPar_t::print()
{
  cerr << "\npseudoCountType:\t\t" << pseudoCountType
       << "\nverbose:\t\t" << verbose << endl;

  cerr << "mcDir:\t\t";
  if ( mcDir )
    cerr << mcDir << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "treeFile:\t\t";
  if ( treeFile )
    cerr << treeFile << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "trgFiles:\t\t";
  int n = trgFiles.size();
  if ( n )
  {
    for ( int i = 0; i < n; ++i )
      cerr << trgFiles[i] << "\t";
    cerr << endl;
  } else {
    cerr << "MISSING" << endl;
  }

  cerr << "kMerLens:\t\t";
  n = kMerLens.size();
  if ( n )
  {
    for ( int i = 0; i < n; ++i )
      cerr << "\t" << kMerLens[i];
    cerr << endl;
  } else {
    cerr << "MISSING" << endl;
  }
}

//============================== local sub-routines =========================
void parseArgs( int argc, char ** argv, inPar_t *p );

//============================== main ======================================
int main(int argc, char **argv)
{
    #define DEBUG_SP 0

    //-- setting up init parameters
    inPar_t *inPar = new inPar_t();

    //-- parsing input parameters
    parseArgs(argc, argv, inPar);

    if ( inPar->verbose )
      inPar->print();

    // Loading model tree
    NewickTree_t nt;
    string trFile = string(inPar->mcDir) + string("/model.tree");
    STRDUP(inPar->treeFile, trFile.c_str());
    if ( !nt.loadTree(inPar->treeFile) )
    {
      fprintf(stderr,"\n\n\tERROR: Could not load Newick tree from %s\n\n", inPar->treeFile);
      printHelp(argv[0]);
      exit(1);
    }

    if ( inPar->verbose )
      fprintf(stderr,"--- Model tree depth: %d\n", nt.getDepth());

    if ( inPar->verbose )
      cerr << "\r--- Reading k-mer frequency tables from " << inPar->mcDir << " ... ";

    // Loading MC models
    int wordLen = 8;
    MarkovChains_t *probModel = new MarkovChains_t(wordLen-1,
                                                   inPar->mcDir,
                                                   inPar->maxNumAmbCodes,
                                                   inPar->pseudoCountType);
    if (inPar->verbose )
      cerr << "done" << endl;

    vector<char *> modelIds = probModel->modelIds();
    vector<string> modelStrIds;
    probModel->modelIds( modelStrIds );
    nt.modelIdx( modelStrIds );

    //
    // Setting up traversal of the model tree using breath first search
    //
    queue<NewickNode_t *> bfs; // breath first search queue
    NewickNode_t *root = nt.root();
    bfs.push(root);

    NewickNode_t *node;
    int numChildren;
    int nodeCount = 1;

    //
    // Traversing the models tree
    //
    while ( !bfs.empty() )
    {
      node = bfs.front();
      bfs.pop();

      if ( inPar->verbose && node->idx < 0 )
        fprintf(stderr, "--- [%d] Visiting %s idx=%d depth=%d\n", nodeCount, node->label.c_str(), node->idx, node->depth_m);

      if ( node != root )
      {
        nodeCount++;
        //numChildren = node->children_m.size();

        #if 0
        fprintf(stderr,"sp_model_seq_lpps(): %s  numChildren=%d  idx=%d\n",
                node->label.c_str(), numChildren, node->idx);
        #endif

        //if ( numChildren==0 ) // the node is a leaf, that is a species
        if ( node->idx >= 0 )
        {
          if ( inPar->verbose )
            fprintf(stderr, "--- [%d] Processing %s idx=%d depth=%d\n", nodeCount, node->label.c_str(), node->idx, node->depth_m);

          probModel->sample_pp( node, inPar->randSampleSize, inPar->mcDir);
          //sleep(1);

        } // END OF if ( numChildren==0 )

      } // END OF if ( node != root )

      //if ( numChildren )
      if ( node->idx < 0 )
      {
        numChildren = node->children_m.size();
        for (int i = 0; i < numChildren; i++)
          bfs.push(node->children_m[i]);
      }

    } // END OF while ( !bfs.empty() )

    delete probModel;

    return EXIT_SUCCESS;
}



//----------------------------------------------------------- parseArgs ----
//! parse command line arguments
void parseArgs( int argc, char ** argv, inPar_t *p )
{
    int c, errflg = 0;
    optarg = NULL;

    static struct option longOptions[] = {
      {"sample-size"        ,required_argument, 0,          's'},
      {"max-num-amb-codes"  ,required_argument, 0,          'b'},
      {"pseudo-count-type"  ,required_argument, 0,          'p'},
      {"help"               ,no_argument,       0,            0},
      {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv,"b:s:p:d:vh",longOptions, NULL)) != -1)
      switch (c) {
        case 'b':
          p->maxNumAmbCodes = atoi(optarg);
          break;

        case 's':
          p->randSampleSize = atoi(optarg);
          break;

        case 'p':
          {
            int pc = atoi(optarg);
            if ( pc == -1 )
            {
              p->pseudoCountType = zeroOffset0;
            }
            else if ( pc == 0 )
            {
              p->pseudoCountType = zeroOffset1;
            }
            else if ( pc == 1 )
            {
              p->pseudoCountType = zeroOffset1;
            }
            else if ( pc == 2 )
            {
              p->pseudoCountType = recPdoCount;
            }
            else
            {
              cerr << "ERROR in " << __FILE__ << " at line " << __LINE__ << ": Undefined pseudo-count type" << endl;
              exit(1);
            }
          }
          break;

        case 'd':
          p->mcDir = strdup(optarg);
          break;

        case 'v':
          p->verbose = true;
          break;

        case 'h':
          printHelp(argv[0]);
          exit (EXIT_SUCCESS);
          break;

        case 0:
          break;

        default:
          cerr << "\n"
               << "=========================================\n"
               << " ERROR: Unrecognized option " << (char)c << "\n" << endl;

          for ( int i=0; i < argc; ++i )
            cerr << argv[i] << " ";
          cerr << endl;

          cerr << "==========================================\n" << endl;
          ++errflg;
          break;
      }

    if ( errflg )
    {
      printUsage(argv[0]);
      cerr << "Try '" << argv[0] << " -h' for more information" << endl;
      exit (EXIT_FAILURE);
    }
}
