/*

  Streamlined version of buildMC that builds MC models. If MC models exist, they will be overwritten.

  It is assumed that buildMC is run after buildModelTree was run.


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

#include <getopt.h>
#include <string>
#include <vector>

#include "CUtilities.h"
#include "CStatUtilities.h"
#include "IOCUtilities.h"
#include "IOCppUtilities.hh"
#include "CppUtilities.hh"
#include "MarkovChains.hh"
#include "StatUtilities.hh"

using namespace std;

//----------------------------------------------------------- printUsage ----
void printUsage( const char *s )
{
    cout << endl

         << "USAGE " << endl << endl
         << s << "-d <MC models directory> [Options]" << endl << endl

         << "\twhere"
         << endl
         << "\tOptions:\n"
         << "\t-d <mcDir>   - directory for MC model files.\n"
         << "\t-t <trgFile> - file containing paths to training fasta files. If not specified, buildMC will use mcDir/tx_fasta_paths.txt file.\n"
         << "\t-k <K>       - K is the k-mer size.\n"
         << "\t--pseudo-count-type, -p <f>  - f=0 for add 1 to all k-mer counts zero-offset.\n"
         << "\t                               f=1 for add 1/4^k to k-mer counts zero-offset.\n"
         << "\t                               f=2 the pseudocounts for a order k+1 model be alpha*probabilities from\n"
         << "\t                                   an order k model, recursively down to pseudocounts of alpha/num_letters\n"
         << "\t                                   for an order 0 model.\n"
         << "\t-v - verbose mode.\n\n"
         << "\t-h|--help - this message\n\n"

         << "\tConditional probabitity tables are store in\n"
         << "\t<file_i>.MC<order>.log10cProb\n\n"

         << "\n\tOutput file format:\n\n"

         << "\tseqId   model1        model2 ...\n"
         << "\tseq_1   log10prob11   log10prob12  ...\n"
         << "\tseq_2   log10prob21   log10prob22  ...\n"
         << "\t...\n"
         << "\n\twhere log10prob_ij is log10 of the prob(seq_i | model_j)\n\n"

         << "\n\tExample: \n"
         << s << " -d sIT_models " << endl << endl;
}


//----------------------------------------------- printHelp ----
void printHelp( const char *s )
{
    cout << endl
         << "Builds Markov chain models using fasta files in sIT_models/tx_fasta_paths.txt" << endl;

    printUsage(s);
}


//================================================= inPar2_t ====
//! holds input parameters
class inPar_t
{
public:
  inPar_t();
  ~inPar_t();

  char *mcDir;              /// input directory for MC model files
  char *trgFile;            /// file containing paths to fasta training files
                            /// for which -log10(prob(seq | model_i)) are to be computed
  vector<char *> trgFiles;  /// list of paths to fasta training files
  int kMerLen;              /// max k-mer length in MC models
  int printCounts;          /// flag initiating print out of word counts
  int maxNumAmbCodes;       /// maximal acceptable number of ambiguity codes for a sequence; above this number log10probIUPAC() returns 1;
  int pseudoCountType;      /// pseudo-count type; see MarkovChains.hh for possible values
  bool verbose;

  void print();
};

//------------------------------------------------- constructor ----
inPar_t::inPar_t()
{
  mcDir           = NULL;
  trgFile         = NULL;
  printCounts     = 0;
  maxNumAmbCodes  = 5;
  pseudoCountType = recPdoCount;
  kMerLen         = 8;
  verbose         = false;
}

//------------------------------------------------- constructor ----
inPar_t::~inPar_t()
{
  if ( mcDir )
    free(mcDir);

  if ( trgFile )
    free(trgFile);

  int n = trgFiles.size();
  for ( int i = 0; i < n; ++i )
    free(trgFiles[i]);
}

//------------------------------------------------------- print ----
void inPar_t::print()
{
  cerr << "printCounts=\t" << printCounts
       << "\npseudoCountType=\t" << pseudoCountType
       << "\nverbose=\t" << verbose
       << "\ntrgFile=\t";

  if ( trgFile )
    cerr << trgFile << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "mcDir=\t\t";
  if ( mcDir )
    cerr << mcDir << endl;
  else
    cerr << "MISSING" << endl;

  int n = (int)trgFiles.size();
  if ( n < 10 ) {

    cerr << "trgFiles:\t";
    for ( int i = 0; i < n; ++i )
      cerr << trgFiles[i] << "\t";
    cerr << endl;
  }

  cerr << "kMerLen: " << kMerLen << endl;
}



//============================== local sub-routines =========================
void parseArgs( int argc, char ** argv, inPar_t *p );

//============================== main ======================================
int main(int argc, char **argv)
{
    //-- setting up init parameters
    inPar_t *inPar = new inPar_t();

    //-- parsing input parameters
    parseArgs(argc, argv, inPar);

    //if ( inPar->verbose )
    //  inPar->print();

    if ( !inPar->mcDir )
    {
      fprintf(stderr, "ERROR: in file %s at line %d: Missing MC modeles directory. Please specify it using -d flag.",
              __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }

    if ( inPar->trgFile )
    {
      readLines(inPar->trgFile, inPar->trgFiles); // path(s) from inPar->trgFile are loaded into inPar->trgFiles
    }
    else
    {
      string trgFile(inPar->mcDir);
      trgFile += string("/tx_fasta_paths.txt");
      inPar->trgFile = strdup( trgFile.c_str() );
      CHECK_FILE(inPar->trgFile);

      if ( inPar->verbose )
        fprintf(stderr, "--- Reading %s ... ", inPar->trgFile);

      readLines(inPar->trgFile, inPar->trgFiles); // path(s) from inPar->trgFile are loaded into inPar->trgFiles

      if ( inPar->verbose )
        fprintf(stderr, "DONE\n");
    }
    free( inPar->trgFile );

    if ( inPar->verbose )
      fprintf(stderr, "--- Generating k-mer (k=%d) frequency tables\n", inPar->kMerLen);

    MarkovChains_t probModel(inPar->kMerLen-1,
                             inPar->trgFiles,
                             inPar->mcDir,
                             inPar->maxNumAmbCodes,
                             inPar->pseudoCountType,
                             inPar->verbose);
    if ( inPar->verbose )
    {
      //fprintf(stderr, " DONE\n");
      fprintf(stderr, "\nMarkov chain models written to %s\n", inPar->mcDir);
    }

    return EXIT_SUCCESS;
}



//----------------------------------------------------------- parseArgs ----
//! parse command line arguments
void parseArgs( int argc, char ** argv, inPar_t *p )
{
    int c, errflg = 0;
    optarg = NULL;

    static struct option longOptions[] = {
      {"max-num-amb-codes"  ,required_argument, 0,          'b'},
      {"pseudo-count-type"  ,required_argument, 0,          'p'},
      {"help"               ,no_argument, 0,                  0},
      {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv,"b:d:t:k:vp:h",longOptions, NULL)) != -1)
      switch (c)
      {
        case 'b':
          p->maxNumAmbCodes = atoi(optarg);
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
              exit(EXIT_FAILURE);
            }
          }
          break;

        case 'd':
          p->mcDir = strdup(optarg);
          break;

        case 't':
          p->trgFile = strdup(optarg);
          break;

        case 'k':
          p->kMerLen = atoi(optarg);
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
