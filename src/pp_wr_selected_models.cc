/*
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
#include <string.h>
#include <string>
#include <vector>
#include <queue>
#include <sys/time.h>

#include "CUtilities.h"
#include "IOCUtilities.h"
#include "IOCppUtilities.hh"
#include "CppUtilities.hh"
#include "MarkovChains2.hh"
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
       << " Computing posterior probabilities with respect of selected models." << endl
       << endl
       << s << " -d < MC models directory> -m <file with one column of selected model names> -i <input fasta file> -o <output file> [Options]" << endl
       << endl
       << "\tOptions:\n"
       << "\t-d <dir>      - directory containing MC model files\n"
       << "\t-o <dir>      - output file\n"
       << "\t-i <inFile>   - input fasta file with sequences for which -log10(prob(seq | model_i)) are to be computed\n"
       << "\t-s <selModel> - name of a model for which posterior probabilities are to be computed\n"
       << "\t-l <selFile>  - file name containing a list models names for which posterior probabilities are to be computed\n"
       << "\t-q|--quiet    - suppers pregress messages\n"
       << "\t-v|--verbose  - verbose mode\n"
       << "\t-h|--help     - this message\n\n"

       << "\n\tExample: \n"

       << "\ncd Firmicutes_group_6_V3V4_dir\n\n"
       << s << " -d Firmicutes_group_6_V3V4_MC_models_dir -i Firmicutes_group_6_V3V4_MC_models_dir/Lactobacillus_plantarum.fa -s Lactobacillus_plantarum -o Firmicutes_group_6_V3V4_MC_models_dir/Lactobacillus_plantarum.postProbs" << endl << endl;
}


//----------------------------------------------- printHelp ----
void printHelp( const char *s )
{
    cout << "\nComputing posterior probabilities with respect of selected models\n";

    printUsage(s);
}

//----------------------------------------------- errTbl_t ----
/// holds errTbl and log10cProb.thld
typedef struct
{
  double **errTbl; ///
  int nrow;
  double thld;  /// the first element of x
  double *x;    /// the first column of errTbl
  double xmax;  /// the max of x
} errTbl_t;

//================================================= inPar2_t ====
//! holds input parameters
class inPar_t
{
public:
  inPar_t();
  ~inPar_t();

  char *inFile;             /// input fasta file with sequences for which -log10(prob(seq | model_i)) are to be computed
  char *outFile;            /// output file with posterior probabilities of the selected models
  char *selModel;           /// name of selected model
  char *selFile;            /// file name of the file containg names of selected models
  vector<string> selModels; /// list of selected models
  char *mcDir;              /// input directory for MC model files
  char *faDir;              /// directory of reference fasta files
  char *trgFile;            /// file containing paths to fasta training files
  char *treeFile;           /// reference tree file
  double thld;              /// threshold for | log( p(x | M_L) / p(x | M_R) | of the competing models
  vector<char *> trgFiles;  /// list of paths to fasta training files
  vector<int> kMerLens;     /// list of word lengths
  int printCounts;          /// flag initiating print out of word counts
  int skipErrThld;          /// ignore classification error condition - with this option on each sequence is classified to the species with the highest p(x|M)
  int maxNumAmbCodes;       /// maximal acceptable number of ambiguity codes for a sequence; above this number log10probIUPAC() returns 1;
  int randSampleSize;       /// number of random sequences of each model (seq length = mean ref seq). If 0, no random samples will be generated.
  int pseudoCountType;      /// pseudo-count type; see MarkovChains2.hh for possible values
  bool verbose;
  bool quiet;
  bool printNCprobs;        /// if true, the program prints to files normalized conditional probabilities for tuning threshold values of taxon assignment
  int dimProbs;             /// max dimension of probs
  bool revComp;             /// reverse-complement query sequences before processing
  int debug;

  void print();
};

//------------------------------------------------- constructor ----
inPar_t::inPar_t()
{
  outFile         = NULL;
  mcDir           = NULL;
  trgFile         = NULL;
  inFile          = NULL;
  treeFile        = NULL;
  thld            = 0.0;
  printCounts     = 0;
  skipErrThld     = 0;
  maxNumAmbCodes  = 5;
  randSampleSize  = 0;
  pseudoCountType = recPdoCount;
  dimProbs        = 0;
  printNCprobs    = false;
  verbose         = false;
  quiet           = false;
  revComp         = false;
  debug           = 0;
}

//------------------------------------------------- constructor ----
inPar_t::~inPar_t()
{
  if ( outFile )
    free(outFile);

  if ( mcDir )
    free(mcDir);

  if ( trgFile )
    free(trgFile);

  if ( inFile )
    free(inFile);

  if ( treeFile )
    free(treeFile);

  int n = trgFiles.size();
  for ( int i = 0; i < n; ++i )
    free(trgFiles[i]);
}

//------------------------------------------------------- print ----
void inPar_t::print()
{
  fprintf(stderr, "printCounts     %3d\n", printCounts);
  fprintf(stderr, "pseudoCountType %3d\n", pseudoCountType);

  if ( trgFile )
    fprintf(stderr, "trgFile         %s\n", trgFile);
  else
    fprintf(stderr, "trgFile         MISSING\n");

  cerr << "mcDir=\t\t";
  if ( mcDir )
    cerr << mcDir << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "outFile=\t\t";
  if ( outFile )
    cerr << outFile << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "inFile=\t\t";
  if ( inFile )
    cerr << inFile << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "treeFile=\t\t";
  if ( treeFile )
    cerr << treeFile << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "selFile=\t\t";
  if ( selFile )
    cerr << selFile << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "selModel=\t\t";
  if ( selModel )
    cerr << selModel << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "trgFiles:\t";
  int n = trgFiles.size();
  for ( int i = 0; i < n; ++i )
    cerr << trgFiles[i] << "\t";
  cerr << endl;

  cerr << "kMerLens:";
  n = kMerLens.size();
  for ( int i = 0; i < n; ++i )
    cerr << "\t" << kMerLens[i];
  cerr << endl;

  cerr << "skipErrThld: " << skipErrThld << endl;
}

//============================== local sub-routines =========================
void parseArgs( int argc, char ** argv, inPar_t *p );
bool dComp (double i, double j) { return (i>j); }

//============================== main ======================================
int main(int argc, char **argv)
{
  #define DEBUGMAIN 0
  #define DEBUGMAIN1 0

  struct timeval  tvStart, tvCurrent;
  gettimeofday(&tvStart, NULL);

  //-- setting up init parameters
  inPar_t *inPar = new inPar_t();

  //-- parsing input parameters
  parseArgs(argc, argv, inPar);

  if ( inPar->verbose )
    inPar->print();

  if ( !inPar->inFile )
  {
    fprintf(stderr, "\n\ntERROR in %s at line %d: Input fasta file is missing. Please specify it with the -i flag.\n\n", __FILE__, __LINE__);
    printHelp(argv[0]);
    exit(1);
  }

  if ( !inPar->outFile )
  {
    fprintf(stderr, "\n\n\tERROR in %s at line %d: Missing output file. Please specify it with the -o flag.\n\n", __FILE__, __LINE__);
    printHelp(argv[0]);
    exit(1);
  }

  if ( inPar->selModel )
  {
    inPar->selModels.push_back( string(inPar->selModel) );
  }
  else if ( inPar->selFile )
  {
    readLines(inPar->selFile, inPar->selModels); // copies selected models from inPar->selFile into inPar->selModels
    free( inPar->selFile );
  }
  else
  {
    fprintf(stderr, "\n\ntERROR in %s at line %d: Missing both selected model and selected models file (one of them needs to be specified)\n", __FILE__, __LINE__);
    fprintf(stderr, "\tPlease specify the first with the -s flag and the second with the -l flag.\n\n");
    printHelp(argv[0]);
    exit(1);
  }

  if ( inPar->trgFile )
  {
    readLines(inPar->trgFile, inPar->trgFiles); // path(s) from inPar->trgFile are loaded into inPar->trgFiles
    free( inPar->trgFile );
  }

  int nModels = 0;
  if ( inPar->trgFiles.size() )
  {
    nModels = inPar->trgFiles.size();
  }

  map<string, errTbl_t *> modelErrTbl;
  map<string, double> thldTbl;
  if ( inPar->mcDir ) // extracting number of models and k-mer size
  {
    string inFile(inPar->mcDir);
    inFile += "/modelIds.txt";
    FILE *in = fopen(inFile.c_str(), "r");

    if ( in )
    {
      vector<char *> modelIds;
      readLines(inFile.c_str(), modelIds);
      nModels = modelIds.size();

      if ( !inPar->skipErrThld )
      {
	// reading ncProbThlds.txt file
	string file = string(inPar->mcDir) + string("/ncProbThlds.txt");
	double **thlds;
	int nrow, ncol;
	char **rowNames;
	char **colNames;
	readTable( file.c_str(), &thlds, &nrow, &ncol, &rowNames, &colNames );
	for ( int i = 0; i < nrow; i++ )
	  thldTbl[string(rowNames[i])] = thlds[i][0];

#if 0
	map<string, double>::iterator it;
	cerr << "thldTbl" << endl;
	for ( it = thldTbl.begin(); it != thldTbl.end(); it++ )
	  cerr << it->first << "\t" << it->second << endl;
	cerr << endl;
	exit(1);
#endif

	// reading _error.txt files
	for ( int i = 0; i < nModels; ++i )
	{
	  string file = string(inPar->mcDir) + string("/") + string(modelIds[i]) + string("_error.txt");
	  double **errTbl;
	  int nrow, ncol;
	  int header = 0;
	  readMatrix( file.c_str(), &errTbl, &nrow, &ncol, header );

	  errTbl_t *errObj = new errTbl_t;
	  errObj->errTbl = errTbl;
	  errObj->nrow = nrow;

	  if ( ncol != 2 )
	  {
	    fprintf(stderr, "\n\n\tERROR: in %s at line %d: errTbl should have two columns and ncol=%d", __FILE__, __LINE__, ncol);
	    fprintf(stderr, "%s: \n", modelIds[i]);
	    printDblTbl(errTbl, nrow, ncol);
	    exit(1);
	  }

	  errObj->thld = errTbl[0][0];

	  MALLOC(errObj->x, double*, nrow * sizeof(double));
	  for ( int j = 0; j < nrow; ++j )
	    errObj->x[j] = errTbl[j][0];

	  errObj->xmax = errTbl[nrow-1][0];
	  modelErrTbl[ modelIds[i] ] = errObj;
	}
      }

      for ( int i = 0; i < nModels; ++i )
	free(modelIds[i]);

      int k = 0;
      char countStr[5];
      sprintf(countStr,"%d",k);
      string file = string(inPar->mcDir) + string("/MC") + string(countStr) + string(".log10cProb");

      while ( exists( file.c_str() ) )
      {
	k++;
	sprintf(countStr,"%d",k);
	file = string(inPar->mcDir) + string("/MC") + string(countStr) + string(".log10cProb");
      }

      if ( (inPar->kMerLens.size() && inPar->kMerLens[0] > k) )
      {
	inPar->kMerLens[0] = k;
      }
      else if ( !inPar->kMerLens.size() )
      {
	//inPar->kMerLens.clear();
	inPar->kMerLens.push_back(k);
      }
    } // end of if ( in )
    else
    {
      fprintf(stderr, "\n\n\tERROR: in %s at line %d: Cannot open %s\n\n", __FILE__, __LINE__, inFile.c_str());
      exit(1);
    }

    fclose(in);
  }
  else if ( !inPar->mcDir && !inPar->trgFile )
  {
    fprintf(stderr, "\n\n\tERROR: in %s at line %d: Please specify a directory with MC model files using -d flag.\n\n", __FILE__, __LINE__);
    printHelp(argv[0]);
    exit(1);
  }

  NewickTree_t nt;
  if ( inPar->treeFile ) // load ref tree
  {
    if ( !nt.loadTree(inPar->treeFile) )
    {
      fprintf(stderr,"\n\n\tERROR: Could not load Newick tree from %s\n\n", inPar->treeFile);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    // lets see if we can find ref tree in mcDir
    // model.tree
    string trFile = string(inPar->mcDir) + "/model.tree";
    STRDUP(inPar->treeFile, trFile.c_str());

    if ( !nt.loadTree(inPar->treeFile) )
    {
      fprintf(stderr, "\n\n\tERROR: in %s at line %d: reference tree file is missing and I have trouble loading %s. Please specify it with the -r flag.\n\n", __FILE__, __LINE__, inPar->treeFile);
      printHelp(argv[0]);
      exit(1);
    }
  }

  int depth = nt.getDepth();
  if ( inPar->verbose )
    cerr << "--- Depth of the reference tree: " << depth << endl;

  if ( inPar->kMerLens.size() == 0 )
  {
    int kMers[] = {3};
    cerr << endl << "WARNING: Setting k-mer size to " << kMers[0] << endl;
    int n = sizeof(kMers) / sizeof(int);
    for ( int i = 0; i < n; ++i )
      inPar->kMerLens.push_back(kMers[i]);
  }

  if ( inPar->verbose )
    cerr << "--- Number of Models: " << nModels << endl;

  int wordLen = inPar->kMerLens[0];

  if ( inPar->verbose )
  {
    cerr << "\rk=" << wordLen << "\n";

    if ( inPar->mcDir && !inPar->trgFiles.size() )
      cerr << "\r--- Reading conditional probabilities tables from " << inPar->mcDir << " ... ";
    else
      cerr << "\r--- Generating k-mer frequency tables for k=1:" << wordLen << " ... ";
  }

  MarkovChains2_t *probModel;
  probModel = new MarkovChains2_t( wordLen-1,
				   inPar->trgFiles,
				   inPar->mcDir,
				   inPar->maxNumAmbCodes,
				   inPar->pseudoCountType );
  if ( inPar->verbose )
    cerr << "done" << endl;

  vector<char *> modelIds = probModel->modelIds();
  vector<string> modelStrIds;
  probModel->modelIds( modelStrIds );

  #if 0
  map<string, errTbl_t *>::iterator itr = modelErrTbl.begin();
  for ( ; itr != modelErrTbl.end(); ++itr )
  {
    errTbl_t *errObj = itr->second;
    fprintf(stderr, "\n\n%s\tthld=%f\n", itr->first.c_str(), errObj->thld);
    printDblTbl(errObj->errTbl, errObj->nrow, errObj->ncol);
  }
  exit(1);
  #endif

  nt.modelIdx( modelStrIds );

  map<string, int> modelIdx;
  //int n = modeStrIds.size();
  for ( int i = 0; i < nModels; i++ )
    modelIdx[ modelStrIds[i] ] = i;

  map<string, int> selModelIdx;
  int nSelModels = inPar->selModels.size();
  for ( int i = 0; i < nSelModels; i++ )
  {
    map<string, int>::iterator it = modelIdx.find( inPar->selModels[i] );
    if ( it != modelIdx.end() )
    {
      selModelIdx[it->first] = it->second;
    }
    else
    {
      fprintf(stderr, "\n\n\tERROR in %s at line %d: Could not find %s in modelIdx\n\n",
	      __FILE__, __LINE__, inPar->selModels[i].c_str());
      exit(1);
    }
  }


  char str[10];
  sprintf(str,"%d",(wordLen-1));

  // ==== computing probabilities of each sequence of inFile to come from each of the MC models ====
  vector<string> path; // decision path tranced from the root to the final node for each query sequence
  char *id;
  int count = 0;
  //int depthCount;
  string currentLabel;
  double x[nModels]; // stores conditional probabilities p(x | M) for children of each node. the root node has 3 children

  FILE *out = fOpen(inPar->outFile, "w");
  FILE *in = fOpen(inPar->inFile, "r");

  int nRecs = numRecordsInFasta( inPar->inFile );
  int q01 = 0;

  if ( nRecs > 1000 )
    q01 = int(0.01 * nRecs);

  pair<string, double> tx2score;
  // vector< pair<string, double> > score; // quality score table taxonomy => score, where
  //                                       // the score is a measure of uncertainty about
  //                                       // classification of a given sequence to 'taxonomy'

  int seqLen;
  size_t alloc = 1024*1024;
  char *data, *seq, *rcseq;
  MALLOC(data, char*, alloc * sizeof(char));
  MALLOC(seq, char*, alloc * sizeof(char));
  MALLOC(rcseq, char*, alloc * sizeof(char));
  //double x1, x2;

  double *probs;
  MALLOC(probs, double*, alloc * sizeof(double));

  if ( inPar->verbose )
    cerr << "--- Number of sequences in " << inPar->inFile << ": " << nRecs << endl;

  int runTime;
  int timeMin = 0;
  int timeSec = 0;
  int perc;
  map<string,int>::iterator it;

  while ( getNextFastaRecord( in, id, data, alloc, seq, seqLen) )
  {
    if ( inPar->verbose && q01 && (count % q01) == 0 )
    {
      perc = (int)( (100.0*count) / nRecs);

      gettimeofday(&tvCurrent, NULL);
      runTime = tvCurrent.tv_sec  - tvStart.tv_sec;

      if ( runTime > 60 )
      {
	timeMin = runTime / 60;
	timeSec = runTime % 60;
      }
      else
      {
	timeSec = runTime;
      }
      fprintf(stderr,"\r%d:%02d  %d [%02d%%]", timeMin, timeSec, count, perc);
    }
    count++;

    if ( inPar->revComp )
    {
      for ( int j = 0; j < seqLen; ++j )
	rcseq[j] = Complement(seq[seqLen-1-j]);
      rcseq[seqLen] = '\0';
    }

    fprintf(out,"%s", id) ;
    //for ( int j = 0; j < nSelModels; j++ )  {
    int i = 0;
    for (it = selModelIdx.begin(); it != selModelIdx.end(); ++it, ++i)
    {
      if ( inPar->revComp )
      {
	x[i] = probModel->normLog10prob(rcseq, seqLen, it->second );
      }
      else
      {
	x[i] = probModel->normLog10prob(seq, seqLen, it->second );
      }
      fprintf(out,"\t%f", x[i]) ;
    }
    //}
    fprintf(out,"\n");

  } // end of   while ( getNextFastaRecord( in, id, data, alloc, seq, seqLen) )

  fclose(in);
  fclose(out);

#if DEBUGMAIN
  fclose(debugout);
#endif

  //fprintf(stderr,"\nNumber of errors: %d / %d (%.2f%%)\n\n", errorCount, nRecs, 100.0*errorCount / (double)nRecs );


  free(seq);
  free(data);

  // It may be a nice idea to report the number of species found

  gettimeofday(&tvCurrent, NULL);
  runTime = tvCurrent.tv_sec  - tvStart.tv_sec;

  if ( runTime > 60 )
  {
    timeMin = runTime / 60;
    timeSec = runTime % 60;
  }
  else
  {
    timeSec = runTime;
  }

  if ( inPar->verbose )
  {
    fprintf(stderr,"\r                                                                       \n");
    fprintf(stderr,"    Elapsed time: %d:%02d                                              \n", timeMin, timeSec);
    fprintf(stderr,"    Output written to %s\n", inPar->outFile);
  }

  #if DEBUGMAIN
  fprintf(stderr,"\n\nDEBUGING Output written to %s\n\n\n", debugFile.c_str());
  #endif

  #if SPPDEBUG
  fprintf(stderr,"\n\nDEBUGING Output written to spp_pprob.csv\n\n");
  #endif

  return EXIT_SUCCESS;
}



//----------------------------------------------------------- parseArgs ----
//! parse command line arguments
void parseArgs( int argc, char ** argv, inPar_t *p )
{
  int c, errflg = 0;
  optarg = NULL;

  static struct option longOptions[] = {
    {"print-counts"       ,no_argument, &p->printCounts,    1},
    //{"skip-err-thld"      ,no_argument, &p->skipErrThld,    1},
    {"skip-err-thld"      ,no_argument,       0, 'x'},
    {"max-num-amb-codes"  ,required_argument, 0, 'b'},
    {"out-file"            ,required_argument, 0, 'o'},
    {"ref-tree"           ,required_argument, 0, 'r'},
    {"pseudo-count-type"  ,required_argument, 0, 'p'},
    {"sel-model"          ,required_argument, 0, 's'},
    {"sel-models-file"    ,required_argument, 0, 'l'},
    {"rev-comp"           ,no_argument,       0, 'c'},
    {"quiet"              ,no_argument,       0, 'q'},
    {"verbose"            ,no_argument,       0, 'v'},
    {"help"               ,no_argument,       0, 'h'},
    {"debug"              ,no_argument, &p->debug, 0},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv,"a:b:c:d:e:g:t:i:k:o:vpq:r:hs:l:y:x",longOptions, NULL)) != -1)
    switch (c)
    {
      case 'a':
	p->dimProbs = atoi(optarg);
	break;

      case 'b':
	p->maxNumAmbCodes = atoi(optarg);
	break;

      case 'c':
	p->revComp = true;
	break;

      case 's':
	p->selModel = strdup(optarg);
	break;

      case 'l':
	p->selFile = strdup(optarg);
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

      case 'o':
	p->outFile = strdup(optarg);
	break;

      case 't':
	p->trgFile = strdup(optarg);
	break;

      case 'g':
	p->faDir = strdup(optarg);
	break;

      case 'i':
	p->inFile = strdup(optarg);
	break;

      case 'k':
	parseCommaList(optarg, p->kMerLens);
	break;

      case 'v':
	p->verbose = true;
	break;

      case 'q':
	p->quiet = true;
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

  for ( ; optind < argc; ++ optind )
    p->trgFiles.push_back( strdup(argv[optind]) );
}
