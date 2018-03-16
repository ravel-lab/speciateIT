/*
  Generating for each non-leaf node of the model tree and all ref sequnces of the
  children nodes a table with each seq's pp w/r each child node.

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
       << "  Generating for each non-leaf node of the model tree and all ref sequnces of the\n"
       << "  children nodes a table with each seq's pp w/r each child node.\n"
       << endl
       << "  " << s << " -d < MC models directory>[Options]\n"
       << endl
       << "\tOptions:\n"
       << "\t-d <dir>       - directory containing MC model files and reference sequences fasta files\n"

       << endl
       << "  The following files are generated in the MC directory:\n"
       << "  - nodeName_lpps.txt\n"
       << "    File format: tab delimited sequnce: seqID cltr lpp's w/r all child nodes\n"
       << endl

       << "\n\tExample: \n"

       << "\t" << s << " -v -d Firmicutes_group_5_V3V4_MC_models_dir -o Firmicutes_group_5_V3V4_pp_embedding_dir" << endl << endl;
}


//----------------------------------------------- printHelp ----
void printHelp( const char *s )
{
    printUsage(s);
}

//================================================= inPar_t ====
//! holds input parameters
class inPar_t
{
public:
  inPar_t();
  ~inPar_t();

  char *mcDir;              /// input directory for MC model files
  char *trgFile;            /// file containing paths to fasta training files
  char *treeFile;           /// reference tree file
  double thld;              /// threshold for | log( p(x | M_L) / p(x | M_R) | of the competing models
  vector<char *> trgFiles;  /// list of paths to fasta training files
  vector<int> kMerLens;     /// list of word lengths
  int printCounts;          /// flag initiating print out of word counts
  int maxNumAmbCodes;       /// maximal acceptable number of ambiguity codes for a sequence; above this number log10probIUPAC() returns 1;
  int pseudoCountType;      /// pseudo-count type; see MarkovChains2.hh for possible values
  bool verbose;
  int debug;

  void print();
};

//------------------------------------------------- constructor ----
inPar_t::inPar_t()
{
  mcDir           = NULL;
  trgFile         = NULL;
  treeFile        = NULL;
  thld            = 0.0;
  printCounts     = 0;
  maxNumAmbCodes  = 5;
  pseudoCountType = recPdoCount;
  verbose         = false;
  debug           = 0;
}

//------------------------------------------------- constructor ----
inPar_t::~inPar_t()
{
  if ( mcDir )
    free(mcDir);

  if ( trgFile )
    free(trgFile);

  if ( treeFile )
    free(treeFile);

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
    cerr << endl;

  cerr << "mcDir=\t\t";
  if ( mcDir )
    cerr << mcDir << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "treeFile=\t\t";
  if ( treeFile )
    cerr << treeFile << endl;
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

  if ( inPar->debug )
    inPar->print();

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
    // lets see if we can find ref tree in mcDir/model.tree
    string trFile = string(inPar->mcDir) + "/model.tree";
    STRDUP(inPar->treeFile, trFile.c_str());

    if ( !nt.loadTree(inPar->treeFile) )
    {
      fprintf(stderr,"\n\n\tERROR: Could not load Newick tree from %s\n\n", inPar->treeFile);
      printHelp(argv[0]);
      exit(1);
    }
  }

  int depth = nt.getDepth();
  if ( inPar->debug )
    cerr << "--- Depth of the tree: " << depth << endl;

  int nModels = 0;

  if ( inPar->mcDir ) // extracting number of models and k-mer size
  {
    string modelFile(inPar->mcDir);
    modelFile += "/modelIds.txt";
    FILE *in = fopen(modelFile.c_str(), "r");
    // this should be moved to readLines()
    if ( !in )
    {
      fprintf(stderr, "\n\n\tERROR in %s at line %d: Cannot read model ids.\n\n", __FILE__, __LINE__);
      exit(1);
    }
    fclose(in);

    vector<char *> modelIds;
    readLines(modelFile.c_str(), modelIds);
    nModels = modelIds.size();

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
      inPar->kMerLens.push_back(k);
    }
  }
  else if ( !inPar->mcDir )
  {
    fprintf(stderr, "\n\n\tERROR: in %s at line %d: Please specify a directory with MC model files using -d flag.\n\n", __FILE__, __LINE__);
    printHelp(argv[0]);
    exit(1);
  }


  if ( inPar->kMerLens.size() == 0 )
  {
    int kMers[] = {8};
    fprintf(stderr, "\n\tWARNING: Setting k-mer size to %d.\n", kMers[0]);
    int n = sizeof(kMers) / sizeof(int);
    for ( int i = 0; i < n; ++i )
      inPar->kMerLens.push_back(kMers[i]);
  }

  if ( inPar->verbose )
    cerr << "nModels: " << nModels << endl;

  int wordLen = inPar->kMerLens[0];

  if ( inPar->debug )
  {
    cerr << "\rk=" << wordLen << "\n";

    if ( inPar->mcDir && !inPar->trgFiles.size() )
      cerr << "\r--- Reading k-mer frequency tables from " << inPar->mcDir << " ... ";
    else
      cerr << "\r--- Generating k-mer frequency tables for k=1:" << wordLen << " ... ";
  }

  size_t alloc = 1024*1024;
  char *data, *seq;
  MALLOC(data, char*, alloc * sizeof(char));
  MALLOC(seq, char*, alloc * sizeof(char));

  MarkovChains2_t *probModel;

  // loading MC models
  probModel = new MarkovChains2_t( wordLen-1,
				   inPar->trgFiles,
				   inPar->mcDir,
				   inPar->maxNumAmbCodes,
				   inPar->pseudoCountType );
  if (inPar->debug )
    cerr << "done" << endl;

  vector<char *> modelIds = probModel->modelIds();
  vector<string> modelStrIds;
  probModel->modelIds( modelStrIds );
  nt.modelIdx( modelStrIds );

  // traverse the reference tree using breath first search
  NewickNode_t *node;
  int nChildren;
  int nodeCount = 1;
  map<string, string>::iterator itr;
  map<string, string> seqRecs; // fasta file sequence records
  double lpp;

  queue<NewickNode_t *> bfs;
  NewickNode_t *root = nt.root();
  bfs.push(root);

  if ( inPar->verbose )
    fprintf(stderr, "\n\n");

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    nChildren = node->children_m.size();

    if ( node != root )
    {
      if ( inPar->verbose )
      {
	fprintf(stderr, "--- [%d] Processing %s\n", nodeCount, node->label.c_str());
	nodeCount++;
      }

      if ( nChildren )
      {
	string sFile = string(inPar->mcDir) + string("/") + node->label + string("_ref_lpps.txt");
	FILE *sOut = fOpen( sFile.c_str(), "w");

	// header
	fprintf(sOut,"seqIDs\tcltr");
	for (int i = 0; i < nChildren; i++)
	  fprintf(sOut,"\t%s",node->children_m[i]->label.c_str());
	fprintf(sOut,"\n");
	// end of header

	// pp loop
	for (int i = 0; i < nChildren; i++)
	{
	  string faFile = string(inPar->mcDir) + string("/") + node->children_m[i]->label + string(".fa");
	  seqRecs.clear();
	  readFasta( faFile.c_str(), seqRecs);

	  for ( itr = seqRecs.begin(); itr != seqRecs.end(); ++itr )
	  {
	    fprintf(sOut,"%s\t%d", itr->first.c_str(),i);
	    for (int j = 0; j < nChildren; j++)
	    {
	      lpp = probModel->normLog10prob(itr->second.c_str(), (int)itr->second.size(), node->children_m[j]->model_idx );
	      fprintf(sOut,"\t%f", lpp);
	    }
	    fprintf(sOut,"\n");
	  }
	} // if ( nChildren )
	fclose(sOut);
      }
    } // if ( node != root )

    if ( nChildren )
    {
      for (int i = 0; i < nChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  } // end of while ( !bfs.empty() ) loop

  if ( inPar->verbose )
    fprintf(stderr,"\r\n\n\tOutput written to %s\n\n", inPar->mcDir);

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
    {"print-counts"       ,no_argument, &p->printCounts,    1},
    {"max-num-amb-codes"  ,required_argument, 0,          'b'},
    {"model-tree"         ,required_argument, 0,          'r'},
    {"pseudo-count-type"  ,required_argument, 0,          'p'},
    {"help"               ,no_argument,       0,            0},
    {"debug"              ,no_argument, &p->debug,          0},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv,"b:d:hk:p:r:t:v",longOptions, NULL)) != -1)
    switch (c)
    {
      case 'b':
	p->maxNumAmbCodes = atoi(optarg);
	break;

      case 'd':
	p->mcDir = strdup(optarg);
	break;

      case 'h':
	printHelp(argv[0]);
	exit (EXIT_SUCCESS);
	break;

      case 'k':
	parseCommaList(optarg, p->kMerLens);
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

      case 'r':
	p->treeFile = strdup(optarg);
	break;

      case 't':
	p->trgFile = strdup(optarg);
	break;

      case 'v':
	p->verbose = true;
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
