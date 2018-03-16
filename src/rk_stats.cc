/*
  Generating rank stats for taxons with at least 50 reference sequences.

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
//#include <regex>

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
       << "  Generating rank stats for taxons with at least 50 reference sequences.\n"
       << endl
       << "  For each taxon and record the highest and the second highest posterior\n"
       << "  probability of the sequence.\n"
       << "  These estimated should give us an idea about the cases of novel species\n"
       << "  and how to avoid classifying them to the known species\n"
       << endl
       << "  " << s << " -d < MC models directory> -o <output directory> [Options]\n"
       << endl
       << "\tOptions:\n"
       << "\t-d <dir>       - directory containing MC model files and reference sequences fasta files\n"
       << "\t-o <dir>       - output directory\n"

       << "\n\tExample: \n"

       << "\t" << s << " -v -d Firmicutes_group_5_V3V4_MC_models_dir -o Firmicutes_group_5_V3V4_rk_stats_dir" << endl << endl;
}


//----------------------------------------------- printHelp ----
void printHelp( const char *s )
{
    printUsage(s);

    cout << endl
	 << "  File generated:  rk.stats\n"
	 << "  File format: <taxonName> followed by\n"
	 << "  tab delimited sequence of log posterior probabilities pairs <highest>\t<second highest>\n"
	 << endl << endl;
}

//================================================= inPar_t ====
//! holds input parameters
class inPar_t
{
public:
  inPar_t();
  ~inPar_t();

  char *outDir;             /// output directory for MC taxonomy files
  char *mcDir;              /// input directory for MC model files
  char *trgFile;            /// file containing paths to fasta training files
  char *faDir;              /// directory of reference fasta files
  char *inFile;             /// input file with path(s) to fasta file(s) containing sequences
                            /// for which -log10(prob(seq | model_i)) are to be computed
  char *seqID;              /// sequence ID of a sequence from the training fasta files that is to be excluded
                            /// from model building and needs to be used for cross validation
  char *treeFile;           /// reference tree file
  double thld;              /// threshold for | log( p(x | M_L) / p(x | M_R) | of the competing models
  vector<char *> trgFiles;  /// list of paths to fasta training files
  vector<int> kMerLens;     /// list of word lengths
  int printCounts;          /// flag initiating print out of word counts
  int maxNumAmbCodes;       /// maximal acceptable number of ambiguity codes for a sequence; above this number log10probIUPAC() returns 1;
  int randSampleSize;       /// number of random sequences of each model (seq length = mean ref seq). If 0, no random samples will be generated.
  int pseudoCountType;      /// pseudo-count type; see MarkovChains2.hh for possible values
  bool verbose;
  int debug;

  void print();
};

//------------------------------------------------- constructor ----
inPar_t::inPar_t()
{
  outDir          = NULL;
  mcDir           = NULL;
  faDir           = NULL;
  trgFile         = NULL;
  inFile          = NULL;
  treeFile        = NULL;
  seqID           = NULL;
  thld            = 0.0;
  printCounts     = 0;
  maxNumAmbCodes  = 5;
  randSampleSize  = 1000;
  pseudoCountType = recPdoCount;
  verbose         = false;
  debug           = 0;
}

//------------------------------------------------- constructor ----
inPar_t::~inPar_t()
{
  if ( outDir )
    free(outDir);

  if ( mcDir )
    free(mcDir);

  if ( faDir )
    free(faDir);

  if ( trgFile )
    free(trgFile);

  if ( inFile )
    free(inFile);

  if ( seqID )
    free(seqID);

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

  cerr << "faDir=\t\t";
  if ( faDir )
    cerr << faDir << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "outDir=\t\t";
  if ( outDir )
    cerr << outDir << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "inFile=\t\t";
  if ( inFile )
    cerr << inFile << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "seqID=\t\t";
  if ( seqID )
    cerr << seqID << endl;
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
bool dComp (double i, double j) { return (i>j); }

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

  if ( inPar->outDir ) // create output directory
  {
    string cmd("rm -rf ");
    cmd += string(inPar->outDir);
    cmd += string("; mkdir ");
    cmd += string(inPar->outDir);
    system(cmd.c_str());
  }
  else
  {
    fprintf(stderr, "\n\n\tERROR in %s at line %d: Missing output directory. Please specify it with the -o flag.\n\n", __FILE__, __LINE__);
    printHelp(argv[0]);
    exit(1);
  }

  int nModels = 0;

  if ( inPar->mcDir ) // extracting number of models and k-mer size
  {
    string inFile(inPar->mcDir);
    inFile += "/modelIds.txt";
    FILE *in = fopen(inFile.c_str(), "r");
    if ( !in )
    {
      fprintf(stderr, "\n\n\tERROR in %s at line %d: Cannot read model ids.\n\n", __FILE__, __LINE__);
      exit(1);
    }
    fclose(in);

    vector<char *> modelIds;
    readLines(inFile.c_str(), modelIds);
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
      //inPar->kMerLens.clear();
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


  string outFile = string(inPar->outDir) + string("/rk.stats");
  FILE *outFH = fOpen( outFile.c_str(), "w");

  // traverse the reference tree using breath first search
  NewickNode_t *node;
  NewickNode_t *sibnode;
  NewickNode_t *pnode;
  int numChildren;
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

    numChildren = node->children_m.size();

    if ( node != root )
    {
      if ( inPar->verbose )
      {
	fprintf(stderr, "--- [%d] Processing %s\n", nodeCount, node->label.c_str());
	nodeCount++;
      }

      // Identify siblings of node
      pnode = node->parent_m;
      vector<NewickNode_t *> siblings;
      int n = pnode->children_m.size();
      for (int i = 0; i < n; i++)
	if ( pnode->children_m[i] != node )
	  siblings.push_back(pnode->children_m[i]);

      int nSiblings = (int)siblings.size();

      #if 0
      //debug
      fprintf(stderr, "\tIdentifying siblings of %s\n", node->label.c_str());
      fprintf(stderr, "\tSiblings:\n");
      for (int i = 0; i < (int)siblings.size(); ++i )
	fprintf(stderr, "\t\t%s\n", siblings[i]->label.c_str());
      fprintf(stderr, "\n");
      #endif

      //
      // Recording the highest and the second highest log posterior probabilities
      // for each reference sequence of 'node'
      //
      string faFile = string(inPar->mcDir) + string("/") + node->label + string(".fa");
      seqRecs.clear();
      readFasta( faFile.c_str(), seqRecs);

      fprintf(outFH,"%s", node->label.c_str());
      for ( itr = seqRecs.begin(); itr != seqRecs.end(); ++itr )
      {
	lpp = probModel->normLog10prob(itr->second.c_str(), (int)itr->second.size(), node->model_idx );
	double refLogPP = lpp;

	int i = 0;
	sibnode = siblings[i];
	lpp = probModel->normLog10prob(itr->second.c_str(), (int)itr->second.size(), sibnode->model_idx );
	double sibMaxLogPP = lpp;
	//int secondHighestIdx = i;

	for (i = 1; i < nSiblings; i++)
	{
	  sibnode = siblings[i];
	  lpp = probModel->normLog10prob(itr->second.c_str(), (int)itr->second.size(), sibnode->model_idx );
	  if ( lpp > sibMaxLogPP )
	  {
	    sibMaxLogPP = lpp;
	    //secondHighestIdx = i;
	  }
	}

	#if 0
	if ( sibMaxLogPP > refLogPP ) // if we were to care about the max pp and the second largest pp
	{
	  refLogPP = sibMaxLogPP;

	  lpp = probModel->normLog10prob(itr->second.c_str(), (int)itr->second.size(), node->model_idx );
	  sibMaxLogPP = lpp;

	  for (i = 0; i < nSiblings; i++)
	  {
	    if ( i != secondHighestIdx )
	    {
	      sibnode = siblings[i];
	      lpp = probModel->normLog10prob(itr->second.c_str(), (int)itr->second.size(), sibnode->model_idx );
	      if ( lpp > sibMaxLogPP )
	      {
		sibMaxLogPP = lpp;
		secondHighestIdx = i;
	      }
	    }
	  }
	}
	#endif

	fprintf(outFH,"\t%.3f\t%.3f", pow(10, refLogPP), pow(10, sibMaxLogPP));
      }
      fprintf(outFH,"\n");
    } // if ( node != root )

    if ( numChildren )
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  } // end of while ( !bfs.empty() ) loop

  fclose(outFH);

  if ( inPar->verbose )
    fprintf(stderr,"\r\n\n\tOutput written to %s\n\n", inPar->outDir);

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
    {"fasta-dir"          ,required_argument, 0,          'f'},
    {"out-dir"            ,required_argument, 0,          'o'},
    {"ref-tree"           ,required_argument, 0,          'r'},
    {"pseudo-count-type"  ,required_argument, 0,          'p'},
    {"sample-size"        ,required_argument, 0,          's'},
    {"help"               ,no_argument,       0,            0},
    {"debug"              ,no_argument, &p->debug,          0},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv,"b:d:e:f:t:i:k:o:vp:r:s:h",longOptions, NULL)) != -1)
    switch (c)
    {
      case 'b':
	p->maxNumAmbCodes = atoi(optarg);
	break;

      case 's':
	p->randSampleSize = atoi(optarg);
	break;

      case 'r':
	p->treeFile = strdup(optarg);
	break;

      case 'e':
	p->seqID = strdup(optarg);
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

      case 'f':
	p->faDir = strdup(optarg);
	break;

      case 'o':
	p->outDir = strdup(optarg);
	break;

      case 't':
	p->trgFile = strdup(optarg);
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
