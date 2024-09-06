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
       << " Using prebuilt MC models to classify sequences of an input fasta file" << endl
       << endl
       << s << " -d < MC models directory> -r <ref tree> -i <input fasta file> -o <output directory> [Options]" << endl
       << endl
       << "\tOptions:\n"
       << "\t-d <dir>      - directory containing MC model files\n"
       << "\t-o <dir>      - output directory for MC taxonomy files\n"
       << "\t-i <inFile>   - input fasta file with sequences for which -log10(prob(seq | model_i)) are to be computed\n"
       << "\t-r <ref tree> - reference tree with node labels corresponding to the names of the model files\n"
       << "\t-t <trgFile>  - file containing paths to training fasta files\n"
       << "\t-f <fullTx>   - fullTx file. Its optional parameter for printing classification output in a long format like in RDP classifier\n"
       << "\t-g <faDir>    - directory with reference fasta files\n"
       << "\t-e <seqID>    - sequence ID of a sequence from the training fasta files that is to be excluded\n"
       << "                  from model building and needs to be used for cross validation\n"

       << "\t--rev-comp, -c          - reverse complement query sequences before computing classification posterior probabilities\n"
       << "\t--max-num-amb-codes <n> - maximal acceptable number of ambiguity codes for a sequence\n"
       << "\t                          above this number sequence's log10prob() is not computed and\n"
       << "\t                          the sequence's id it appended to <genus>_more_than_<n>_amb_codes_reads.txt file.\n"
       << "\t                          Default value: 5\n\n"
       << "\t--pseudo-count-type, -p <f>  - f=0 for add 1 to all k-mer counts zero-offset\n"
       << "\t                               f=1 for add 1/4^k to k-mer counts zero-offset\n"
       << "\t                               f=2 the pseudocounts for a order k+1 model be alpha*probabilities from\n"
       << "\t                                   an order k model, recursively down to pseudocounts of alpha/num_letters\n"
       << "\t                                   for an order 0 model.\n"
       << "\t--print-nc-probs, -s - print to files <tx>_true_ncProbs.txt, <tx>_false_ncProbs.txt, where <tx> are all taxons present in reference data,\n"
       << "\t                       normalized conditional probabilities for tuning threshold values of taxon assignment\n"
       << "\t-v                   - verbose mode\n\n"
       << "\t-h|--help            - this message\n\n"

       << "\tConditional probabitity tables are store in\n"
       << "\t<file_i>.MC<order>.log10cProb\n\n"

       << "\n\tOutput file format:\n\n"

       << "\tseqId   model1        model2 ...\n"
       << "\tseq_1   log10prob11   log10prob12  ...\n"
       << "\tseq_2   log10prob21   log10prob22  ...\n"
       << "\t...\n"
       << "\n\twhere log10prob_ij is log10 of the prob(seq_i | model_j)\n\n"

       << "\n\tExample: \n"

       << s << " -d vaginal_v2_MCdir -r vaginal_v2_dir/refTx.tree -f vaginal_v2.fullTx -i vaginal_v2.1.fa -o testDir" << endl << endl
       << s << " -d vaginal_v2_MCdir -r vaginal_v2_dir/refTx.tree -i vaginal_v2.1.fa -o testDir" << endl << endl
       << s << " -e 2BVBACT-97 -t vaginal_v2_dir/spp_paths.txt -k 8 -r vaginal_sppCondensed_v2i.tree -o testDir" << endl << endl;
}


//----------------------------------------------- printHelp ----
void printHelp( const char *s )
{
    cout << endl
	 << "Given a directory of MC model files, reference tree and a fasta file of query sequences,\n"
	 << "classify each sequence of the fasta file to a taxonomic rank corresponding to model\n"
	 << "with the highest probability given that the | log( p(x | M_L) / p(x | M_R) | > thld \n\n";

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
class inPar2_t
{
public:
  inPar2_t();
  ~inPar2_t();

  char *outDir;             /// output directory for MC taxonomy files
  char *mcDir;              /// input directory for MC model files
  char *faDir;              /// directory of reference fasta files
  char *coreErrRFile;       /// core clError R file
  char *trgFile;            /// file containing paths to fasta training files
  char *fullTxFile;         /// fullTx file for printing classification ouput in a long format as in RDP classifier's fixrank
  char *inFile;             /// input file with path(s) to fasta file(s) containing sequences
                            /// for which -log10(prob(seq | model_i)) are to be computed
  char *seqID;              /// sequence ID of a sequence from the training fasta files that is to be excluded
                            /// from model building and needs to be used for cross validation
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
  bool printNCprobs;        /// if true, the program prints to files normalized conditional probabilities for tuning threshold values of taxon assignment
  int dimProbs;             /// max dimension of probs
  bool revComp;             /// reverse-complement query sequences before processing

  void print();
};

//------------------------------------------------- constructor ----
inPar2_t::inPar2_t()
{
  outDir          = NULL;
  mcDir           = NULL;
  trgFile         = NULL;
  fullTxFile      = NULL;
  inFile          = NULL;
  treeFile        = NULL;
  seqID           = NULL;
  thld            = 0.0;
  printCounts     = 0;
  skipErrThld     = 0;
  maxNumAmbCodes  = 5;
  randSampleSize  = 0;
  pseudoCountType = recPdoCount;
  dimProbs        = 0;
  printNCprobs    = false;
  verbose         = false;
  revComp         = false;
}

//------------------------------------------------- constructor ----
inPar2_t::~inPar2_t()
{
  if ( outDir )
    free(outDir);

  if ( mcDir )
    free(mcDir);

  if ( trgFile )
    free(trgFile);

  if ( inFile )
    free(inFile);

  if ( fullTxFile )
    free(fullTxFile);

  if ( seqID )
    free(seqID);

  if ( treeFile )
    free(treeFile);

  int n = trgFiles.size();
  for ( int i = 0; i < n; ++i )
    free(trgFiles[i]);
}

//------------------------------------------------------- print ----
void inPar2_t::print()
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

  cerr << "fullTxFile=\t\t";
  if ( fullTxFile )
    cerr << fullTxFile << endl;
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
void parseArgs( int argc, char ** argv, inPar2_t *p );
bool dComp (double i, double j) { return (i>j); }

//============================== main ======================================
int main(int argc, char **argv)
{
  #define DEBUGMAIN 0
  #define DEBUGMAIN1 0

  struct timeval  tvStart, tvCurrent;
  gettimeofday(&tvStart, NULL);

  //-- setting up init parameters
  inPar2_t *inPar = new inPar2_t();

  //-- parsing input parameters
  parseArgs(argc, argv, inPar);

  if ( inPar->verbose )
    inPar->print();

  if ( inPar->trgFile )
  {
    readLines(inPar->trgFile, inPar->trgFiles); // path(s) from inPar->trgFile are loaded into inPar->trgFiles
    free( inPar->trgFile );
  }

  if ( !inPar->inFile && !inPar->seqID )
  {
    cout << endl << "ERROR in "<< __FILE__ << " at line " << __LINE__ << ": Input fasta file is missing. Please specify it with the -i flag." << endl;
    printHelp(argv[0]);
    exit(1);
  }

  if ( inPar->outDir )
  {
    string cmd("mkdir -p ");
    cmd += string(inPar->outDir);
    system(cmd.c_str());
  }
  else
  {
    cout << endl << "ERROR in "<< __FILE__ << " at line " << __LINE__ << ": Output directory is missing. Please specify it with the -o flag." << endl;
    printHelp(argv[0]);
    exit(1);
  }


  // creating spID => vector of higher taxonomic ranks tbl
  // Ex
  // BVAB1	g_Shuttleworthia	f_Lachnospiraceae	o_Clostridiales	c_Clostridia	p_Firmicutes	d_Bacteria
  // corresponds to
  // fullTx[BVAB1] = (g_Shuttleworthia, f_Lachnospiraceae, o_Clostridiales, c_Clostridia, p_Firmicutes, d_Bacteria)
  map<string, vector<string> > fullTx;

  if ( inPar->fullTxFile )
  {
    char ***tbl;
    int nRows, nCols;
    readCharTbl( inPar->fullTxFile, &tbl, &nRows, &nCols );

    for ( int i = 0; i < nRows; ++i )
      for ( int j = 1; j < nCols; ++j )
      {
	fullTx[ string(tbl[i][0]) ].push_back( string(tbl[i][j]) );
      }
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
    // if ( !in )
    // {
    //   cerr << "Cannot read model ids in " << __FILE__ << " at line " << __LINE__ << endl;
    //   exit(1);
    // }
    // fclose(in);

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
	    fprintf(stderr, "ERROR in %s at line %d: errTbl should have two columns and ncol=%d", __FILE__, __LINE__, ncol);
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

	  // fprintf(stderr, "%s: \n", modelIds[i]);
	  // printDblTbl(errTbl, nrow, ncol);
	  // fprintf(stderr, "\nthld=%f\n", errObj->thld);
	  // exit(1);
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
    fclose(in);
  }
  else if ( !inPar->mcDir && !inPar->trgFile )
  {
    cout << endl << "ERROR in "<< __FILE__ << " at line " << __LINE__ << ": Please specify a directory with MC model files using -d flag." << endl;
    printHelp(argv[0]);
    exit(1);
  }

  NewickTree_t nt;
  if ( inPar->treeFile ) // load ref tree
  {
    if ( !nt.loadTree(inPar->treeFile) )
    {
      fprintf(stderr,"Could not load Newick tree from %s\n", inPar->treeFile);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    // lets see if we can find ref tree in mcDir
    // refTx.tree
    string trFile = string(inPar->mcDir) + "/refTx.tree";
    STRDUP(inPar->treeFile, trFile.c_str());

    if ( !nt.loadTree(inPar->treeFile) )
    {
      cout << endl << "ERROR in "<< __FILE__ << " at line " << __LINE__ << ": reference tree Newick format file is missing. Please specify it with the -r flag." << endl;
      printHelp(argv[0]);
      exit(1);
    }
  }


  int depth = nt.getDepth();
  cerr << "--- Depth of the reference tree: " << depth << endl;


  if ( inPar->kMerLens.size() == 0 )
  {
    int kMers[] = {3};
    cerr << endl << "WARNING: Setting k-mer size to " << kMers[0] << endl;
    int n = sizeof(kMers) / sizeof(int);
    for ( int i = 0; i < n; ++i )
      inPar->kMerLens.push_back(kMers[i]);
  }

  cerr << "--- Number of Models: " << nModels << endl;

  int wordLen = inPar->kMerLens[0];

  if ( inPar->verbose )
    cerr << "\rk=" << wordLen << "\n";

  if ( inPar->mcDir && !inPar->trgFiles.size() )
    cerr << "\r--- Reading conditional probabilities tables from " << inPar->mcDir << " ... ";
  else
    cerr << "\r--- Generating k-mer frequency tables for k=1:" << wordLen << " ... ";


  MarkovChains2_t *probModel;
  probModel = new MarkovChains2_t( wordLen-1,
				   inPar->trgFiles,
				   inPar->mcDir,
				   inPar->maxNumAmbCodes,
				   inPar->pseudoCountType );
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

  char str[10];
  sprintf(str,"%d",(wordLen-1));

  // ==== computing probabilities of each sequence of inFile to come from each of the MC models ====
  string outFile = string(inPar->outDir) + string("/") + string("MC_order") + string(str) + string("_results.txt");

  #if DEBUGMAIN
  string debugFile = string(inPar->outDir) + string("/") + string("debug_log.txt");
  FILE *debugout = fOpen(debugFile.c_str(), "w");
  #endif

  #if DEBUGMAIN1
  FILE *debugout = stderr;
  #endif

  // double *aLogOdds;
  // MALLOC(aLogOdds, double*, depth * sizeof(double));

  vector<string> path; // decision path tranced from the root to the final node for each query sequence
  char *id;
  int count = 0;
  //int depthCount;
  string currentLabel;
  double x[nModels]; // stores conditional probabilities p(x | M) for children of each node. the root node has 3 children

  FILE *out = fOpen(outFile.c_str(), "w");
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

  map<string, vector<double> > txTrueNCProb;  // hash table assigning to each
                                              // taxon a vector of normalized
                                              // conditional probabilities that
                                              // quary sequences achieve using
                                              // max p(x|M) algorithm =
                                              // classification to the taxon with
                                              // the highest posterior
                                              // probability

  map<string, vector<double> > txFalseNCProb; // normalized conditional
                                              // probabilities in cases when
                                              // the taxon does not have the
                                              // highest normalized conditional
                                              // probability

  map<string, vector<char *> > txTrueID;      // hash table assigning to each
                                              // taxon a vector of sequence IDs
                                              // for sequnces with p(x|M) max
                                              // for M the model associated
                                              // with the given taxon
  map<string, vector<char *> > txFalseID;     // same as above but for other sequences


  //cerr << "nRecs=" << nRecs << "\tq01=" << q01 << endl;
  cerr << "--- Number of sequences in " << inPar->inFile << ": " << nRecs << endl;

  //int rcseqCount = 0; // number of times rcseq had higher probabitity than seq
  //int seqCount = 0;   // number of times seq had higher probabitity than rcseq

  int currentModelIdx = 0; // model index of the model, M, with the highest p( x | M )
  int rank = probModel->order() + 1;

  FILE *probsOut = NULL;
  if ( inPar->dimProbs )
  {
    string probsFile = string(inPar->outDir) + string("/") + string("condProbs.csv");
    probsOut = fOpen(probsFile.c_str(), "w");
    inPar->dimProbs = inPar->dimProbs - rank;
  }

  int runTime;
  int timeMin = 0;
  int timeSec = 0;
  int perc;

  while ( getNextFastaRecord( in, id, data, alloc, seq, seqLen) )
  {
    if ( q01 && (count % q01) == 0 )
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

    // traverse the reference tree at each node making a choice of a model
    // and checking log odds of the best model, M, against 'not-M' model

    NewickNode_t *node = nt.root();
    int numChildren = node->children_m.size();
    //path.clear();
    //path.push_back( node->label );
    double err = 0;
    int breakLoop = 0;


#if DEBUGMAIN
    fprintf(debugout,"---- depth %d\n",depthCount++) ;
    for ( int i = 0; i < numChildren; i++ )
      //fprintf(debugout,"\t%s\t%f\t%f\n", node->children_m[i]->label.c_str(), x[i], x2[i]) ;
      fprintf(debugout,"\t%s\t%f\n", node->children_m[i]->label.c_str(), x[i]) ;
#endif

#if DEBUGMAIN1
    fprintf(debugout,"\n---- Processing %s\n",id) ;
    fprintf(debugout,"---- Current node %s\n", node->label.c_str()) ;
    fprintf(debugout,"---- Number of children: %d\n", numChildren) ;
    fprintf(debugout,"---- Children:\n") ;
    for ( int i = 0; i < numChildren; i++ )
	fprintf(debugout,"\t%s\n", node->children_m[i]->label.c_str()) ;
#endif

    //score.clear();
    //int depthCount = 1;
    while ( numChildren && !breakLoop )
    {
      // compute model probabilities for seq and rcseq
      // NOTE: after a few iterations only seq or rcseq should be processed !!!
      for ( int i = 0; i < numChildren; i++ )
      {
	#if 0
	double x1 = probModel->normLog10prob(seq, seqLen, (node->children_m[i])->model_idx );
	double x2 = probModel->normLog10prob(rcseq, seqLen, (node->children_m[i])->model_idx );
	x[i] = ( x1 > x2 ) ? x1 : x2;
	if ( x2 > x1 ) rcseqCount++; else seqCount++;
	#endif

	if ( inPar->revComp )
	{
	  x[i] = probModel->normLog10prob(rcseq, seqLen, (node->children_m[i])->model_idx );
	}
	else
	{
	  x[i] = probModel->normLog10prob(seq, seqLen, (node->children_m[i])->model_idx );
	}

	#if DEBUGMAIN1
	errTbl_t *errObj = modelErrTbl[ (node->children_m[i])->label ];
	fprintf(debugout,"---- Evaluating MC model of %s\tLogPostProb: %f\terrThld: %f\n", (node->children_m[i])->label.c_str(), x[i], errObj->thld);
	//fprintf(debugout,"---- Evaluating MC model of %s\tclErr: %f\tOrientation used: %s\n", (node->children_m[i])->label.c_str(), x[i], ( x1 > x2 ) ? "seq" : "rcSeq") ;
	#endif
      }

      int imax = which_max( x, numChildren );
      currentModelIdx = (node->children_m[imax])->model_idx;

      if ( inPar->printNCprobs )
      {
	txTrueNCProb[node->children_m[imax]->label].push_back(x[imax]);
	txTrueID[node->children_m[imax]->label].push_back(id);
	for ( int i = 0; i < numChildren; i++ )
	  if ( i != imax )
	  {
	    txFalseNCProb[node->children_m[i]->label].push_back(x[i]);
	    txFalseID[node->children_m[i]->label].push_back(id);
	  }

	string file1 = string(inPar->outDir) + string("/") + node->children_m[imax]->label + string("_true_seq.ids");
	FILE *out1 = fOpen(file1.c_str(), "a");
	fprintf(out1, "%s\n", id);
	fclose(out1);

	for ( int i = 0; i < numChildren; i++ )
	  if ( i != imax )
	  {
	    string file2 = string(inPar->outDir) + string("/") + node->children_m[i]->label + string("_true_seq.ids");
	    FILE *out2 = fOpen(file2.c_str(), "a");
	    fprintf(out2, "%s\n", id);
	    fclose(out2);
	  }
      }

      node = node->children_m[imax];

      if ( !inPar->skipErrThld )
      {
	errTbl_t *errObj = modelErrTbl[ node->label ];

	if ( x[imax] > errObj->thld )
	{
	  if ( x[imax] > errObj->xmax )
	  {
	    err = 0;
	  }
	  else
	  {
	    int ierr = bsearchDbl( errObj->x, errObj->nrow, x[imax] );
	    err = errObj->errTbl[ierr][1];
	  }
	}
	else
	{
	  node = node->parent_m;
	  breakLoop = 1;

	  if ( node->label=="d_Bacteria" )
	  break;
	}

        #if DEBUGMAIN1
	fprintf(debugout,"xmax: %s\tLogPostProb=%f\tthld=%f\txmax=%f\terr=%f\n",
		node->label.c_str(), x[imax], errObj->thld, errObj->xmax, err);
        #endif


	//score.push_back( tx2score );
      }

      // tx2score.first = node->label;
      // tx2score.second = err;

      //currentLabel = node->label;
      //path.push_back( node->label );
      numChildren = node->children_m.size();
    }

    fprintf(out,"%s\t%s\t%.4f\n", id, node->label.c_str(), err);
    //fprintf(out,"%s\t%s\n", id, node->label.c_str());
    //fprintf(out,"%s\t%s\t%.2f\n", id, tx2score.first.c_str(), tx2score.second);


    // -----------------------------------------
    // Printing conditional probabilities
    // -----------------------------------------
    if ( inPar->dimProbs )
    {
      // probs = vector of conditional probabilities at each position of the sequence, rcseq, given the modelIdx-th model
      int k = probModel->log10probVect( rcseq, seqLen, currentModelIdx, probs );

      if ( k > inPar->dimProbs )
	k = inPar->dimProbs;

      int k1 = k-1;
      fprintf(probsOut, "%s,", id);
      for ( int i = 0; i < k1; i++ )
	fprintf(probsOut, "%f,", pow(10, probs[i]));
      fprintf(probsOut, "%.10f", pow(10, probs[k1]));

      if ( k < inPar->dimProbs )
	for ( int i = k; i < inPar->dimProbs; i++ )
	  fprintf(probsOut, ",0");

      fprintf(probsOut, "\n");

      // fprintf(stderr, "\n\nk=%d\tdimProbs=%d\n", k, inPar->dimProbs);
      // break;
    }


  } // end of   while ( getNextFastaRecord( in, id, data, alloc, seq, seqLen) )

  if ( inPar->printNCprobs )
  {
    map<string, vector<double> >::iterator it;
    for ( it = txTrueNCProb.begin(); it != txTrueNCProb.end(); it++ )
    {
      string file1 = string(inPar->outDir) + string("/") + it->first + string("_true_ncProbs.txt");
      FILE *out1 = fOpen(file1.c_str(), "w");
      int n = it->second.size();
      for ( int i = 0; i < n; i++ )
	fprintf(out1, "%f\n", it->second[i]);
      fclose(out1);
    }

    for ( it = txFalseNCProb.begin(); it != txFalseNCProb.end(); it++ )
    {
      string file1 = string(inPar->outDir) + string("/") + it->first + string("_false_ncProbs.txt");
      FILE *out1 = fOpen(file1.c_str(), "w");
      int n = it->second.size();
      for ( int i = 0; i < n; i++ )
	fprintf(out1, "%f\n", it->second[i]);
      fclose(out1);
    }

    #if 0
    map<string, vector<char *> >::iterator itr;
    for ( itr = txTrueID.begin(); itr != txTrueID.end(); itr++ )
    {
      string file1 = string(inPar->outDir) + string("/") + itr->first + string("_true_seq.ids");
      FILE *out1 = fOpen(file1.c_str(), "w");
      int n = itr->second.size();
      for ( int i = 0; i < n; i++ )
	fprintf(out1, "%s\n", itr->second[i]);
      fclose(out1);
    }

    for ( itr = txFalseID.begin(); itr != txFalseID.end(); itr++ )
    {
      string file1 = string(inPar->outDir) + string("/") + itr->first + string("_false_seq.ids");
      FILE *out1 = fOpen(file1.c_str(), "w");
      int n = itr->second.size();
      for ( int i = 0; i < n; i++ )
	fprintf(out1, "%s\n", itr->second[i]);
      fclose(out1);
    }
    #endif
  }


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
  fprintf(stderr,"\r                                                                       \n");
  fprintf(stderr,"    Elapsed time: %d:%02d                                              \n", timeMin, timeSec);

  // fprintf(stderr,"\r--- Number of processed sequences: %d                                  \n", count);
  // fprintf(stderr,"    Number of times rcseq had higher probabitity than seq: %d\n", rcseqCount);
  // fprintf(stderr,"    Number of times rcseq had lower probabitity than seq: %d\n", seqCount);
  //fprintf(stderr,"Output written to %s\n", outFile.c_str());
  fprintf(stderr,"    Output written to %s\n", inPar->outDir);

  fprintf(stderr,"\n    To create a sample x phylotype count table, run\n");
  fprintf(stderr,"\n        count_tbl.pl -i %s -o %s/spp_count_tbl.txt\n\n", outFile.c_str(), inPar->outDir);

  #if DEBUGMAIN
  fprintf(stderr,"\n\nDEBUGING Output written to %s\n\n\n", debugFile.c_str());
  #endif

  return EXIT_SUCCESS;
}



//----------------------------------------------------------- parseArgs ----
//! parse command line arguments
void parseArgs( int argc, char ** argv, inPar2_t *p )
{
  int c, errflg = 0;
  optarg = NULL;

  static struct option longOptions[] = {
    {"print-counts"       ,no_argument, &p->printCounts,    1},
    {"skip-err-thld"      ,no_argument, &p->skipErrThld,    1},
    {"max-num-amb-codes"  ,required_argument, 0,          'b'},
    {"fullTx-file"        ,required_argument, 0,          'f'},
    {"out-dir"            ,required_argument, 0,          'o'},
    {"ref-tree"           ,required_argument, 0,          'r'},
    {"pseudo-count-type"  ,required_argument, 0,          'p'},
    {"print-nc-probs"     ,no_argument, 0,                's'},
    {"rev-comp"           ,no_argument, 0,                'c'},
    {"help"               ,no_argument, 0,                  0},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv,"a:b:c:d:e:f:g:t:i:k:o:vp:r:hsy:",longOptions, NULL)) != -1)
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

      case 'y':
	p->coreErrRFile = strdup(optarg);
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

      case 'o':
	p->outDir = strdup(optarg);
	break;

      case 't':
	p->trgFile = strdup(optarg);
	break;

      case 'f':
	p->fullTxFile = strdup(optarg);
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

      case 's':
	p->printNCprobs = true;
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
