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
       << " Using prebuilt MC models to classify sequences of an input fasta file" << endl
       << endl
       << s << " -d < MC models directory> -i <input fasta file> -o <output directory> [Options]" << endl
       << endl
       << "\tOptions:\n"
       << "\t-d <mcDir>      - directory containing MC model files\n"
       << "\t-o <outDir>     - output directory for MC taxonomy files\n"
       << "\t-i <inFile>     - input fasta file with sequences for which -log10(prob(seq | model_i)) are to be computed\n"
       << "\t-r <model tree> - model tree with node labels corresponding to the names of the model files\n"
       << "\t-t <trgFile>    - file containing paths to training fasta files\n"
       << "\t-f <fullTx>     - fullTx file. Its optional parameter for printing classification output in a long format like in RDP classifier\n"
       << "\t-g <faDir>      - directory with reference fasta files\n"
       << "\t--rev-comp, -c          - reverse complement query sequences before computing classification posterior probabilities\n"
       << "\t--skip-err-thld         - classify all sequences to the species level\n"
       << "\t--pp-embedding          - for each internal node report pp's of all children on the given sequence.\n"
       << "                            Each internal node's table is written to a file <node name>_ref_lpps.txt (log posterior probabilities)\n"
       << "\t--max-num-amb-codes <n> - maximal acceptable number of ambiguity codes for a sequence\n"
       << "\t                          above this number sequence's log10prob() is not computed and\n"
       << "\t                          the sequence's id it appended to <genus>_more_than_<n>_amb_codes_reads.txt file.\n"
       << "\t                          Default value: 5\n\n"
       << "\t--pseudo-count-type, -p <f>  - f=0 for add 1 to all k-mer counts zero-offset\n"
       << "\t                               f=1 for add 1/4^k to k-mer counts zero-offset\n"
       << "\t                               f=2 the pseudocounts for a order k+1 model be alpha*probabilities from\n"
       << "\t                                   an order k model, recursively down to pseudocounts of alpha/num_letters\n"
       << "\t                                   for an order 0 model.\n"
       << "\t-q|--quiet           - suppers pregress messages\n"
       << "\t-v|--verbose         - verbose mode\n"
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

       << s << " -d vaginal_v2_MCdir -f vaginal_v2.fullTx -i vaginal_v2.1.fa -o testDir" << endl << endl
       << s << " -d vaginal_v2_MCdir -r vaginal_v2_dir/model.tree -i vaginal_v2.1.fa -o testDir" << endl << endl
       << s << " -e 2BVBACT-97 -t vaginal_v2_dir/tx_fasta_paths.txt -k 8 -r vaginal_sppCondensed_v2i.tree -o testDir" << endl << endl;
}


//----------------------------------------------- printHelp ----
void printHelp( const char *s )
{
    cout << endl
	 << "Given a fasta file of query sequences, a directory of MC model files and the reference tree\n"
	 << "classify each sequence of the fasta file to a taxonomic rank corresponding to model\n"
	 << "with the highest probability given that the | log( p(x | M_L) / p(x | M_R) | > thld (obsolete) \n\n";

    printUsage(s);
}

//----------------------------------------------- errTbl_t ----
/// holds errTbl and log10cProb.thld
typedef struct
{
  double **errTbl; ///
  int    nrow;
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

  char *outDir;             /// output directory for MC taxonomy files
  char *mcDir;              /// input directory for MC model files
  char *faDir;              /// directory of reference fasta files
  char *coreErrRFile;       /// core clError R file
  char *trgFile;            /// file containing paths to fasta training files
  char *fullTxFile;         /// fullTx file for printing classification ouput in a long format as in RDP classifier's fixrank
  char *inFile;             /// input file with path(s) to fasta file(s) containing sequences
                            /// for which -log10(prob(seq | model_i)) are to be computed
  char *treeFile;           /// reference tree file
  double thld;              /// threshold for | log( p(x | M_L) / p(x | M_R) | of the competing models
  vector<char *> trgFiles;  /// list of paths to fasta training files
  int kMerLen;              /// max k-mer length in MC models
  int printCounts;          /// flag initiating print out of word counts
  int skipErrThld;          /// ignore classification error condition - with this option on each sequence is classified to the species with the highest p(x|M)
  int maxNumAmbCodes;       /// maximal acceptable number of ambiguity codes for a sequence; above this number log10probIUPAC() returns 1;
  int randSampleSize;       /// number of random sequences of each model (seq length = mean ref seq). If 0, no random samples will be generated.
  int pseudoCountType;      /// pseudo-count type; see MarkovChains.hh for possible values
  bool verbose;
  bool quiet;
  bool revComp;             /// reverse-complement query sequences before processing
  bool ppEmbedding;

  void print();
};

//------------------------------------------------- constructor ----
inPar_t::inPar_t()
{
  outDir          = NULL;
  mcDir           = NULL;
  trgFile         = NULL;
  fullTxFile      = NULL;
  inFile          = NULL;
  treeFile        = NULL;
  kMerLen         = 8;
  thld            = 0.0;
  printCounts     = 0;
  skipErrThld     = 0;
  maxNumAmbCodes  = 5;
  randSampleSize  = 0;
  pseudoCountType = recPdoCount;
  verbose         = false;
  quiet           = false;
  revComp         = false;
  ppEmbedding     = false;
}

//------------------------------------------------- constructor ----
inPar_t::~inPar_t()
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

  cerr << "kMerLen: " << kMerLen << endl;

  cerr << "skipErrThld: " << skipErrThld << endl;
}

//============================== local sub-routines =========================
void parseArgs( int argc, char ** argv, inPar_t *p );
bool dComp (double i, double j) { return (i>j); }

//============================== main ======================================
int main(int argc, char **argv)
{
    #define DEBUG_CLASSIFY 0

    struct timeval  tvStart, tvCurrent;
    gettimeofday(&tvStart, NULL);

    //-- setting up init parameters
    inPar_t *inPar = new inPar_t();

    //-- parsing input parameters
    parseArgs(argc, argv, inPar);

    if ( inPar->verbose )
      inPar->print();

    #if 0
    if ( inPar->trgFile )
    {
      readLines(inPar->trgFile, inPar->trgFiles); // path(s) from inPar->trgFile are loaded into inPar->trgFiles
      free( inPar->trgFile );
    }
    #endif

    if ( !inPar->inFile )
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

    //map<string, errTbl_t *> modelErrTbl;
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
          // Reading error thresholds file
          string file = string(inPar->mcDir) + string("/error_thlds.txt");
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

        if ( k != inPar->kMerLen )
        {
          fprintf(stderr,"WARNING Resetting kMerLen to %d\n", k);
          inPar->kMerLen = k;
        }
      } // END OF if ( in ); READING modelIds.txt, error_thlds.txt
      fclose(in);
    }
    else if ( !inPar->mcDir && !inPar->trgFile )
    {
      cout << endl << "ERROR in "<< __FILE__ << " at line " << __LINE__ << ": Please specify a directory with MC model files using -d flag." << endl;
      printHelp(argv[0]);
      exit(1);
    }

    if ( inPar->verbose )
      fprintf(stderr,"--- Number of Models: %d\n", nModels);

    NewickTree_t *nt;
    if ( inPar->treeFile ) // load ref tree
    {
      nt = readNewickTree( inPar->treeFile );
      if ( !nt )
      {
        fprintf(stderr,"ERROR in %s at line %d: Could not load Newick tree from %s\n", __FILE__, __LINE__, inPar->treeFile);
        exit(EXIT_FAILURE);
      }
    }
    else
    {
      // Lets see if we can find a model tree in mcDir
      string trFile = string(inPar->mcDir) + "/model.tree";
      STRDUP(inPar->treeFile, trFile.c_str());

      nt = readNewickTree( inPar->treeFile );
      if ( !nt )
      {
        fprintf(stderr,"ERROR in %s at line %d: Could not load Newick tree from %s\n", __FILE__, __LINE__, inPar->treeFile);
        exit(EXIT_FAILURE);
      }
    }

    int depth = nt->getDepth();
    if ( inPar->verbose )
      fprintf(stderr,"--- Depth of the model tree: %d\n", depth);


    if ( inPar->verbose )
    {
      fprintf(stderr,"\r--- Rank of MC models: %d\n", inPar->kMerLen);
      fprintf(stderr,"\r--- Reading conditional probabilities tables from %s\n", inPar->mcDir);
    }
    MarkovChains_t *probModel = new MarkovChains_t(inPar->kMerLen-1,
                                                   inPar->mcDir,
                                                   inPar->maxNumAmbCodes,
                                                   inPar->pseudoCountType,
                                                   inPar->verbose );
    if ( inPar->verbose )
      fprintf(stderr,"DONE\n");

    vector<char *> modelIds = probModel->modelIds();
    vector<string> modelStrIds;
    probModel->modelIds( modelStrIds );

    nt->modelIdx( modelStrIds );

    // ==== computing probabilities of each sequence of inFile to come from each of the MC models ====

    // the output (classification results) file
    char str[10];
    sprintf(str,"%d",(inPar->kMerLen-1));
    string outFile = string(inPar->outDir) + string("/") + string("MC_order") + string(str) + string("_results.txt");
    FILE *out      = fOpen(outFile.c_str(), "w");

    // input fasta file
    FILE *in = fOpen(inPar->inFile, "r");

    // This is used only to report progress of the main loop every q01 number of sequences
    //int count = 0;
    int nRecs = numRecordsInFasta( inPar->inFile );
    //int q01 = 0;
    //if ( nRecs > 1000 )
    //  q01 = int(0.01 * nRecs);

    // This is for reading fasta records
    char *id;
    int seqLen;
    size_t alloc = 1024*1024;
    char *data, *seq, *rcseq;
    MALLOC(data, char*, alloc * sizeof(char));
    MALLOC(seq, char*, alloc * sizeof(char));

    // This is used only when sequences are reverse complemented
    if ( inPar->revComp )
      MALLOC(rcseq, char*, alloc * sizeof(char));

    //double lpp[nModels]; // stores conditional log10 posterior probabilities p(x | M) for children of each node. The root node has 3 children
    double *lpp;
    MALLOC(lpp, double*, nModels * sizeof(double));

    if ( inPar->verbose )
      cerr << "--- Number of sequences in " << inPar->inFile << ": " << nRecs << endl;

    // // this is only for time reporting (should be inside some el_time() routine )
    // int runTime;
    // int timeMin = 0;
    // int timeSec = 0;
    // int perc;

    while ( getNextFastaRecord(in, id, data, alloc, seq, seqLen) )
    {
      if ( inPar->revComp )
      {
        for ( int j = 0; j < seqLen; ++j )
          rcseq[j] = Complement(seq[seqLen-1-j]);
        rcseq[seqLen] = '\0';
      }

      #if DEBUG_CLASSIFY
      fprintf(stderr,"\ncount: %d  seqID: %s\n", count, id);
      //fprintf(stderr,"seq: %s\n", seq);
      #endif

      //
      // Traversing the model tree
      //
      // At each node select a model with the highest posterior probability of
      // the sequence coming from the model, given the posterior probability is
      // above the error threshold of the model.
      //
      NewickNode_t *node = nt->root();

      double finalPP = 0.0; // posterior probability of the sequence w/r to the winner model
      int nDecissions = 0;  // for each sequence count the number of decisions the
                            // classifier has to make to get to the final node.

      #if DEBUG_CLASSIFY
      fprintf(stderr,"Before entering model tree\nNumber of root's children: %d\n", nChildren);
      #endif

      //int depthCount = 1;
      while ( node->idx < 0 ) // internal node
      {
        int nChildren = node->children_m.size();
        nDecissions += nChildren;

        // Computing model log10 posterior probabilities for seq (or rcseq)
        for ( int i = 0; i < nChildren; i++ )
        {
          #if DEBUG_CLASSIFY
          fprintf(stderr, "Processing child %d\t%s\n", i, (node->children_m[i])->label.c_str());
          //fprintf(stderr, "seq: %s\n", seq);
          fprintf(stderr, "model_idx: %d\n", (node->children_m[i])->model_idx);
          #endif

          if ( !inPar->revComp )
          {
            lpp[i] = probModel->normLog10prob(seq, seqLen, (node->children_m[i])->model_idx );
          }
          else
          {
            lpp[i] = probModel->normLog10prob(rcseq, seqLen, (node->children_m[i])->model_idx );
          }

          #if DEBUG_CLASSIFY
          fprintf(stderr, "lpp[%d]: %f\n", i, lpp[i] );
          #endif
        }

        int imax = which_max(lpp, nChildren);
        finalPP = pow(10, lpp[imax]);
        node = node->children_m[imax];

        #if DEBUG_CLASSIFY
        fprintf(stderr,"maxModel: %s\tlpp: %f\tthld: %.4f\n",
                node->label.c_str(), lpp[imax], thldTbl[ node->label ]);
        #endif

        if ( !inPar->skipErrThld && nChildren==0 && lpp[imax] < thldTbl[ node->label ] )
        {
          #if 0
          fprintf(stderr,"\n---- Processing %s\n",id) ;
          fprintf(stderr,"maxModel: %s\tlpp: %f\tthld: %.4f\t",
                  node->label.c_str(), lpp[imax], thldTbl[ node->label ]);

          map<string, string>::iterator itr;
          map<string, string> seqRecs; // fasta file sequence records

          string faFile = string(inPar->mcDir) + string("/fasta_files/") + node->label + string(".fa");
          seqRecs.clear();
          readFasta( faFile.c_str(), seqRecs);
          int nRefSeqs = seqRecs.size();

          vector<double> lpp(nRefSeqs);
          double minLpp = 1.0;
          int i = 0;
          for ( itr = seqRecs.begin(); itr != seqRecs.end(); ++itr )
          {
            lpp[i] = probModel->normLog10prob(itr->second.c_str(), (int)itr->second.size(), node->model_idx );
            if ( lpp[i] < minLpp )
              minLpp = lpp[i];
            i++;
          }

          fprintf(stderr,"min lpp: %f\n", minLpp);
          #endif

          node = node->parent_m;
          //breakLoop = 1;

          if ( node->label=="d_Bacteria" )
            break;

        } // END OF if ( !inPar->skipErrThld && node->children_m.size()==0 && lpp[imax] < thldTbl[ node->label ] )

      } // END OF while ( node->idx < 0 )


    fprintf(out,"%s\t%s\t%f\t%d\n", id, node->label.c_str(), finalPP, nDecissions);

    #if DEBUG_CLASSIFY
    fprintf(stderr,"%s\t%s\t%f\t%d\n", id, node->label.c_str(), finalPP, nDecissions);
    exit(0);
    #endif

  } // END OF  while ( getNextFastaRecord( in, id, data, alloc, seq, seqLen) )

  fclose(in);
  fclose(out);

  free(lpp);
  free(seq);
  free(data);

  if ( inPar->revComp )
    free(rcseq);

  // It may be a nice idea to report the number of species found

  gettimeofday(&tvCurrent, NULL);
  //runTime = tvCurrent.tv_sec  - tvStart.tv_sec;

  #if 0
  if ( runTime > 60 )
  {
    timeMin = runTime / 60;
    timeSec = runTime % 60;
  }
  else
  {
    timeSec = runTime;
  }
  #endif

  if ( inPar->verbose )
  {
    //fprintf(stderr,"\r                                                                       \n");
    //fprintf(stderr,"    Elapsed time: %d:%02d                                              \n", timeMin, timeSec);

    // fprintf(stderr,"\r--- Number of processed sequences: %d                                  \n", count);
    // fprintf(stderr,"    Number of times rcseq had higher probabitity than seq: %d\n", rcseqCount);
    // fprintf(stderr,"    Number of times rcseq had lower probabitity than seq: %d\n", seqCount);
    //fprintf(stderr,"Output written to %s\n", outFile.c_str());
    fprintf(stderr,"    Output written to %s\n", inPar->outDir);

    fprintf(stderr,"\n    To create a sample x phylotype count table, run\n");
    fprintf(stderr,"\n        count_tbl.pl -i %s -o %s/spp_count_tbl.txt\n\n", outFile.c_str(), inPar->outDir);
  }

#if 0
  fprintf(stderr,"\n\nDEBUGING Output written to %s\n\n\n", debugFile.c_str());
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
    {"max-num-amb-codes"  ,required_argument, 0, 'b'},
    {"rev-comp"           ,no_argument,       0, 'c'},
    {"pp-embedding"       ,no_argument,       0, 'e'},
    {"fullTx-file"        ,required_argument, 0, 'f'},
    {"out-dir"            ,required_argument, 0, 'o'},
    {"pseudo-count-type"  ,required_argument, 0, 'p'},
    {"ref-tree"           ,required_argument, 0, 'r'},
    {"quiet"              ,no_argument,       0, 'q'},
    {"verbose"            ,no_argument,       0, 'v'},
    {"skip-err-thld"      ,no_argument,       0, 'x'},
    {"help"               ,no_argument,       0,   0},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv,"b:c:d:ef:g:hi:k:o:p:q:r:t:vy:x",longOptions, NULL)) != -1)
    switch (c)
    {
      case 'b':
        p->maxNumAmbCodes = atoi(optarg);
        break;

      case 'c':
        p->revComp = true;
        break;

      case 'd':
        p->mcDir = strdup(optarg);
        break;

      case 'e':
        p->ppEmbedding = true;
        break;

      case 'f':
        p->fullTxFile = strdup(optarg);
        break;

      case 'g':
        p->faDir = strdup(optarg);
        break;

      case 'h':
        printHelp(argv[0]);
        exit (EXIT_SUCCESS);
        break;

      case 'i':
        p->inFile = strdup(optarg);
        break;

      case 'k':
        p->kMerLen = atoi(optarg);
        break;

      case 'o':
        p->outDir = strdup(optarg);
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

      case 'q':
        p->quiet = true;
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

      case 'x':
        p->skipErrThld = 1;
        //cerr << "Setting p->skipErrThld to " << p->skipErrThld << endl;
        break;

      case 'y':
        p->coreErrRFile = strdup(optarg);
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
