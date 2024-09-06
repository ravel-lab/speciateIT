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
#include "CUtilities.h"
#include "CStatUtilities.h"
#include "IOCUtilities.h"
#include "IOCppUtilities.hh"
#include "CppUtilities.hh"
#include "MarkovChains2.hh"
#include "StatUtilities.hh"

using namespace std;

//----------------------------------------------------------- printUsage ----
void printUsage( const char *s )
{
    cout << endl

	 << "USAGE " << endl << endl

	 << "1. Building MC models on the fly" << endl << endl
	 << s << "-t <training file's paths file> -k <k-mer size> -i <input fasta file> -o <output directory> [Options]" << endl << endl

	 << "2. Building MC models only" << endl << endl
	 << s << "-t <training file's paths file> -k <k-mer size> -d < MC models directory> [Options]" << endl << endl

	 << "3. Using prebuilt MC models" << endl << endl
	 << s << "-d < MC models directory> -i <input fasta file> -o <output directory> [Options]" << endl << endl

	 << "\twhere <file_1.fa> ... <file_n.fa> are fasta files to be used for building MC models\n"
	 << endl
	 << "\tOptions:\n"
	 << "\t-d <dir>     - directory for MC model files\n"
	 << "\t-o <dir>     - output directory for MC taxonomy files\n"
	 << "\t-i <inFile>  - input fasta file with sequences for which -log10(prob(seq | model_i)) are to be computed\n"
	 << "\t-t <trgFile> - file containing paths to training fasta files\n"
	 << "\t-k <K>       - K is the k-mer size\n"
	 << "\t--random-sample-size, -r <n> - number of random sequences to be generated for each MC model\n"
	 << "\t--pseudo-count-type, -p <f>  - f=0 for add 1 to all k-mer counts zero-offset\n"
	 << "\t                               f=1 for add 1/4^k to k-mer counts zero-offset\n"
	 << "\t                               f=2 the pseudocounts for a order k+1 model be alpha*probabilities from\n"
	 << "\t                                   an order k model, recursively down to pseudocounts of alpha/num_letters\n"
	 << "\t                                   for an order 0 model.\n"
	 << "\t-v - verbose mode\n\n"
	 << "\t-h|--help      - this message\n\n"

	 // << "\t--max-num-amb-codes <n> - maximal acceptable number of ambiguity codes for a sequence\n"
	 // << "\t                          above this number sequence's log10prob() is not computed and\n"
	 // << "\t                          the sequence's id it appended to <genus>_more_than_<n>_amb_codes_reads.txt file.\n"
	 // << "\t                          Default value: 5\n\n"

	 << "\tConditional probabitity tables are store in\n"
	 << "\t<file_i>.MC<order>.log10cProb\n\n"

	 << "\n\tOutput file format:\n\n"

	 << "\tseqId   model1        model2 ...\n"
	 << "\tseq_1   log10prob11   log10prob12  ...\n"
	 << "\tseq_2   log10prob21   log10prob22  ...\n"
	 << "\t...\n"
	 << "\n\twhere log10prob_ij is log10 of the prob(seq_i | model_j)\n\n"

	 << "\n\tExample: \n"

	 << "1. Building MC models on the fly" << endl << endl
	 << s << " -t vaginal_319F_806R_nr_dir/spp_paths.txt -k 3 -i test.fa -o test_dir" << endl << endl

	 << "2. Building MC models only" << endl << endl
	 << s <<" -t vaginal_319F_806R_nr_dir/spp_paths.txt -k 3 -d vaginal_319F_806R_nr_MCdir" << endl << endl

	 << "3. Using prebuilt MC models" << endl << endl
	 << s << " -d vaginal_319F_806R_nr_MCdir -i test.fa -o test_dir" << endl << endl
	 << endl;
}


//----------------------------------------------- printHelp ----
void printHelp( const char *s )
{
    cout << endl
	 << "Given fasta file(s) of training sequences, a sequence of non-negative integers,\n"
	 << "and a fasta file of query sequences, build Markov chain models for sequences in each fasta file and each order\n\n";

    printUsage(s);
}

// TODO
// 1. replace  vector<int> kMerLens
//    by in kMerSize
// 2. -log10(prob(x|M)) calculation enclose in a method

//================================================= inPar2_t ====
//! holds input parameters
class inPar2_t
{
public:
  inPar2_t();
  ~inPar2_t();

  char *outDir;             /// output directory for MC taxonomy files
  char *mcDir;              /// input directory for MC model files
  char *trgFile;            /// file containing paths to fasta training files
  char *inFile;             /// input file with path(s) to fasta file(s) containing sequences
                            /// for which -log10(prob(seq | model_i)) are to be computed
  vector<char *> trgFiles;  /// list of paths to fasta training files
  vector<int> kMerLens;     /// list of word lengths
  int printCounts;          /// flag initiating print out of word counts
  int maxNumAmbCodes;       /// maximal acceptable number of ambiguity codes for a sequence; above this number log10probIUPAC() returns 1;
  int randSampleSize;       /// number of random sequences of each model (seq length = mean ref seq). If 0, no random samples will be generated.
  int pseudoCountType;      /// pseudo-count type; see MarkovChains2.hh for possible values
  bool verbose;

  void print();
};

//------------------------------------------------- constructor ----
inPar2_t::inPar2_t()
{
  outDir          = NULL;
  mcDir           = NULL;
  trgFile         = NULL;
  inFile          = NULL;
  printCounts     = 0;
  maxNumAmbCodes  = 5;
  randSampleSize  = 0;
  pseudoCountType = recPdoCount;
  verbose         = false;
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

  int n = trgFiles.size();
  for ( int i = 0; i < n; ++i )
    free(trgFiles[i]);
}

//------------------------------------------------------- print ----
void inPar2_t::print()
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

//============================== main ======================================
int main(int argc, char **argv)
{
  //-- setting up init parameters
  inPar2_t *inPar = new inPar2_t();

  //-- parsing input parameters
  parseArgs(argc, argv, inPar);

  if ( inPar->trgFile )
  {
    readLines(inPar->trgFile, inPar->trgFiles); // path(s) from inPar->trgFile are loaded into inPar->trgFiles
    free( inPar->trgFile );
  }

  if ( inPar->verbose )
    inPar->print();

  if ( !inPar->mcDir && !inPar->trgFiles.size() && !inPar->trgFile )
  {
    cout << endl
	 << "ERROR: Please specify either" << endl
	 << "\tA directory with MC model files using -d flag." << endl
	 << "\tor" << endl
	 << "\tA file with path(s) to fasta training file(s) using -t flag." << endl;
    printHelp(argv[0]);
    exit(1);
  }


  if ( inPar->mcDir && inPar->trgFiles.size() && inPar->inFile )
  {
    cout << endl
	 << "ERROR: Please use -d or -t flag but not both at the same time." << endl;
    printHelp(argv[0]);
    exit(1);
  }


  if ( inPar->inFile && !inPar->outDir )
  {
    cout << endl << "ERROR: Output directory is missing. Please specify it with the -o flag." << endl;
    printHelp(argv[0]);
    exit(1);
  }

  if ( inPar->outDir )
  {
    string cmd("mkdir -p ");
    cmd += string(inPar->outDir);
    system(cmd.c_str());
  }

  int nModels = 0;

  if ( inPar->trgFiles.size() )
  {
    nModels = inPar->trgFiles.size();
  }
  else if ( inPar->mcDir ) // extracting number of models and k-mer size
  {
    string inFile(inPar->mcDir);
    inFile += "/modelIds.txt";
    FILE *in = fopen(inFile.c_str(), "r");
    if ( !in )
    {
      cerr << "Cannot read model ids in " << __FILE__ << " at line " << __LINE__ << endl;
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

  if ( inPar->kMerLens.size() == 0 )
  {
    int kMers[] = {3};
    cerr << endl << "WARNING: Setting k-mer size to " << kMers[0] << endl;
    int n = sizeof(kMers) / sizeof(int);
    for ( int i = 0; i < n; ++i )
      inPar->kMerLens.push_back(kMers[i]);
  }

  //int nWordSizes = inPar->kMerLens.size();
  double *probs;
  MALLOC(probs, double*, nModels * sizeof(double));

  char *id, *seq;
  int seqLen;

  cerr << endl << "--- nModels=" << nModels << endl;
  //<< "\tnWordSizes="  << nWordSizes << endl;



  // ==== generating random samples from MC models ====
  if ( inPar->randSampleSize && !( (inPar->mcDir || inPar->trgFiles.size()) && inPar->outDir ) )
  {
    if ( !inPar->mcDir && !inPar->trgFiles.size() )
    {
      cout << endl
	   << "ERROR: When generating random samples from MC models either mcDir or sppFile need to be specified." << endl;
    }
    else
    {
      cout << endl
	   << "ERROR: When generating random samples from MC models output directory needs to be specified." << endl;
    }
    printHelp(argv[0]);
    exit(1);
  }
  else if ( inPar->randSampleSize )
  {
    int wordLen = inPar->kMerLens[0];

    if ( inPar->mcDir && !inPar->trgFiles.size() )
      cerr << "\r--- Reading k-mer frequency tables from " << inPar->mcDir << " ... ";
    else
      cerr << "\r--- Generating k-mer frequency tables for k=1:" << wordLen << " ... ";

    MarkovChains2_t probModel( wordLen-1,
			       inPar->trgFiles,
			       inPar->mcDir,
			       inPar->maxNumAmbCodes,
			       inPar->pseudoCountType);
    cerr << "done" << endl;

    string faFile = string(inPar->outDir) + string("/") + string("rsample.fa");
    string txFile = string(inPar->outDir) + string("/") + string("rsample.tx");

    probModel.sample( faFile.c_str(), txFile.c_str(), inPar->randSampleSize );

    cout << endl << "Random sequences for all MC models written to " << faFile << endl << endl;

    return EXIT_SUCCESS;
  }


  //for ( int kIdx = 0; kIdx < nWordSizes; ++kIdx )
  {
    //cerr << "kIdx=" << kIdx << endl;
    int wordLen = inPar->kMerLens[0];
    //cerr << "wordLen=" << wordLen << endl;

    if ( inPar->verbose )
      cerr << "\rk=" << wordLen << "\n";

    if ( inPar->mcDir && !inPar->trgFiles.size() )
      cerr << "\r--- Reading k-mer frequency tables from " << inPar->mcDir << " ... ";
    else
      cerr << "\r--- Generating k-mer frequency tables for k=1:" << wordLen << " ... ";

    MarkovChains2_t probModel( wordLen-1,
			       inPar->trgFiles,
			       inPar->mcDir,
			       inPar->maxNumAmbCodes,
			       inPar->pseudoCountType );

    vector<char *> modelIds = probModel.modelIds();

    cerr << "done" << endl;


    if ( inPar->verbose )
    {
      cout << "modelIds.size()=" << modelIds.size() << "\nmodelIds:\t";
      printVector(modelIds);

      cout << "words[0]:\t";
      vector<vector<char *> > words = probModel.wordStrgs();
      printVector(words[0]);
    }

    char str[10];
    sprintf(str,"%d",(wordLen-1));

    if ( !inPar->inFile && inPar->mcDir )
    {
      cout << "\nMarkov chain models written to " << inPar->mcDir << endl << endl;
      return EXIT_SUCCESS;
    }


    // ==== computing probabilities of each sequence of inFile to come from each of the MC models ====
    string outFile = string(inPar->outDir) + string("/") + string(baseFileName(inPar->inFile));//string(pathWithoutSuffix(inPar->inFile));
    outFile += string(".MC.order") + string(str);
    outFile += string(".nLog10prob");

    // cerr << "outFile=" << outFile.c_str() << endl;
    // output consists of log10 probabilities; file format
    //
    //   seqId  model1         model2        ... modelN
    //   seq_1  log10(prob11)  log10(prob12) ... log10(prob1N)
    //   seq_2  log10(prob21)  log10(prob22) ... log10(prob2N)
    //   ....

    FILE *out = fOpen(outFile.c_str(), "w");
    //writeHeader( out, modelIds );

    // low quality (more than 7 ambiguity codes) reads are
    // not processed and their ids are printed to inFile_low_quality_reads.txt
    string outFile2 = string(inPar->outDir) + string("/") + (baseFileName(inPar->inFile));
    sprintf(str,"%d",inPar->maxNumAmbCodes);
    outFile2 += string("_more_than_") + string(str) + string("_amb_codes_reads.txt");
    FILE *out2 = fOpen(outFile2.c_str(), "w");
    fclose(out2);
    out2 = fOpen(outFile2.c_str(), "a");

    FILE *in = fOpen(inPar->inFile, "r");

    int nRecs = numRecordsInFasta( inPar->inFile );
    int q01 = 0;
    int count = 0;

    if ( nRecs > 1000 )
      q01 = int(0.01 * nRecs);

    //cerr << "nRecs=" << nRecs << "\tq01=" << q01 << endl;
    cerr << "--- Number of sequences in " << inPar->inFile << " = " << nRecs << endl;

    size_t alloc = 1024*1024;
    char *data, *rcseq;
    MALLOC(data, char*, alloc * sizeof(char));
    MALLOC(seq, char*, alloc * sizeof(char));
    MALLOC(rcseq, char*, alloc * sizeof(char));

    while ( getNextFastaRecord( in, id, data, alloc, seq, seqLen) )
    {
      if ( q01 && (count % q01) == 0 )
      {
	int q = int( (100.0*count) / nRecs);
	cerr << "\r" << "Number of sequences processed = " << count << "\t" << q << "%";
      }

      count++;

      int nAmbCodes = 0;

      for ( int j = 0; j < seqLen; ++j )
	if ( intACGTLookup[int(seq[j])] == -1 )
	  nAmbCodes++;

      if ( nAmbCodes > inPar->maxNumAmbCodes )
      {
	fprintf(out2,"%s\t%d\n",id,nAmbCodes);
	continue;
      }


      for ( int i = 0; i < nModels; ++i )
      {
	double x1 = probModel.log10prob(seq, seqLen, i);

#if 0
	// proces reverse complement as well and pick the one with smaller log10prob()
	for ( int j = 0; j < seqLen; ++j )
	  rcseq[j] = Complement(seq[seqLen-1-j]);
	rcseq[seqLen] = '\0';

	double x2 = probModel.log10prob(rcseq, seqLen, i);
	probs[i] = ( x1 > x2 ) ? x1 : x2;
#endif
	probs[i] = x1;
      }

      int imax = which_max( probs, nModels );
      fprintf(out,"%s\t%s\n", id, modelIds[imax]);
      //writeProbs(out, id, probs, nModels);
    }

    free(seq);
    free(rcseq);
    free(data);

    fclose(in);
    fclose(out);
    fclose(out2);

    cerr << "\r\nOutput written to " << outFile.c_str() << endl;
    cerr << "Low quality read ids written to " << outFile2.c_str() << endl << endl;
  }

  free(probs);

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
    {"max-num-amb-codes"  ,required_argument, 0,          'b'},
    {"out-dir"            ,required_argument, 0,          'o'},
    {"random-sample-size" ,required_argument, 0,          'r'},
    {"pseudo-count-type"  ,required_argument, 0,          'p'},
    {"help"               ,no_argument, 0,                  0},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv,"b:d:t:i:k:o:vp:r:h",longOptions, NULL)) != -1)
    switch (c)
    {
      case 'b':
	p->maxNumAmbCodes = atoi(optarg);
	break;

      case 'r':
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
