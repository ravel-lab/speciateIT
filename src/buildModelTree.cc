//
// buildModelTree.cc
//
// Input
//    - reference .tx file (seqID => sppID)
//    - full taxon file, .fullTx
//    - reference fasta file
//    - output directory (outDir)

// Output
//    - For each node of the ref tree (except root) a fasta file of all ref
//      seq's corresopnding to the node's subtree
//    - A file of absolute paths to just created fasta files
//    - A file of taxonomic assignments of internal nodes

// Given
//    - a species lineage file that has the following structure of a tab delimited table

//  BVAB1	g_Shuttleworthia	f_Lachnospiraceae	o_Clostridiales	c_Clostridia	p_Firmicutes	d_Bacteria
//  BVAB2	g_Acetivibrio	f_Ruminococcaceae	o_Clostridiales	c_Clostridia	p_Firmicutes	d_Bacteria
//  BVAB3	g_Acetivibrio	f_Ruminococcaceae	o_Clostridiales	c_Clostridia	p_Firmicutes	d_Bacteria
//  Dialister_sp._type_1	g_Dialister	f_Veillonellaceae	o_Clostridiales	c_Clostridia	p_Firmicutes	d_Bacteria

//    - a reference fasta file of sequences representing species in the lineage file
//    - a species taxon file <seq ID> => <species name>
//    - an output directory

// The program builds

// - a model tree reflecting the species lineage data
//   parent/child structure with leaf lables being species names and internal
//   nodes corresponding to higher taxonomic ranks

// - For each node of the model tree (except root) a fasta file of all ref seq's
//    corresopnding to the node's subtree

// - A file of absolute paths to just created fasta files 3. A file of taxonomic
//    assignments of internal nodes

// Pawel Gajer
//   - Major rewrite: March 22, 2017
//   - Initial version: October 14, 2013

// Usage
//    buildModelTree -l <species lineage file> -i <fasta file> -t <species taxon file> -o <output dir>

// Example
//    cd ~/devel/packages/vaginal_species_oct18_2013

//    buildModelTree -l Firmicutes_group_6_V3V4_final.spLineage -i Firmicutes_group_6_V3V4_final.fa -t Firmicutes_group_6_V3V4_final.tx -o Firmicutes_group_6_V3V4_MC_models_dir


#define PATH_MAX 1024

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <set>
#include <algorithm>
#include <iterator>

#include "Newick2.hh"
#include "IOCUtilities.h"
#include "IOCppUtilities.hh"

using namespace std;

//----------------------------------------------------------- printUsage ----
void printUsage( const char *s )
{
    cout << endl

	 << "USAGE " << endl << endl

	 << "Create model tree and reference fasta files" << endl << endl
	 << s << "-l <lineage file> -i <fasta file> -o <output directory> [Options]" << endl << endl

	 << "\tOptions:\n"
	 << "\t-o <dir>         - output directory for MC taxonomy files\n"
	 << "\t-i <faFile>      - input fasta file that will be used to construct reference fasta files for MC models\n"
	 << "\t-l <lineageFile> - species lineage file\n"
	 << "\t-t <txFile>      - species taxon file\n"
	 << "\t-q|--quiet       - suppers pregress messages\n"
	 << "\t-v|--verbose     - verbose mode\n\n"
	 << "\t-h|--help        - this message\n\n"

	 << "\n\tExample\n"

	 << "\n\tcd Firmicutes_group_6_V3V4_dir" << endl

	 << "\t" << s << " -l Firmicutes_group_6_V3V4.spLineage -i Firmicutes_group_6_V3V4.fa -t Firmicutes_group_6_V3V4.tx -o MC_models_dir" << endl << endl

	 << endl;
}


//----------------------------------------------- printHelp ----
void printHelp( const char *s )
{
    cout << endl
	 << "Given fasta and lineage files, construct reference tree and fasta files for MC models\n\n";

    printUsage(s);
}

//================================================= inPar_t ====
//! holds input parameters
class inPar_t
{
public:
  inPar_t();
  ~inPar_t();

  char *outDir;             /// output directory for MC taxonomy files
  char *lineageFile;        /// species lineage file
  char *faFile;             /// fasta file
  char *txFile;             /// species taxon file
  bool verbose;
  bool quiet;

  void print();
};

//------------------------------------------------- constructor ----
inPar_t::inPar_t()
{
  outDir      = NULL;
  lineageFile = NULL;
  faFile      = NULL;
  txFile      = NULL;
  verbose     = false;
  quiet       = false;
}

//------------------------------------------------- constructor ----
inPar_t::~inPar_t()
{
  if ( outDir )
    free(outDir);

  if ( lineageFile )
    free(lineageFile);

  if ( faFile )
    free(faFile);

  if ( txFile )
    free(txFile);
}

//------------------------------------------------------- print ----
void inPar_t::print()
{
  if ( lineageFile )
    cerr << lineageFile << endl;
  else
    cerr << endl;

  cerr << "faFile=\t\t";
  if ( faFile )
    cerr << faFile << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "txFile=\t\t";
  if ( txFile )
    cerr << txFile << endl;
  else
    cerr << "MISSING" << endl;

  cerr << "outDir=\t\t";
  if ( outDir )
    cerr << outDir << endl;
  else
    cerr << "MISSING" << endl;
}


//============================== local sub-routines =========================
void parseArgs( int argc, char ** argv, inPar_t *p );


// ===================================================================
//                             main
// ===================================================================

int main(int argc, char **argv)
{
  //-- setting up init parameters
  inPar_t *inPar = new inPar_t();

  //-- parsing input parameters
  parseArgs(argc, argv, inPar);

  if ( inPar->verbose )
    inPar->print();

  if ( !inPar->lineageFile )
  {
    cerr << endl << "ERROR: Missing lineage file. Please specify it with the -l flag." << endl;
    printHelp(argv[0]);
    exit(1);
  }

  if ( !inPar->faFile )
  {
    cerr << endl << "ERROR: Missing reference seq's fasta file. Please specify it with the -i flag." << endl;
    printHelp(argv[0]);
    exit(1);
  }

  if ( !inPar->txFile )
  {
    cerr << endl << "ERROR: Missing species taxon file. Please specify it with the -t flag." << endl;
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
    cerr << endl << "ERROR: Missing Output directory. Please specify it with the -o flag." << endl;
    printHelp(argv[0]);
    exit(1);
  }


  char *lineageFile = inPar->lineageFile;
  char *faFile      = inPar->faFile;
  char *txFile      = inPar->txFile;
  char *outDir      = inPar->outDir;

  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Generating model tree from %s\n",lineageFile);
  NewickTree_t nt;
  nt.loadFullTxTree( lineageFile );

  //nt.printTree();

  string treeFile = string(outDir) + string("/model.tree");
  FILE *out = fOpen(treeFile.c_str(), "w");
  nt.writeTree(out);
  #if 0
  fprintf(stderr,"\n\tModel tree written to %s\n\n", treeFile.c_str());
  #endif

  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Loading ref tx file %s\n", txFile);
  char ***txTbl;
  int nRows, nCols;
  readCharTbl( txFile, &txTbl, &nRows, &nCols );

  // printf("nRows: %d\n", nRows);
  // printCharTbl(txTbl, 10, nCols); // test

  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Creating (tx => set of seqIDs) map\n");
  // typedef map<string, set<string> > strSet_t; // defined in IOCppUtilities.hh
  strSet_t tx2seqIDs;
  txTbl2txSet( txTbl, nRows, tx2seqIDs);

  #if 0
  strSet_t::iterator it;
  for ( it = tx2seqIDs.begin(); it != tx2seqIDs.end(); it++ )
  {
    cout << it->first << ": ";
    printStringSet(it->second);
  }
  #endif


  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Extending tx2seqIDs to internal nodes\n");
  nt.txSet2txTree( tx2seqIDs );

  #if 0
  printf("--- Printing key values of updated tx2seqIDs\n");
  strSet_t::iterator itr;
  for ( itr = tx2seqIDs.begin(); itr != tx2seqIDs.end(); itr++ )
    cout << itr->first << endl;
  cout << endl;
  #endif

  #if 0
  strSet_t::iterator it;
  for ( it = tx2seqIDs.begin(); it != tx2seqIDs.end(); it++ )
  {
    cout << it->first << ": ";
    printStringSet(it->second);
  }
  #endif

  #if 0
  // in order to test txSet2txTree() I am going to update tree labels so they
  // include the number of sequences associated with each node
  nt.updateLabels( tx2seqIDs );

  FILE *out1 = fOpen("testTree.tree","w");
  nt.writeTree(out1);
  fprintf(stderr,"--- Test tree written to testTree.tree\n");
  exit(1);
  #endif


  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Creating <interna node> => <taxonomy> table inodeTx\n");
  typedef map<string, string> str2str_t;
  map<string, string> inodeTx;
  nt.inodeTx( lineageFile, inodeTx );

  if ( !inPar->quiet )
    printf("--- Writing inodeTx to a file\n");
  string inodeTxFile = string(outDir) + string("/inode.tx");
  out = fOpen(inodeTxFile.c_str(), "w");
  map<string, string>::iterator it1;
  for ( it1 = inodeTx.begin(); it1 != inodeTx.end(); it1++ )
    fprintf(out, "%s\t%s\n", (it1->first).c_str(), (it1->second).c_str());
  fclose(out);

  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Loading ref fasta file %s\n",faFile);
  str2str_t seqTbl;
  readFasta( faFile, seqTbl);

  #if 0
  int i, nSeq = 10;
  str2str_t::iterator it2;
  for ( it2 = seqTbl.begin(), i = 0; it2 != seqTbl.end() && i < nSeq; it2++, i++ )
    cout << it2->first << "\t" << it2->second << endl;
  #endif

  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Create fasta file for each species\n");
  map<string, set<string> >::iterator it3;
  for ( it3 = tx2seqIDs.begin(); it3 != tx2seqIDs.end(); it3++ )
  {
    if ( !inPar->quiet )
      printf("\r--- processing %s", (it3->first).c_str());
    string outFile = string(outDir) + string("/") + it3->first + string(".fa");
    FILE *out = fOpen(outFile.c_str(), "w");

    set<string> seqIDs = it3->second;
    set<string>::iterator it4;
    for ( it4 = seqIDs.begin(); it4 != seqIDs.end(); it4++ )
      fprintf(out, ">%s\n%s\n", (*it4).c_str(), seqTbl[ *it4 ].c_str() );

    fclose(out);
  }

  // ==================================================================
  if ( !inPar->quiet )
    printf("\r--- Creating a file with absolute paths to just created fasta files\n");
  string pathsFile = string(outDir) + string("/spp_paths.txt");
  out = fOpen(pathsFile.c_str(), "w");
  for ( it3 = tx2seqIDs.begin(); it3 != tx2seqIDs.end(); it3++ )
  {
    string outFile = string(outDir) + string("/") + it3->first + string(".fa");
    char fullpath[PATH_MAX];
    realpath(outFile.c_str(), fullpath);
    fprintf(out, "%s\n", fullpath);
  }
  fclose(out);

  if ( !inPar->quiet )
    printf("\r--- Output written to %s\n", outDir);

  return EXIT_SUCCESS;
}

//----------------------------------------------------------- parseArgs ----
//! parse command line arguments
void parseArgs( int argc, char ** argv, inPar_t *p )
{
  int c, errflg = 0;
  optarg = NULL;

  static struct option longOptions[] = {
    {"fasta-file"   ,required_argument, 0, 'i'},
    {"taxon-file"   ,required_argument, 0, 't'},
    {"lineage-file" ,required_argument, 0, 'l'},
    {"out-dir"      ,required_argument, 0, 'o'},
    {"help"         ,no_argument,       0,   0},
    {"quiet"        ,no_argument,       0, 'q'},
    {"verbose"      ,no_argument,       0,   0},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv,"i:l:o:t:hvq",longOptions, NULL)) != -1)
    switch (c)
    {
      case 'o':
	p->outDir = strdup(optarg);
	break;

      case 'l':
	p->lineageFile = strdup(optarg);
	break;

      case 'i':
	p->faFile = strdup(optarg);
	break;

      case 't':
	p->txFile = strdup(optarg);
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
}
