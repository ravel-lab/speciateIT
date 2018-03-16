//
// tree_stitcher.cc
//
// Description

// Given a tree, T, and a table:

//    <leaf name> => <tree path>

// of trees labelled by the leaves of T, join the root of each tree in that
// table to the corresponding leaf of T.

// Syntax
//    tree_stitcher -i <tree file> -j <children tree file> -o <stitched tree file name> [Options]

// Example

//    tree_stitcher -i exMaster.tree -j exChnTrees.txt -t exStitched.tree


#define PATH_MAX 1000

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

#include "Newick.hh"
#include "IOCUtilities.h"
#include "CUtilities.h"
#include "IOCppUtilities.hh"

using namespace std;

//----------------------------------------------------------- printUsage ----
void printUsage()
{
  cout << endl
       << "USAGE " << endl << endl

       << "Create model tree and reference fasta files" << endl << endl
       << "\ttree_stitcher -i <tree file> -m <children trees mapping file> -o <stitched tree file name> [Options]" << endl << endl

       << "\tOptions:\n"
       << "\t--tree-file|-i <treeFile>   - master tree file\n"
       << "\t--mapping-file|-m <mapFile> - children trees file\n"
       << "\t--out-file|-o <outFile>     - stitched tree file\n"
       << "\t--quiet|-q                  - suppers pregress messages\n"
       << "\t--verbose|-v                - verbose mode\n"
       << "\t--debug                     - debug mode\n"
       << "\t--help|-h                   - this message\n\n"

       << "\n\tExample\n"

       << "\ttree_stitcher -i exMaster.tree -m exChnTrees.txt -o exStitched.tree" << endl << endl

       << endl;
}


//----------------------------------------------- printHelp ----
void printHelp()
{
  cout << "Given a tree, T, and a table:\n\n"
       << "   <leaf name> => <tree path>\n\n"
       << "of trees labelled by the leaves of T, the root of each tree in the\n"
       << " table is joined to the corresponding leaf of T\n\n";

  printUsage();
}

//================================================= inPar_t ====
//! holds input parameters
class inPar_t
{
public:
  inPar_t();
  ~inPar_t();

  char *treeFile;  /// master tree file
  char *mapFile;   /// children trees mapping file
  char *outFile;   /// stitched tree output file
  int quiet;
  int debug;
  int verbose;

  void print();
};

//------------------------------------------------- constructor ----
inPar_t::inPar_t()
{
  treeFile    = NULL;
  mapFile     = NULL;
  outFile     = NULL;
  verbose     = 0;
  quiet       = 0;
  debug       = 0;
}

//------------------------------------------------- constructor ----
inPar_t::~inPar_t()
{
  if ( treeFile )
    free(treeFile);

  if ( mapFile )
    free(mapFile);

  if ( outFile )
    free(outFile);
}

//------------------------------------------------------- print ----
void inPar_t::print()
{
  cerr << "treeFile=\t\t";
  if ( treeFile )
    cerr << treeFile << endl;

  cerr << "mapFile=\t\t";
  if ( mapFile )
    cerr << mapFile << endl;

  cerr << "outFile=\t\t";
  if ( outFile )
    cerr << outFile << endl;
}

//============================== local sub-routines =========================
void parseArgs( int argc, char ** argv, inPar_t *p );

void stitch_tree( NewickNode_t *node, // master tree node
                  map<string, NewickTree_t *> &trMap );

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

  if ( !inPar->treeFile )
  {
    errorMsg("Missing tree file");
    printHelp();
    exit(1);
  }

  if ( !inPar->mapFile )
  {
    cerr << endl << "ERROR: Missing children tree's file." << endl;
    printHelp();
    exit(1);
  }

  if ( !inPar->outFile )
  {
    cerr << endl << "ERROR: Missing output file name." << endl;
    printHelp();
    exit(1);
  }

  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Reading master tree\n");

  NewickTree_t *mTr = readNewickTree( inPar->treeFile );
  if ( !mTr )
  {
    fprintf(stderr,"\n\n\tERROR: Could not load master tree from %s\n\n", inPar->treeFile);
    exit(EXIT_FAILURE);
  }

  if ( 0 && inPar->debug )
    mTr->printTree();

  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Extracting leaves of the master tree\n");

  vector<string> mTrLeaves;
  mTr->leafLabels( mTr->root(), mTrLeaves);


  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Loading tree mapping file\n");

  map<string, string> trFileMap;
  read_2s_tbl( inPar->mapFile, trFileMap );

  //int nChTrees = trFileMap.size();

  if ( inPar->debug )
  {
    map<string, string>::iterator it;
    for ( it = trFileMap.begin(); it != trFileMap.end(); it++ )
    {
      fprintf(stderr,"|%s|\t|%s|\n", it->first.c_str(), it->second.c_str());
    }
    fprintf(stderr,"\n\n");
  }

  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Creating mapping: leafID => treeFile\n");

  map<string, NewickTree_t *> trMap;
  map<string, string>::iterator it;
  for ( it = trFileMap.begin(); it != trFileMap.end(); it++ )
  {
    // check if it->first is a leaf label of the master tree
    if ( ! exists_in_vector( mTrLeaves, it->first ) )
    {
      fprintf(stderr,"\n\n\tERROR: %s is not a leaf of the master tree\n\n", it->first.c_str() );
      exit(EXIT_FAILURE);
    }

    if ( inPar->debug )
      fprintf(stderr,"Processing %s with %s\n", it->first.c_str(), it->second.c_str() );

    NewickTree_t *tr = readNewickTree( it->second.c_str() );
	if ( !tr )
    {
      fprintf(stderr,"\n\n\tERROR: Could not load %s\n\n", it->second.c_str() );
      exit(EXIT_FAILURE);
    }

    if ( inPar->debug )
    {
      tr->printTree();
      fprintf( stderr,"Tree loaded\n" );
    }

    trMap[ it->first ] = tr;

    if ( inPar->debug )
      fprintf( stderr,"trMap assignment made\n\n" );
  }

  if ( inPar->debug )
  {
    fprintf( stderr,"Attempt to plot all trees\n\n" );
    map<string, NewickTree_t *>::iterator it2;
    for ( it2 = trMap.begin(); it2 != trMap.end(); it2++ )
    {
      cout << it2->first << ": ";
      (it2->second)->printTree();
    }
  }

  // ==================================================================
  if ( !inPar->quiet )
    printf("--- Traversing leaves of the master tree and stitching to them children trees\n");
  stitch_tree( mTr->root(), trMap );

  if ( inPar->debug )
    mTr->printTree();

  FILE *out = fOpen( inPar->outFile, "w" );
  mTr->writeTree( out );
  fclose(out);

  if ( !inPar->quiet )
    printf("\r--- Output written to %s\n", inPar->outFile);

  return EXIT_SUCCESS;
}

// ===================================================================
//                   subroutine definitions
//===================================================================

//----------------------------------------------- stitch_tree ----
  // post order traversal of a tree
void stitch_tree( NewickNode_t *node, // master tree node
                  map<string, NewickTree_t *> &trMap )
{
  int numChildren = node->children_m.size();
  if ( numChildren==0 ) // leaf
  {
    // check if node->label is a key of trMap
    if ( exists_in_map( trMap, node->label ) )
    {
      // stitch  tr = trMap[ node->label ] to node
      NewickTree_t * tr = trMap[ node->label ];
      node->stitch( tr->root() );
    }
  }
  else
  {
    for (int i = 0; i < numChildren; i++)
      stitch_tree( node->children_m[i], trMap );
  }
}

//----------------------------------------------------------- parseArgs ----
//! parse command line arguments
void parseArgs( int argc, char ** argv, inPar_t *p )
{
  int c, errflg = 0;
  optarg = NULL;

  static struct option longOptions[] = {
    {"tree-file"    ,required_argument, 0, 'i'},
    {"mapping-file" ,required_argument, 0, 'm'},
    {"out-file"     ,required_argument, 0, 'o'},
    {"help"         ,no_argument,       0, 'h'},
    {"quiet"        ,no_argument,       0, 'q'},
    {"debug"        ,no_argument, &p->debug, 'd'},
    {"verbose"      ,no_argument,       0, 'v'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv,"hi:m:o:qt:v",longOptions, NULL)) != -1)
    switch (c)
    {
      case 'd':
        p->debug = 1;
        break;

      case 'h':
        printHelp();
        exit (EXIT_SUCCESS);
        break;

      case 'i':
        p->treeFile = strdup(optarg);
        break;

      case 'm':
        p->mapFile = strdup(optarg);
        break;

      case 'o':
        p->outFile = strdup(optarg);
        break;

      case 'q':
        p->quiet = 1;
        break;

      case 'v':
        p->verbose = 1;
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
    printUsage();
    cerr << "Try '" << argv[0] << " -h' for more information" << endl;
    exit (EXIT_FAILURE);
  }
}
