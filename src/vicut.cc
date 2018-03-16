/*
  vicut source file

  March 26, 2016: changed the algorithm

  When internal cut node is found and there is more then one annotation in the
  corresponding cluster, then majoriy vote is used to modify the annotation, so
  that all leaves have the same one

  Nov 2, 2012: changing interface to

  vicut -t <tree file> -a <annotation file> -o <output dir> [Options]

  where <tree file> is in the newick format.

  I am dropping for now the -f options and associated options: -l, -d.

  One would have to create a tree outside of vicut. I am not going to use cluster.{c,h} anymore.


  NOTE: Taxonomy is determined in NewickTree_t::saveNAtaxonomy()
  See description of this routine for more details.

  Example:

  [~/devel/16SMC2/HMP]$ vicut -t gg12_8_seqs_nr99_align_105_2231_nr.filter_lbld_rooted4.tree -a ../GreenGenesDB_2012_8/Versioned_gg_12_8_phylum.taxonomy  -o gg12_8_seqs_nr99_align_105_2231_nr.filter_lbld_rooted4_phylum_vicut_Dir


  [~/devel/16SMC2/HMP]$ vicut -c 5 -x 10 -t gg12_8_seqs_nr99_align_105_2231_nr.filter_lbld_rooted4.tree -q seqs_nr99_nochimera_nr_105_2231.seqID -a ../GreenGenesDB_2012_8/Versioned_gg_12_8_genus2.taxonomy -o gg12_8_seqs_nr99_align_105_2231_nr.filter_lbld_rooted4_genus_vicut_Dir

*/

#include <getopt.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <queue>

#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>
//#include <tr1/unordered_map>

using namespace std;

#include "Newick.hh"
#include "IOUtilities.h"
#include "IOCppUtilities.hh"


//---------------------------------- printUsage ------------------
void printUsage( const char *s )
{
    cout << endl
         << "USAGE: " << s << " -t <tree file> -a <annotation file> -q <query seq ID file> -o <output dir> [Options]\n\n"
         << endl
         << "\tOptions:\n"
         << "\t-t,--tree-file <tree file>          - tree file in the newick format\n"
         << "\t-a,--ann-file <annotation file>     - annotation file\n"
         << "\t-u,--tax-table-file <tax tbl file>  - tax table file with two columns <tax_name>\t<parent_name> \n"
         << "\t-q,--query-file <query seq ID file> - file containing sequence IDs of query sequences. Used only for cluster stats.\n"
         << "\t-o,--out-dir <output dir>           - output directory will include minNodeCut.idx, minNodeCut.cltrs and minNodeCut.cltrs2 files\n"
         << "\t-x,--min-taxon <minT>               - min number of reference sequences of a given taxon in a given cluster for the taxon to be reported in minNodeCut_NAge<minT>.cltrsStats file; default value = 1\n"
         << "\t-c,--min-query <minQ>               - min number of query sequences in a cluaster for the cluster to be reported in minNodeCut_NAge<minQ>_TXge<minT>_taxons.txt file; default value = 1\n"
         << "\t--quiet                             - suppress progress messages\n"
         << "\t-p,--print-min-node-cut-cltrs       - print min mode cut clusters\n\n"

         << "If run without -q option, the output directory contains files\n"
         << "\tminNodeCut.cltrs\n"
         << "\tminNodeCut.cltrs2\n"
         << "\tminNodeCut.cltrsStats\n\n"
         << "With -q option, the following files are also generated\n"
         << "\tminNodeCut_NAge1.cltrsStats\n"
         << "\tminNodeCut_NAge1_TXge1_querySeqs.taxonomy\n"
         << "\tminNodeCut_NAge1_TXge1_refSeqIDs.txt\n"
         << "\tminNodeCut_NAge1_TXge1_taxons.txt\n"
         << "\tminNodeCut_NAl1.cltrsStats\n\n"
         << "Description of the output files is in the README file of the output directory.\n\n"

         << "\n\tExamples: \n"
         << "\tcd examples\n"
         << "\t" << s << " -x 1 -c 1 -t fig1v2.nw -a fig1v2_ann.txt -o fig1v2Dir -v\n"
         << "\t" << s << " -x 1 -c 1 -t fig1v2.nw -a fig1v2_ann.txt -q fig1v2_query.ids -o fig1v2Dir -v\n"
         << endl;
}

//--------------------------------- printHelp -----------------
void printHelp( const char *s )
{
    cout << endl
         << "vicut performs vi-cut clustering of a non-necessary binary tree with anannotation data.\n"
         << "The input must consists of a tree file (newick format) with an annotation file."
         << "For more info on vi-cut clustering see http://www.cbcb.umd.edu/VICut\n\n";

    printUsage(s);
}

// generate README file with short description of the output
void generateREADMEfile(const char *outDir)
{
    int n1 = strlen(outDir);
    char readmeFileName[] = "/README";
    int n = strlen(readmeFileName);
    char *routFile = (char*)malloc((n1+n+1)*sizeof(char));
    strcpy(routFile,outDir);
    strcat(routFile,readmeFileName);
    //writeREADME(routFile);
    FILE *fh = fOpen(routFile,"w");
    const char str[] = "The output of vicut consists of the following files.\n\n"

      "minNodeCut.cltrs: file listing reference and query sequence IDs, their cluster\n"
      "assignment and annotation (only for reference sequences). Annotation field of\n"
      "query sequences is 'NA'. Here is a fragment of such a file\n\n"

      "readId	clstr	annot\n"
      "993695	312	Actinobacteria\n"
      "691953	312	Actinobacteria\n"
      "UAB130.10.5.06142010_138096	312	NA\n"
      "UAB017.10.6.12162009_361558	312	NA\n"
      "UAB060.8.5.02262010_500708	312	NA\n"
      "UAB130.2.1.04152010_127238	312	NA\n"
      "143824	148	Synergistetes\n"
      "200701	148	Synergistetes\n"
      "UAB135.5.6.05182010_306482	148	NA\n"
      "495024	148	Synergistetes\n"
      "...\n"
      "...\n\n"

      "minNodeCut.cltrs2: file listing clustes and reads within each cluster with their\n"
      "annotation (for reference sequences). For example\n\n"

      "Cluster 109:\n"
      "\t256248	Proteobacteria\n"
      "Cluster 110:\n"
      "\t200370	OP11\n"
      "\t254839	OP11\n"
      "\t201943	OP11\n"
      "Cluster 111:\n"
      "\t593456	OP11\n"
      "\t255767	OP11\n"
      "\t154380	OP11\n"
      "\t...\n"
      "\t...\n\n"

      "minNodeCut.cltrsStats: file lists for each cluster annotations of all sequences in\n"
      "this cluster together with number of reads and percentage of the reads in the cluster\n"

      "minNodeCut.NAannStats: file listing annotation strings (taxons) present in the\n"
      "vicut clusters containg query sequences together with the maximum percentage\n"
      "share of this taxon and the corresponding sequence count\n\n"

      //"minNodeCut.NAcltrs: "

      "minNodeCut.NAcltrsStats: file listing clusters containing query sequences and taxons\n"
      "present in them with the number of sequences per taxon and the percentage share\n"
      "of the taxon in the whole cluster. For example: \n\n"
      "Cluster 148:\n"
      "\tSynergistetes	769	98.2%\n"
      "\tNA	5	0.6%\n"
      "\tFirmicutes	4	0.5%\n"
      "\tArmatimonadetes	2	0.3%\n"
      "\tNKB19	1	0.1%\n"
      "\tTenericutes	1	0.1%\n"
      "\tThermotogae	1	0.1%\n"
      "\tTOTAL	783	100%\n"
      "Cluster 160:\n"
      "\tFirmicutes	123831	82.0%\n"
      "\tNA	18396	12.2%\n"
      "\tCyanobacteria	5348	3.5%\n"
      "\tFusobacteria	1957	1.3%\n"
      "\tTenericutes	1214	0.8%\n"
      "\t...\n"
      "\t...\n\n";

    fprintf(fh, "%s", str);
    fclose(fh);
}

double vicut(vector<int> &nodeCut, map<int, NewickNode_t *> &idx2node, int *annIdx, int nIds,
             map<string,int> &uAnn, map<string, int> &leafIdx, bool verbose);

// ------------------------------------------------------------
//                          main
// ------------------------------------------------------------
int main(int argc, char **argv)
{
    bool verbose      = 0;
    bool quiet        = 0;
    bool printMinNodeCutCltrs = 0;
    char *outDir      = 0;
    char *annFile     = 0;
    char *taxTblFile  = 0;
    char *queryFile   = 0;
    char *treeFile    = 0;
    optarg            = 0;
    int minNA         = 1; // min number of query elements in a cluster
    int minAnn        = 1; // min number of reference sequnces of a given taxon for a taxon to be selected
    int c, n, errflg  = 0;

    static struct option longOptions[] = {
      {"verbose"          ,no_argument,       0, 'v'},
      {"ann-file"         ,required_argument, 0, 'a'},
      {"tax-table-file"   ,required_argument, 0, 'u'},
      {"query-file"       ,required_argument, 0, 'q'},
      {"tree-file"        ,required_argument, 0, 't'},
      {"out-dir"          ,required_argument, 0, 'o'},
      {"min-query"        ,required_argument, 0, 'c'},
      {"min-taxon"        ,required_argument, 0, 'x'},
      {"print-min-node-cut-cltrs", no_argument, 0, 'p'},
      {"quiet",            no_argument,       0, 'e'},
      {"help"             ,no_argument,       0, 'h'},
      {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv,"a:c:q:o:t:x:u:vpqh",longOptions, NULL)) != -1)
      switch (c)
      {
        case 'a':
          n = strlen(optarg) + 1;
          annFile = (char*)malloc(n*sizeof(char));
          strncpy(annFile, optarg, n);
          break;

        case 'u':
          n = strlen(optarg) + 1;
          taxTblFile = (char*)malloc(n*sizeof(char));
          strncpy(taxTblFile, optarg, n);
          break;

        case 'c':
          minNA = atoi(optarg);
          if ( minNA == 0 )
          {
            fprintf(stderr,"ERROR: The argument, %s, of the '-c' flag cannot be converted to integer.\n",optarg);
            exit( EXIT_FAILURE);
          }
          break;

        case 'e':
          quiet = true;
          break;

        case 'x':
          minAnn = atoi(optarg);
          if ( minAnn == 0 )
          {
            fprintf(stderr,"ERROR: The argument, %s, of the '-x' flag cannot be converted to integer.\n",optarg);
            exit( EXIT_FAILURE);
          }
          break;

        case 'q':
          n = strlen(optarg) + 1;
          queryFile = (char*)malloc(n*sizeof(char));
          strncpy(queryFile, optarg, n);
          break;

        case 't':
          n = strlen(optarg) + 1;
          treeFile = (char*)malloc(n*sizeof(char));
          strncpy(treeFile, optarg, n);
          break;

        case 'o':
          n = strlen(optarg) + 1;
          outDir = (char*)malloc(n*sizeof(char));
          strncpy(outDir, optarg, n);
          break;

        case 'v':
          verbose = true;
          break;

        case 'p':
          printMinNodeCutCltrs = true;
          break;

        case 'h':
        case 0:
          printHelp(argv[0]);
          exit (EXIT_SUCCESS);
          break;

        default:
          cerr << "\n========================================\n"
               << "ERROR: Unrecognized option " << (char)c << "\n"
               << "==========================================\n" << endl;
          ++errflg;
          break;
      }

    if ( errflg )
    {
      printUsage(argv[0]);
      exit (EXIT_FAILURE);
    }

    if ( !annFile )
    {
      fprintf(stderr,"ERROR: annotation file missing!\n");
      printUsage(argv[0]);
      exit (EXIT_FAILURE);
    }

    // if ( !queryFile )
    // {
    //   fprintf(stderr,"ERROR: query file missing!\n");
    //   printUsage(argv[0]);
    //   exit (EXIT_FAILURE);
    // }

    if ( !treeFile )
    {
      fprintf(stderr,"ERROR: tree file missing!\n");
      printUsage(argv[0]);
      exit (EXIT_FAILURE);
    }

    if ( !outDir )
    {
      fprintf(stderr,"ERROR: output directory missing!\n");
      printUsage(argv[0]);
      exit (EXIT_FAILURE);
    }

    clock_t tic = clock();

    // reading tree using readTree() defined in Newick.c
    if ( !quiet )
      fprintf(stderr, "-- Reading tree data from %s\n", treeFile);

    NewickTree_t *nt = readNewickTree( treeFile );
    if ( !nt )
    {
      fprintf(stderr,"Could not load Newick tree from %s\n", treeFile);
      exit(EXIT_FAILURE);
    }

    if ( verbose )
    {
      // fprintf(stderr,"leafLables\n");
      // for ( int i = 0; i < nLeaves; i++ )
      //   fprintf(stderr, "i: %d, %s\n",i, leafLabel[i]);
      // fprintf(stderr,"\n\n");

      nt->printTree(1);
    }

    // hash table: <node index> => <node pointer>
    map<int, NewickNode_t *> idx2node;
    nt->indexNewickNodes(idx2node);

    int nLeaves =  nt->getNleaves();
    char **leafLabel = nt->leafLabels();

    map<string, int> leafIdx; // maps leaf label to the corresponding index of leafLabel
    for ( int j = 0; j < nLeaves; j++ )
      leafIdx[string(leafLabel[j])] = j;

    // parsing annotation data
    if ( !quiet )
      fprintf(stderr, "-- Parsing annotation data from %s\n", annFile);
    map<string,int> annToIdx; // table of unique annotation strings and their indices
    int *annIdx = parseAnnotation(annFile, leafLabel, nLeaves, annToIdx, verbose);

    // fprintf(stderr,"annToIdx.size()=%d\n",(int)annToIdx.size());
    // map<string,int>::iterator itr2;
    // for (itr2 = annToIdx.begin(); itr2 != annToIdx.end(); itr2++)
    //   fprintf(stderr,"%s\t%d\n",(*itr2).first.c_str(),(*itr2).second);


    // parsing query sequence IDs
    if ( queryFile )
    {
      if ( !quiet )
        fprintf(stderr, "-- Parsing query sequence IDs from %s\n", queryFile);
      map<string, bool> queryIDs; // table of query sequence IDs; it assigns to each query ID 'true'. The purpose of the table is to quickly find out if a sequence is a query sequence.
      parseQueryIDs(queryFile, queryIDs);

      if ( verbose )
      {
        printf("queryIDs (%d):\n", (int)queryIDs.size());
        map<string,bool>::iterator itr;
        int count = 0;
        for (itr = queryIDs.begin(); itr != queryIDs.end(); itr++)
        {
          printf("%s\n",(*itr).first.c_str());

          if ( count == 10 )
            break;
          count++;
        }
      }

      // update annIdx and annToIdx, so that annToIdx has indices for 'NA' and
      // 'Unclassified' strings with 'NA' corresponding to query sequences and
      // 'Unclassified' to reference sequences without annotation.

      // annToIdx['Unclassified'] = -1
      // annToIdx['NA'] = -2

      // Thus, annIdx on query sequences is -2 and on Unclassified sequences is -1.
      queryIDsToAnnIdx ( queryIDs, leafLabel, nLeaves, annToIdx, annIdx);
    }
    else
    {
      annToIdx[string("NA")] = -1;
    }

    map<int, string> idxToAnn;
    map<string,int>::iterator itr;
    for (itr = annToIdx.begin(); itr != annToIdx.end(); itr++)
      idxToAnn[(*itr).second] = (*itr).first;

    if ( verbose )
    {
      // print idxToAnn array
      //int count = nLeaves;
      int count = 1000;
      count = (count < nLeaves ? count : nLeaves);
      printf("--- annIdx and idxToAnn - first %d elements\n",count);
      printf("i\tlabel\tannIdx\tidxToAnn\n");
      for ( int i = 0; i < count; ++i )
        printf("%d\t%s\t%d\t%s\n",i,(idx2node[i]->label).c_str(),annIdx[i],idxToAnn[annIdx[i]].c_str());
      //if ( annIdx[i] > -1 )
      //printf("%d\t%d\t%s\n",i,annIdx[i],idxToAnn[annIdx[i]].c_str());
      // else
      // 	printf("%d\t%d\t-\n",i,annIdx[i]);
      printf("\n");
    }

    // vi-cut
    if ( !quiet )
      fprintf(stderr, "-- Running vicut\n");
    vector<int> nodeCut;
    double viDist = vicut(nodeCut, idx2node, annIdx, nLeaves, annToIdx, leafIdx, verbose); // modeType mode


    // print min-node-cut clusters
    if ( printMinNodeCutCltrs )
    {
      // print nodeCut;
      printf("min-node-cut:\n");
      int m = nodeCut.size();
      printf("min-node-cut clusters:\n");
      for ( int i = 0; i < m; ++i )
      {
        printf("Cluster %d:",i);
        //printLeaves( tree, nodeCut[i], leafLabel);
      }
    }

    if ( !quiet )
      fprintf(stderr, "-- Writing output files\n");

    char s[1024];
    sprintf(s,"mkdir -p %s",outDir);
    system(s);

    int n1 = strlen(outDir);
    char fileName[] = "/minNodeCut.cltrs";
    int n2 = strlen(fileName);
    char *outFile = (char*)malloc((n1+n2+1)*sizeof(char));
    strcpy(outFile,outDir);
    strcat(outFile,fileName);
    nt->saveCltrMemb( outFile, nodeCut, annIdx, idxToAnn);

    // save cluster data in another format
    n = strlen(outFile);
    char *outFile2 = (char*)malloc((n+2)*sizeof(char));
    strcpy(outFile2,outFile);
    outFile2[n] = '2';
    outFile2[n+1] = '\0';
    nt->saveCltrMemb2( outFile2, nodeCut, annIdx, idxToAnn);

    // save cluster stats (annotation counts)
    char fileName3[] = "/minNodeCut.cltrsStats";
    int n3 = strlen(fileName3);
    char *outFile3 = (char*)malloc((n1+n3+1)*sizeof(char));
    strcpy(outFile3,outDir);
    strcat(outFile3,fileName3);
    nt->saveCltrAnnStats( outFile3, nodeCut, annIdx, idxToAnn);

    int nNAs = 0;
    int nNAs_with_tx = 0;
    int tx_changed = 0;
    int nClades_modified = 0;

    if ( queryFile )
    {
      // save cluster stats only for clusters containg query sequences
      vector<string> selTx; // selected taxons (present in clusters with query sequences)

      if ( minNA == 0 )
      {
        char fileName4[] = "/minNodeCut.NAcltrsStats";
        int n4 = strlen(fileName4);
        char *outFile4 = (char*)malloc((n1+n4+1)*sizeof(char));
        strcpy(outFile4,outDir);
        strcat(outFile4,fileName4);
        selTx = nt->saveNAcltrAnnStats(outFile4, nodeCut, annIdx, idxToAnn);
      }
      else
      {
        char fileName4a[256]; // minNodeCut.NAcltrsStats file with clusters containing at least minNA query elements
        char fileName4b[256]; // minNodeCut.NAcltrsStats file with clusters less than minNA query elements
        char fileName4c[256]; // file containing sequence IDs from clusters with at least minNA query elements and taxons with at least minAnn elements.
        char fileName4d[256]; // file containing taxons with at least minAnn elements from clusters with at least minNA query elements
        char fileName4e[256]; // taxonomy assignments of query sequences from clusters containing either only query sequences
        // (then its just Unclassified, later OTU Id) or query sequences with at least minNA annotation sequences
        // and minAnn reference sequences; taxonomy is based on the most abundant reference taxon
        char fileName4f[256]; // log file for Unassigned clusters
        char fileName4g[256]; // log file of genus OTUs

        int n4a = sprintf(fileName4a, "/minNodeCut_NAge%d.cltrsStats",minNA);
        int n4b = sprintf(fileName4b, "/minNodeCut_NAl%d.cltrsStats",minNA);
        int n4c = sprintf(fileName4c, "/minNodeCut_NAge%d_TXge%d_refSeqIDs.txt",minNA, minAnn);
        int n4d = sprintf(fileName4d, "/minNodeCut_NAge%d_TXge%d_taxons.txt",minNA, minAnn);
        int n4e = sprintf(fileName4e, "/minNodeCut_NAge%d_TXge%d_querySeqs.taxonomy",minNA, minAnn);
        int n4f = sprintf(fileName4f, "/minNodeCut_NAge%d_TXge%d_unassigned.log",minNA, minAnn);
        int n4g = sprintf(fileName4g, "/minNodeCut_NAge%d_TXge%d_genus_OTUs.txt",minNA, minAnn);

        char *outFile4a = (char*)malloc((n1+n4a+1)*sizeof(char));
        strcpy(outFile4a,outDir);
        strcat(outFile4a,fileName4a);

        char *outFile4b = (char*)malloc((n1+n4b+1)*sizeof(char));
        strcpy(outFile4b,outDir);
        strcat(outFile4b,fileName4b);

        char *outFile4c = (char*)malloc((n1+n4c+1)*sizeof(char));
        strcpy(outFile4c,outDir);
        strcat(outFile4c,fileName4c);

        char *outFile4d = (char*)malloc((n1+n4d+1)*sizeof(char));
        strcpy(outFile4d,outDir);
        strcat(outFile4d,fileName4d);

        nt->saveNAcltrAnnStats( outFile4a, outFile4b, outFile4c, outFile4d,
                                nodeCut, annIdx, idxToAnn, minNA, minAnn );

        char *outFile4e = (char*)malloc((n1+n4e+1)*sizeof(char));
        strcpy(outFile4e,outDir);
        strcat(outFile4e,fileName4e);
        string logFile4e = string(outDir) + string("/minNodeCut_sizes_of_2_dominant_taxons.txt");

        nt->saveNAtaxonomy0( outFile4e, logFile4e.c_str(), nodeCut, annIdx, idxToAnn, minAnn,
                             nNAs, nNAs_with_tx, tx_changed, nClades_modified );

        if ( taxTblFile )
        {
          char *outFile4f = (char*)malloc((n1+n4f+1)*sizeof(char));
          strcpy(outFile4f,outDir);
          strcat(outFile4f,fileName4f);

          char *outFile4g = (char*)malloc((n1+n4g+1)*sizeof(char));
          strcpy(outFile4g,outDir);
          strcat(outFile4g,fileName4g);

          map<string, string> txTbl;
          read2cols(taxTblFile, txTbl);
          //nt->saveNAtaxonomy( outFile4e, outFile4f, outFile4g, nodeCut, annIdx, idxToAnn, minNA, txTbl);

          char dirName[256];
          int ndir = sprintf(dirName, "/genusTrees/");
          char *genusDir = (char*)malloc((n1+ndir+1)*sizeof(char));
          strcpy(genusDir,outDir);
          strcat(genusDir,dirName);

          //nt->printGenusTrees( genusDir, nodeCut, annIdx, minNA, idxToAnn, txTbl );

          //
          char fileName[256];
          int n = sprintf(fileName, "/minNodeCut_NAge%d_TXge%d.taxonomy",minNA, minAnn);
          char *outFile = (char*)malloc((n1+n+1)*sizeof(char));
          strcpy(outFile,outDir);
          strcat(outFile,fileName);

          vector<sppProf_t*> sppProfs;
          nt->getSppProfs( sppProfs, nodeCut, annIdx, minNA, idxToAnn, txTbl);

          vector<sppProf_t*> ancRootProfs;
          nt->getAncRoots( sppProfs, ancRootProfs);

          nt->printTx(outFile, ancRootProfs);
        }
      }


#if 0
      // write selTx to a file
      char fileName5[] = "/minNodeCut.NAcltrs";
      int n5 = strlen(fileName5);
      char *outFile5 = (char*)malloc((n1+n5+1)*sizeof(char));
      strcpy(outFile5,outDir);
      strcat(outFile5,fileName5);
      writeStrVector(outFile5, selTx);

      // write to a file annotation strings present in the clusters containg query sequences
      // together with the maximum percentage share of this taxon and the count of sequences
      char fileName6[] = "/minNodeCut.NAannStats";
      int n6 = strlen(fileName6);
      char *outFile6 = (char*)malloc((n1+n6+1)*sizeof(char));
      strcpy(outFile6,outDir);
      strcat(outFile6,fileName6);
      writeNAannStats(outFile6, tree, nodeCut, annIdx, idxToAnn);
#endif

      // generate README file with short description of the output
      generateREADMEfile(outDir);

      if ( queryFile  && !quiet )
      {
        printf("\nNumber of vicut clusters with taxonomy modified      = %d\n"
               "Number of query IDs                                  = %d\n"
               "Number of query IDs to which vicut assigned taxonomy = %d\n"
               "Number of non-query IDs which changed taxonomy       = %d\n\n",
               nClades_modified, nNAs, nNAs_with_tx, tx_changed);
      }


      if ( !quiet )
      {
        printf("\nNumber of clusters          = %d\n"
               "Number of annotation labels = %d\n"
               "--------------------------------------\n"
               "Difference                  = %d\n\n"
               "VI-distance                 = %.6g\n",
             (int)nodeCut.size(),
             (int)annToIdx.size() - 2,  // subtracting "NA' and "Unclassified"
             (int)(nodeCut.size() - (annToIdx.size() - 2)),
             viDist);

        printf("\nOutput written to\n%s\n",outDir);
      }

      // fprintf(stderr,"annToIdx.size()=%d\n",(int)annToIdx.size());
      // map<string,int>::iterator itr2;
      // for (itr2 = annToIdx.begin(); itr2 != annToIdx.end(); itr2++)
      //   fprintf(stderr,"%s\t%d\n",(*itr2).first.c_str(),(*itr2).second);

      free(outFile);
      free(outFile2);
    }

    // cleanup
    for ( int i = 0; i < nLeaves; ++i )
      free(leafLabel[i]);
    free(leafLabel);

    //free(tree);
    free(annIdx);

    if (annFile)
      free(annFile);

    if (treeFile)
      free(treeFile);

    clock_t toc = clock();
    double runTime =(double)(toc - tic) / CLOCKS_PER_SEC;

    if ( runTime > 60 )
    {
      int timeMin = (int)round(runTime / 60);
      int timeSec = (int)runTime % 60;
      if ( !quiet )
        printf("Completed in %dmin %dsec\n",  timeMin, timeSec);
    }
    else
    {
      if ( !quiet )
        printf("Completed in %f seconds\n",  runTime);
    }

    return 0;
}



// ---------------------------- resetMinCut --------------------------------
void resetMinCut(NewickNode_t *selNode, map<int, bool> &minCut)
/*
  traverse tree from node tree[i] and set minCut at all internal nodes of
  this sub-tree to 0
*/
{
    queue<NewickNode_t *> bfs;
    bfs.push(selNode);
    NewickNode_t *node;

    while ( !bfs.empty() )
    {
      node = bfs.front();
      bfs.pop();

      if ( node != selNode )
        minCut[node->idx] = 0;

      int numChildren = node->children_m.size();
      if ( numChildren != 0 )
      {
        for (int i = 0; i < numChildren; i++)
        {
          bfs.push(node->children_m[i]);
        }
      }
    }
}

// ------------------------------ H ----------------------------------------
inline double H(double p)
{
    double x = 0;
    if ( p != 0 )
      x = -p*log2(p);

    return x;
}

// ---------------------------- modeType -----------------------------------
enum modeType {
  low=0, // low min-node-cut
  high   // high min-node-cut
};


// ---------------------------- rankArray ----------------------------------
// rank elements of an int table
// lifted from http://www.cplusplus.com/forum/beginner/34046/
int* rankArray(int table[],int n)
{
    //int result[n];
    int *result = (int*)calloc(n, sizeof(int));

    for(int i=0;i<n;i++)
    {
      int rnk=0;
      for(int z=0;z<n;z++)
      {
        if(table[z]<table[i])
          rnk++;
      }
      result[i]=rnk;
    }
    return result;
}


// ------------------------------ vicut -------------------------------------
double vicut(vector<int> &nodeCut, map<int, NewickNode_t *> &idx2node, int *annIdx, int nLeaves,
             map<string,int> &uAnn, map<string, int> &leafIdx, bool verbose)
/*
  Given a not necessary binary tree 'nt' an annotation index array 'annIdx' and
  the number of tree leaves 'nLeaves' (equal also to the number of elements in
  'annIdx'; annIdx of unannotated elements is < 0), this routine computes a
  node-cut that minimizes variation of information distance between the
  clustering induced by annIdx and clusterings of the leaves of 'nt' induced by
  the node-cut (restricted to nodes derived from the annotated nodes). This
  minimum of variation of information node-cut is called min-mode-cut

  Parameters:

  nodeCut  - vector of tree node indices that constitute a min-node-cut

  Recall, that leaves have indices 0, 1, ... , (nLeaves-1) and internal nodes
  have indices: -1, -2, ... , -(nLeaves-1).

  idx2node - hash table mapping node index to a NewickNote_t pointer of that node

  annIdx   - array of annotation indices assigned to leaves of the tree; unannotated elements have annIdx of -1 or -2

  nLeaves  - number of leavs of the tree (equal to the number of elements of annIdx)

  uAnn     - map: <annotation string> => <int index>

  leafIdx  - maps leaf labels to the corresponding indices defined in leafLabel array (see: char **leafLabel = nt.leafLabels() )

  verbose  - if 1, print diagnostic messages

  Output: VI distance between annotation and min-node-cut clustering

*/
{
    // compute q(x) = p(x)*log2(p(x)) - 2* \sum_{d \in D} p(x,d)*log2(p(x,d))
    // where x is a node of tree, D is the set of unique annotation labels,
    // p(x) is the probability that an element with a known annotation would
    // fall into the cluster induced by x and p(x,d) is the probability
    // that a random annotated element falls into cluster x and has annotation d.
    // The formulas for p(x) and p(x,d) are
    // p(x) = |L(x)|/n
    // p(x,d) = |L(x)\cap A(d)|/n
    // where n is the number of elements (leavs of tree) that have known annotation
    // L(x) is the set of leaves in the subtree rooted at x, that have known annotation
    // and A(d) is the set of leaves (in the whole tree) that are known to have annotation d.

    BOOST_STATIC_ASSERT(true) __attribute__((unused));

#define VICUT_DEBUG 0

    // computing nAnnLeaves and A(d)
    int nAnnLeaves = 0;
    int nuAnn = uAnn.size(); // number of unique annotation strings
    //printf("nuAnn=%d\n",nuAnn);

    int *Acount = (int*)calloc(nuAnn, sizeof(int)); // Acount[i] = A(d) for i-th annotation label d
    for ( int i = 0; i < nLeaves; ++i )
    {
      if ( annIdx[i] > -1 )
      {
        Acount[annIdx[i]]++;
        nAnnLeaves++;
      }
    }

    // the entropy of the partition induced by annotation labels
    double HA = 0;
    for ( int i = 0; i < nuAnn; ++i )
      HA += H( Acount[i]/(double)nAnnLeaves );

#if VICUT_DEBUG
    fprintf(stderr,"\n\tNumber of leaves: %d\n\tNumber of annotated leaves: %d\n\tNumber of unique labels: %d\n\n", nLeaves, nAnnLeaves, nuAnn);
#endif

    // first lets compute q(x) at the leaves, which in the tree are represented
    // by integers 0, 1, ... , (nLeaves-1)
    // later we are going to compute q for internal nodes that are indexed by
    // -1, -2, ... , -(nLeaves-1)

    map<int, double> qMap;    // q(x) = 2*H(p(x,d)) - H(p(x))
    map<int, int*> AcountMap; // AcountMap[i][j] = number of leaves of the i-th node annotated
    // using the j-th annotation label
    map<int, int> countMap;   // countMap[i] = number of annotated leaves of the i-th node
    // countMap[i] = \sum_{j} AcountMap[i][j]
    double p = 1.0/(double)nAnnLeaves; // p(x)
    double Hp = H(p);         // entropy of p(x)
    //double jp = p;            // join probability p(x,d)

    // min-node-cut memebership map
    map<int, bool> minCut; // minCut[i] = 1 if the i-th node is a member of min-node-cut; 0 otherwise

    // initialize all leave nodes to 1 and internal nodes to 0
    for ( int i = 0; i < nLeaves; ++i )
      minCut[i] = 1;

    int minIdx = 0;
    pair<int, NewickNode_t *> pr;
    BOOST_FOREACH(pr, idx2node)
    {
      if ( pr.first < minIdx )
        minIdx = pr.first;
    }

#if VICUT_DEBUG
    fprintf(stderr,"\tminIdx: %d\n", minIdx);
    fprintf(stderr,"\tH(p): %.3f\n", Hp);
#endif

    for( int i = minIdx; i < 0; i++ )
      minCut[i] = 0;

    // min-node-cut internal node index
    //vector<vector<int> > minCutIdx; // minCutIdx[i] vector of indices of elements of i-th internal node subtree that are members of min-node-cut
    // in high mode we don't have to keep it
    // as a node is promoted to min-node-cut when
    // q(x) <= q(x.left) + q(x.right)
    // and in high mode it means demoting leaves of x and promoting x to min-node-cut

    // populating countMap, AcountMap and qMap for leaves
    for ( int i = 0; i < nLeaves; ++i )
    {
      if ( annIdx[i] > -1 )
      {
        countMap[i] = 1;

        AcountMap[i] = (int*)calloc(nuAnn, sizeof(int));
        AcountMap[i][annIdx[i]] = 1;

        qMap[i] = Hp; // = 2*H(jp) - Hp; // which is Hp as H(jp) = H(p)
      }
      else
      {
        countMap[i]  = 0;
        AcountMap[i] = NULL;
        qMap[i] = 0;
      }
    }

    // initializing countMap of internal nodes to -1
    for( int i = minIdx; i < 0; i++ )
      countMap[i] = -1;

    // computing qMap of internal nodes
    NewickNode_t *node;

    //for( int i = -1; i >= minIdx; i-- )
    int round = 0;
    int i = -1;
    while ( i >= minIdx )
    {
      node = idx2node[i];

      if ( node==NULL )
      {
        fprintf(stderr,"ERROR: i=%d node=NULL\n",i);
        exit(EXIT_FAILURE);
      }

      int nChildren = node->children_m.size();

      NewickTree_t nt;
      vector<string> leaves;
      nt.leafLabels(node, leaves);

#if VICUT_DEBUG
      fprintf(stderr,"\n-----------------------------\nRound %d  Internal Node %d\nChildren: ",round, i);
      for (int j = 0; j < (nChildren-1); j++)
        if ( (node->children_m[j])->idx < 0 )
          fprintf(stderr,"%d, ", (node->children_m[j])->idx);
        else
          fprintf(stderr,"%d (%s), ", (node->children_m[j])->idx, (node->children_m[j])->label.c_str());

      if ( (node->children_m[nChildren-1])->idx < 0 )
        fprintf(stderr,"%d\n", (node->children_m[nChildren-1])->idx);
      else
        fprintf(stderr,"%d (%s)\n", (node->children_m[nChildren-1])->idx, (node->children_m[nChildren-1])->label.c_str());

      // print leaf labels of the subtree of the given internal node
      fprintf(stderr,"Leaves: ");
      printVector(leaves, " ");
      //fprintf(stderr,"\n");

      fprintf(stderr,"-----------------------------\n");
#endif

      //#if VICUT_DEBUG
#if 0
      fprintf(stderr,"nChildren=%d\n",nChildren);
      for (int j = 0; j < nChildren; j++)
      {
        if ( (node->children_m[j])->idx >= 0 )
          fprintf(stderr,"nChildren(child%d(%d,%s))=%d\n",j, (node->children_m[j])->idx, (node->children_m[j])->label.c_str(), (int)(node->children_m[j])->children_m.size());
        else
          fprintf(stderr,"nChildren(child%d(%d))=%d\n",j, (node->children_m[j])->idx, (int)(node->children_m[j])->children_m.size());
      }
#endif

      bool hasCounts = false;
      for (int j = 0; j < nChildren; j++)
        hasCounts |= countMap[(node->children_m[j])->idx] >= 0;

      //#if VICUT_DEBUG
#if 0
      fprintf(stderr,"\nhasCounts=%d\n",(int)hasCounts);
      for (int j = 0; j < nChildren; j++)
      {
        if ( (node->children_m[j])->idx >= 0 )
          fprintf(stderr,"countMap[child%d(%d,%s)]=%d\n",j, (node->children_m[j])->idx, (node->children_m[j])->label.c_str(), countMap[(node->children_m[j])->idx]);
        else
          fprintf(stderr,"countMap[child%d(%d)]=%d\n",j, (node->children_m[j])->idx, countMap[(node->children_m[j])->idx]);
      }
      fprintf(stderr,"\n");
#endif

      if ( hasCounts )
      {
        countMap[i] = 0;
        int x;
        for (int j = 0; j < nChildren; j++)
          if ( (x=countMap[(node->children_m[j])->idx]) > 0 )
            countMap[i] += x;

        //#if VICUT_DEBUG
#if 0
        fprintf(stderr,"countMap[%d]=%d\n",i,countMap[i]);
#endif

        bool allCountsLEZ = true; // counts of all children <= 0
        for (int j = 0; j < nChildren; j++)
          allCountsLEZ &= countMap[(node->children_m[j])->idx] == 0;

        //#if VICUT_DEBUG
#if 0
        fprintf(stderr,"\nallCountsLEZ=%d\n",(int)allCountsLEZ);
#endif

        if (allCountsLEZ)
        {
          AcountMap[i] = NULL;
        }
        else
        {
          //#if VICUT_DEBUG
#if 0
          for (int j = 0; j < nChildren; j++)
          {
            fprintf(stderr,"AcountMap[node->children_m[%d]]: ",j);
            if ( AcountMap[(node->children_m[j])->idx] != NULL )
              for ( int k = 0; k < nuAnn; ++k )
                fprintf(stderr,"%d ",AcountMap[(node->children_m[j])->idx][k]);
            else
              fprintf(stderr,"NULL");
            fprintf(stderr,"\n");
          }
#endif

          AcountMap[i] = (int*)calloc(nuAnn, sizeof(int));

          for (int j = 0; j < nChildren; j++)
          {
            int chIdx = (node->children_m[j])->idx;
            //fprintf(stderr,"chIdx: %d\n", chIdx);
            for ( int k = 0; k < nuAnn; ++k )
            {
              if ( countMap[chIdx] > 0 )
              {
                //fprintf(stderr,"k: %d\n", k);
                AcountMap[i][k] += AcountMap[chIdx][k];
              }
            }
          }

          //#if VICUT_DEBUG
#if 0
          fprintf(stderr,"AcountMap[i]: ");
          for ( int k = 0; k < nuAnn; ++k )
            fprintf(stderr,"%d ",AcountMap[i][k]);
          fprintf(stderr,"\n\n");
#endif

          for (int j = 0; j < nChildren; j++)
          {
            int chIdx = (node->children_m[j])->idx;
            if ( AcountMap[chIdx] )
            {
              free( AcountMap[chIdx] );
              AcountMap[chIdx] = NULL;
            }
          }
        }
      }
      else
      {
        fprintf(stderr,"Error in %s at line %d: All children are missing count data\n",
                __FILE__,__LINE__);
        exit(1);
      }

      int nSpp = 0; // counting how many species are present in the current cluster
      if ( countMap[i] )
      {
        p = countMap[i]/(double)nAnnLeaves;
        qMap[i] = -H(p);

#if VICUT_DEBUG
        printf("p=%d/%d=%.3f\t-H(p)=%.3f\n", countMap[i], nAnnLeaves, p, qMap[i]);
        //       printf("%*s\tFreq\tTotal\t2H(Freq/nAnn)\n", maxAnnStrLen, "Label");
#endif

        for ( int j = 0; j < nuAnn; ++j )
        {
          if ( AcountMap[i][j] )
          {
            nSpp++;
            double h = 2*H( AcountMap[i][j] / (double)nAnnLeaves );
            qMap[i] += h;

#if VICUT_DEBUG
            // label frequency total 2H(freq/total)
            fprintf(stderr,"%d\t%.3f\n", AcountMap[i][j], Acount[j], h);
            //printf("p[%d]=A[%d]/n=%d/%d=%.3f\t2*H(A[%d]/n)=%.3f\tq[%d]=%.3f\n",
            //	 j,j,AcountMap[i][j], nAnnLeaves,(AcountMap[i][j] / (double)nAnnLeaves),j,(2*H( AcountMap[i][j] / (double)nAnnLeaves )),i,qMap[i]);
#endif
          }
          else
          {
#if VICUT_DEBUG
            // label frequency total 2H(freq/total)
            fprintf(stderr,"%d\t0\n", Acount[j]);
#endif
          }
        }

#if VICUT_DEBUG
        fprintf(stderr,"%d\t%.3f\t[= -H(p)]\n", countMap[i], nAnnLeaves, -H(p));
#endif
      }
      else
      {
        qMap[i] = 0;
      }

      double qChildren = 0;
      for (int j = 0; j < nChildren; j++)
        qChildren += qMap[(node->children_m[j])->idx];

#if VICUT_DEBUG
      fprintf(stderr,"\nq[%d]=%.3f\tqChildren=%.3f\n",i,qMap[i],qChildren);
#endif

      if ( qMap[i] <= qChildren ) // high mode
        //if ( qMap[i] < qChildren )
      {
        minCut[i] = 1;
        node = idx2node[i];
        resetMinCut(node, minCut); // set minCut on all internal nodes of the subtree rooted at nt[i] to 0
      }
      else
      {
        qMap[i] = qChildren;
      }

#if VICUT_DEBUG
      fprintf(stderr,"minCut[%d]=%d\tqMap[%d]=%.3f\n",i, minCut[i], i, qMap[i]);
#endif

      // finding the highest frequency annotation index
      int j1 = -1, j2 = -1;
      int *r = NULL;
      if ( AcountMap[i] )
      {
        r = rankArray( AcountMap[i], nuAnn );

        for ( int j = 0; j < nuAnn; ++j )
        {
          if ( r[j] == (nuAnn-1) )
            j1 = j;
          else if ( r[j] == (nuAnn-2) )
            j2 = j;
        }
      }

      i--;
    } // end of while ( ) loop

    free(AcountMap[minIdx]);

    if ( verbose )
      printf("\n");

#if VICUT_DEBUG
    fprintf(stderr,"minCut: ");
    pair<int, bool> pi;
    BOOST_FOREACH(pi, minCut)
    {
      if ( pi.second )
        fprintf(stderr,"%d ",pi.first);
    }
    fprintf(stderr,"\n\n");
#endif

    double viDist = 0; // VI distance
    for ( int i = 0; i < nLeaves; ++i )
    {
      if ( minCut[i] )
      {
        nodeCut.push_back(i);
        viDist += qMap[i];
      }
    }

#if VICUT_DEBUG
    fprintf(stderr,"leaf viDist=%.3f\n",viDist);
#endif

    for( int i = minIdx; i < 0; i++ )
    {
      if ( minCut[i] )
      {
        nodeCut.push_back(i);
        viDist += qMap[i];
      }
    }
#if VICUT_DEBUG
    fprintf(stderr,"all viDist=%.3f\tHA=%.3f\tfinal viDist=%.3f\n",viDist, HA, viDist - HA);
#endif

    // clean up
    free(Acount);

    return viDist - HA;
}
