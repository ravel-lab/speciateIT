#ifndef NEWICK
#define NEWICK 1

// Code for loading and manipulating Newick trees
// See: http://evolution.gs.washington.edu/phylip/newick_doc.html

// Author: Adam Phillippy
// Modifications: Pawel Gajer

/*
Copyright (C) 2016 Pawel Gajer (pgajer@gmail.com), Adam Phillippy and Jacques Ravel jravel@som.umaryland.edu

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

#include <vector>
#include <map>
#include <string>

#include "IOCppUtilities.hh"

using namespace std;

class NewickNode_t
{
public:
  NewickNode_t(NewickNode_t * parent = NULL);
  ~NewickNode_t();

  /// routines
  NewickNode_t * addChild();

  /// variables
  NewickNode_t * parent_m;
  vector<NewickNode_t *> children_m;
  double branch_length;
  std::string label;
  int idx;          // numeric index; positive for leaves, negative for internal nodes.
  int depth_m;      // node's depth in a tree; Root's depth is 0
  int model_idx;    // when node labels correspond to model labels model_idx is the index of the node label in MarkovChains2_t's modelIds_m
};

// this is used in getSppProfs() and getTx() routines for taxonomic assignment
class sppProf_t
{
public:
  NewickNode_t *node;      // cut-node
  NewickNode_t *ancNode;   // ancestral node with reference sequences (can be the same as 'node')
  map<string,int> sppFreq; // frequencies of sppecies found in the ancestral subtree
};

class NewickTree_t
{
public:
  NewickTree_t();
  ~NewickTree_t();

  bool loadTree(const char *file);
  void loadFullTxTree(const char *file);
  void loadFullTxTree2(const char *file);

  void rmNodesWith1child();
  void rmLeaf( string &s );

  void writeTree(FILE * fp);
  void writeTree(FILE * fp, NewickNode_t *node);
  void writeTree(FILE * fp, const string &nodeLabel);
  void printTree(bool withIdx=false, const char *indStr="  "); // prints tree to stdout; different depths are differen indentation level

  void modelIdx( vector<string> &modelIds );

  void indexNewickNodes(map<int, NewickNode_t*> &idx2node);

  char **leafLabels();
  void leafLabels(NewickNode_t *_node, vector<string> &leaves);
  void leafLabels(NewickNode_t *_node, set<string> &leaves);
  void leafLabels(NewickNode_t *_node, vector<NewickNode_t *> &leaves);
  void leafLabels(NewickNode_t *_node, vector<NewickNode_t *> &leaves, NewickNode_t *selNode);

  void assignIntNodeNames();

  int getDepth();
  int getNleaves() { return nLeaves_m; }
  int getMinIdx() { return minIdx_m; }
  int leafCount() { return nLeaves_m; }
  void decrementMinIdx() { minIdx_m--; }
  void incrementLeafCount() { nLeaves_m++; }

  void saveCltrMemb( const char *outFile,
		     vector<int> &nodeCut,
		     int *annIdx,
		     map<int, string> &idxToAnn);

  void saveCltrMemb2( const char *outFile,
		      vector<int> &nodeCut,
		      int *annIdx,
		      map<int, string> &idxToAnn);
  void saveCltrAnnStats( const char *outFile,
			 vector<int> &nodeCut,
			 int *annIdx,
			 map<int, string> &idxToAnn);

  vector<string> saveNAcltrAnnStats( const char *outFile,
				     vector<int> &nodeCut,
				     int *annIdx,
				     map<int, string> &idxToAnn);

  void saveNAcltrAnnStats( const char *outFileA, // file to with cluster stats of clusters with at least minNA query sequences
			   const char *outFileB, // file to with cluster stats of clusters with less than minNA query sequences
			   const char *outFileC, // file with reference IDs of sequences for clusters in outFileA
			   const char *outFileD, // file with taxons, for clusters in outFileA, with the number of sequences >= minAnn
			   vector<int> &nodeCut,
			   int *annIdx,
			   map<int, string> &idxToAnn,
			   int minNA,
			   int minAnn);

  void saveNAtaxonomy( const char *outFile,
		       const char *logFile,
		       const char *genusOTUsFile,
		       vector<int> &nodeCut,
		       int *annIdx,
		       map<int, string> &idxToAnn,
		       int minNA,
		       map<string, string> &txParent);

  void printGenusTrees( const char *outDir,
			vector<int> &nodeCut,
			int *annIdx,
			int minNA,
			map<int, string> &idxToAnn,
			map<string, string> &txParent);

  void getSppProfs( vector<sppProf_t*> &sppProfs,
                    vector<int> &nodeCut,
		    int *annIdx,
		    int minNA,
		    map<int, string> &idxToAnn,
		    map<string, string> &txParent);

  void getAncRoots(vector<sppProf_t*> &sppProfs,
		   vector<sppProf_t*> &ancRootNodes);

  void rmTxErrors( vector<sppProf_t*> &ancRootProfs,
		   int depth,
		   int *annIdx,
		   map<int, string> &idxToAnn);

  void printTx( const char *outFile,
                vector<sppProf_t*> &ancRootProfs);

  NewickNode_t * root() { return root_m; }
  void setRoot( NewickNode_t *node ) { root_m = node; }

  void txSet2txTree( strSet_t &tx2seqIDs );

  void updateLabels( strSet_t &tx2seqIDs );

  void inodeTx( const char *fullTxFile, map<string, string> &inodeTx );


private:
  string strMap2Str(map<string,int> &sppFreq);

  void readTaxtable(const char *taxtableFile,
		    map<string, string> &txParent);

  int cltrSpp(NewickNode_t *_node,
	      int *annIdx,
	      map<int, string> &idxToAnn,
	      map<string,int> &sppFreq);

  string maxAnn( NewickNode_t *node,
		 int *annIdx,
		 map<int, string> &idxToAnn);

  NewickNode_t* findAnnNode(std::string label);

  bool hasSelSpp(NewickNode_t *_node,
		 vector<NewickNode_t*> &v);

  NewickNode_t * root_m;
  int nLeaves_m; // number of leaves
  int minIdx_m;  // minimum value of node index
  map<int, NewickNode_t*> idx2node_m;
};

NewickTree_t *readNewickTree( const char *filename );

#endif
