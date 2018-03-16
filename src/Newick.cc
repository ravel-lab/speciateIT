// Author: Adam Phillippy

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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <queue>
#include <deque>
#include <string>

#include "IOCUtilities.h"
#include "IOCppUtilities.hh"
#include "Newick.hh"

using namespace std;


#define DEBUG 0


//--------------------------------------------- NewickNode_t() -----
NewickNode_t::NewickNode_t(NewickNode_t * parent)
  : parent_m(parent)
{
  #if DEBUG
  fprintf(stderr,"in NewickTree_t(parent)\t(parent==NULL)=%d\n",(int)(parent==NULL));
  #endif

  if (parent==NULL)
  {
    #if DEBUG
    fprintf(stderr,"in if(parent==NULL)\n");
    #endif
    depth_m = -1;
  }
  else
  {
    #if DEBUG
    fprintf(stderr,"in else\n");
    #endif
    depth_m = parent->depth_m+1;
  }

  #if DEBUG
  fprintf(stderr,"leaving NewickTree_t(parent)\tdepth_m: %d\n", depth_m);
  #endif

}

//--------------------------------------------- ~NewickNode_t() -----
NewickNode_t::~NewickNode_t()
{
  // for (unsigned int i = 0; i < children_m.size(); i++)
  // {
  //   delete children_m[i];
  // }
}


//--------------------------------------------- stitch -----
// make children on the input node the children of the current node
// the children of the input node have the current node as a parent
void NewickNode_t::stitch( NewickNode_t *node )
{
  #if DEBUG
  fprintf(stderr,"in NewickNode_t::stitch()\n");
  #endif

  children_m = node->children_m;

  int nChildren = children_m.size();
  for ( int i = 0; i < nChildren; i++ )
  {
    children_m[i]->parent_m = this;
  }

  #if DEBUG
  fprintf(stderr,"leaving NewickNode_t::stitch()\n");
  #endif
}

//--------------------------------------------- addChild -----
NewickNode_t * NewickNode_t::addChild()
{
  #if DEBUG
  fprintf(stderr,"in addChild()\n");
  #endif

  NewickNode_t * child = new NewickNode_t(this);


  #if DEBUG
  fprintf(stderr,"before children_m.size()\n");
  //fprintf(stderr,"children_m.size()=%d\n",(int)children_m.size());
  #endif

  children_m.push_back(child);

  #if DEBUG
  fprintf(stderr,"leaving addChild()\n");
  #endif

  return child;
}

//--------------------------------------------- NewickTree_t -----
NewickTree_t::NewickTree_t()
  : root_m(NULL), nLeaves_m(0), minIdx_m(-1)
{
  #if DEBUG
  fprintf(stderr,"in NewickTree_t()\n");
  #endif
}

//--------------------------------------------- ~NewickTree_t -----
NewickTree_t::~NewickTree_t()
{
  delete root_m;
}

//--------------------------------------------- rmLeaf -----
// remove leaf with label s
void NewickTree_t::rmLeaf( string &s )
{
  queue<NewickNode_t *> bfs2;
  bfs2.push(root_m);
  int numChildren;
  NewickNode_t *pnode;
  NewickNode_t *node;
  bool go = true;

  while ( !bfs2.empty() && go )
  {
    node = bfs2.front();
    bfs2.pop();

    numChildren = node->children_m.size();
    if ( numChildren )
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs2.push(node->children_m[i]);
      }
    }
    else if ( node->label == s )
    {
      pnode = node->parent_m;
      // identify which element of pnode->children_m vector is s
      int n = pnode->children_m.size();
      for ( int i = 0; i < n; ++i )
	if ( pnode->children_m[i]->label == s )
	{
	  pnode->children_m.erase ( pnode->children_m.begin() + i );
	  go = false;
	  break;
	}
    }

  } // end of while() loop
}


//--------------------------------------------- rmNodesWith1child -----
// remove nodes with only one child
void NewickTree_t::rmNodesWith1child()
{
  queue<NewickNode_t *> bfs2;
  bfs2.push(root_m);
  int numChildren;
  NewickNode_t *pnode;
  NewickNode_t *node;

  while ( !bfs2.empty() )
  {
    node = bfs2.front();
    bfs2.pop();

    while ( node->children_m.size() == 1 )
    {
      // finding node in a children array of the parent node
      pnode = node->parent_m;
      numChildren = pnode->children_m.size();

      #if DEBUG_LFTT
      fprintf(stderr, "\n\nNode %s has only one child %s; %s's parent is %s with %d children\n",
	      node->label.c_str(), node->children_m[0]->label.c_str(), node->label.c_str(),
	      pnode->label.c_str(), numChildren);
      #endif
      int i;
      for ( i = 0; i < numChildren; i++)
      {
	if ( pnode->children_m[i] == node )
	  break;
      }

      if ( i == numChildren )
      {
	fprintf(stderr, "ERROR in %s at line %d: node %s cannot be found in %s\n",
		__FILE__, __LINE__,(node->label).c_str(), (pnode->label).c_str());
	exit(1);
      }

      #if DEBUG_LFTT
      fprintf(stderr, "%s is the %d-th child of %s\n",
	      node->label.c_str(), i, pnode->label.c_str());

      fprintf(stderr, "%s children BEFORE change: ", pnode->label.c_str());
      for ( int j = 0; j < numChildren; j++)
	fprintf(stderr, "%s  ", pnode->children_m[j]->label.c_str());
      #endif

      node = node->children_m[0];
      pnode->children_m[i] = node;
      node->parent_m = pnode;

      #if DEBUG_LFTT
      fprintf(stderr, "\n%s children AFTER  change: ", pnode->label.c_str());
      for ( int j = 0; j < numChildren; j++)
	fprintf(stderr, "%s  ", pnode->children_m[j]->label.c_str());
      #endif
    }

    numChildren = node->children_m.size();
    if ( numChildren )
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs2.push(node->children_m[i]);
      }
    }
  }
}


//--------------------------------------------- loadTree -----
bool NewickTree_t::loadTree(const char *file)
{
  #define DEBUGLT 0

  #if DEBUGLT
  fprintf(stderr, "in NewickTree_t::loadTree()\n");
  #endif

  FILE * fp = fOpen(file, "r");

  char str[256];
  NewickNode_t * cur = new NewickNode_t();
  int LINE = 1;
  char c;
  bool DONE = false;
  int leafIdx = 0;
  int count = 0;

  while ((c = fgetc(fp)) != EOF)
  {
    #if DEBUGLT
    fprintf(stderr,"c=|%c|\tcount=%d\t(cur==NULL)=%d\n", c, count,(int)(cur==NULL));
    #endif
    count++;

    if (c == ';')  { DONE = true;}
    else if (c == '\n') { LINE++; }
    else if (c == ' ') { }
    else if (DONE)
    {
      cerr << "ERROR: Unexpected data after ; on line " << LINE << endl;
      return false;
    }
    else if (c == '(')
    {
      cur = cur->addChild();
      if (root_m == NULL) { root_m = cur; }
    }
    else if (c == ')') // '):bl,' or '):bl' or ')lable:bl'
    {
      c = fgetc(fp);

      #if DEBUGLT
      fprintf(stderr,"in ) next c=%c\n", c);
      fprintf(stderr,"in ) minIdx_m=%d\n", minIdx_m);
      #endif

      cur->idx = minIdx_m;
      minIdx_m--;

      if ( c == ';' )
      {
		DONE = 1;
		break;
      }
      else if ( c == ':' )
      {
		fscanf(fp, "%lf", &cur->branch_length);
      }
      else if ( isdigit(c) )
      {
		double bootVal;
		fscanf(fp, "%lf:%lf", &bootVal, &cur->branch_length);
      }
      else if ( isalpha(c) )  // label:digit
      {
        #if DEBUGLT
		fprintf(stderr,"in )alpha=%c\n", c);
        #endif

		char *brkt;
		ungetc(c, fp);
		int i = 0;
		while ( isalnum(c) || c==':' || c=='.' || c=='_' || c=='-' )
		{
		  c = fgetc(fp);
		  str[i++] = c;
		}
		str[--i] = '\0';
		if (c != ',') { ungetc(c, fp); }

		strtok_r(str, ":", &brkt);
		cur->label = string(str);
		cur->branch_length = 0; //atof(str);

        #if DEBUGLT
		fprintf(stderr, "label=%s\tbranch_length=%.2f\n", cur->label.c_str(), cur->branch_length);
        #endif
      }
      else //if (c != ':' && !isdigit(c) && c !=';')
      {
		fprintf(stderr, "ERROR in %s at %d: Unexpected data, %c, after ')' in %s on line %d\n",
				__FILE__, __LINE__, c, file, LINE);
		return 0;
      }

      cur = cur->parent_m;

      c = fgetc(fp);
      if (c != ',') { ungetc(c, fp); }
    }
    else // 'name:bl' or 'name:bl,' - leaf
    {
      #if DEBUGLT
      fprintf(stderr,"in name:bl=%c\n", c);
      #endif

      cur = cur->addChild();
      nLeaves_m++;
      cur->label += c;

      while ((c = fgetc(fp)) != ':')
      {
        cur->label += c;
      }

      #if DEBUGLT
      fprintf(stderr,"in name:bl label=%s\n", cur->label.c_str());
      #endif

      fscanf(fp, "%lf", &cur->branch_length);

      cur->idx = leafIdx;
      leafIdx++;

      cur = cur->parent_m;

      c = fgetc(fp);
      if (c != ',') { ungetc(c, fp); }
    }
  } // end of while()


  minIdx_m++;

  if (!root_m)
  {
    fprintf(stderr,"ERROR in %s at line %d: Newick root is null!\n",__FILE__,__LINE__);
    exit(EXIT_FAILURE);
  }

  return true;
}

//--------------------------------------------- writeNewickNode -----
static void writeNewickNode(NewickNode_t * node,
                            FILE * fp,
                            int depth)
{
  if (!node) { return; }

  int numChildren = node->children_m.size();

  if (numChildren == 0)
  {
    if ( node->branch_length )
      fprintf(fp, "%s:%0.20lf", node->label.c_str(), node->branch_length);
    else
      fprintf(fp, "%s:0.1", node->label.c_str());
  }
  else
  {
    //fprintf(fp, "(\n");
    fprintf(fp, "(");

    vector<NewickNode_t *> goodChildren;

    for (int i = 0; i < numChildren; i++)
    {
      goodChildren.push_back(node->children_m[i]);
    }

    int numGoodChildren = goodChildren.size();
    for (int i = 0; i < numGoodChildren; i++)
    {
      writeNewickNode(goodChildren[i], fp, depth+1);

      if (i != numGoodChildren-1) { fprintf(fp, ","); }
      //fprintf(fp, "\n");
    }

    //for (int i = 0; i < depth; i++) { fprintf(fp, " "); }

    if ( node->label.empty() )
    {
      if ( node->branch_length)
	fprintf(fp, "):%0.20lf", node->branch_length);
      else
	fprintf(fp, "):0.1");
    }
    else
    {
      if ( node->branch_length)
	fprintf(fp, ")%s:%0.20lf", node->label.c_str(), node->branch_length);
      else
	fprintf(fp, ")%s:0.1", node->label.c_str());
    }

    if (depth == 0) { fprintf(fp, ";\n"); }
  }
}

//--------------------------------------------- writeTree -----
void NewickTree_t::writeTree(FILE * fp)
{
  if (root_m == NULL)
  {
    fprintf(stderr, "ERROR: Root is NULL!\n");
  }
  writeNewickNode(root_m, fp, 0);
  fflush(fp);
}

//--------------------------------------------- writeTree -----
void NewickTree_t::writeTree(FILE * fp, NewickNode_t *node)
{
  if (root_m == NULL)
  {
    fprintf(stderr, "ERROR: Root is NULL!\n");
  }
  writeNewickNode(node, fp, 0);
  fflush(fp);
}

//--------------------------------------------- writeTree -----
void NewickTree_t::writeTree(FILE * fp, const string &nodeLabel)
{
  if (root_m == NULL) fprintf(stderr, "ERROR: Root is NULL!\n");

  queue<NewickNode_t *> bfs;
  bfs.push(root_m);
  NewickNode_t *node;
  int numChildren;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    if ( node->label==nodeLabel )
    {
      writeNewickNode(node, fp, 0);
      fflush(fp);
      break;
    }

    numChildren = node->children_m.size();
    if ( numChildren )
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }
}

//--------------------------------------------- getDepth -----
// Determine the depth of a tree
int NewickTree_t::getDepth()
{
  deque<NewickNode_t *> dfs;
  dfs.push_front(root_m);
  NewickNode_t *node;

  deque<int> depth;
  int maxDepth = 0;
  depth.push_front(maxDepth);

  while ( !dfs.empty() )
  {
    node = dfs.front();
    dfs.pop_front();

    int d = depth.front();
    depth.pop_front();

    if ( maxDepth < d )
      maxDepth = d;

    int numChildren = node->children_m.size();
    if ( numChildren>0 )
    {
      for (int i = 0; i < numChildren; i++)
      {
	dfs.push_front(node->children_m[i]);
	depth.push_front(d+1);
      }

      //if ( maxDepth < d+1 )
      // maxDepth = d+1;
    }
  }

  return maxDepth;
}

//--------------------------------------------- printTree -----
//
// prints tree to stdout with different depths at differen indentation levels
// indStr - indentation string
void NewickTree_t::printTree(bool withIdx, const char *indStr)
{
    if (root_m == NULL)
    {
      fprintf(stderr, "ERROR: Root is NULL!\n");
    }

    // Determine the depth of the tree
    int maxDepth = getDepth();

    printf("\nmaxDepth: %d\n\n", maxDepth);

    maxDepth++;

    // Set up indent array of indentation strings
    //const char *indStr="  ";
    int indStrLen = strlen(indStr);

    char **indent = (char **)malloc(maxDepth * sizeof(char*));
    indent[0] = (char*)malloc(sizeof(char)); // empty string
    indent[0][0] = '\0';
    for ( int i = 1; i < maxDepth; i++ )
    {
      indent[i] = (char*)malloc(indStrLen * (i+1) * sizeof(char));
      for ( int j = 0; j < i; j++ )
        for ( int k = 0; k < indStrLen; k++ )
          indent[i][j * indStrLen + k] = indStr[k];
      indent[i][indStrLen * i] = '\0';
    }

    //Depth first search
    deque<NewickNode_t *> dfs;
    dfs.push_front(root_m);
    NewickNode_t *node;

    while ( !dfs.empty() )
    {
      node = dfs.front();
      dfs.pop_front();

      int numChildren = node->children_m.size();
      if ( numChildren==0 ) // leaf
      {
        if ( withIdx )
          printf("%s%s(%d)\n",indent[node->depth_m],node->label.c_str(), node->idx);
        else
          //printf("%s%s\n",indent[node->depth_m],node->label.c_str());
          printf("(%d)%s%s\n",node->depth_m,indent[node->depth_m],node->label.c_str());
      }
      else
      {
        if ( withIdx )
        {
          if ( node->label.empty() )
            printf("%s*(%d)\n",indent[node->depth_m],node->idx);
          else
            printf("%s%s(%d)\n",indent[node->depth_m],node->label.c_str(), node->idx);
        }
        else
        {
          if ( node->label.empty() )
            printf("%s*\n",indent[node->depth_m]);
          else
            //printf("%s%s\n",indent[node->depth_m],node->label.c_str());
            printf("(%d)%s%s\n",node->depth_m,indent[node->depth_m],node->label.c_str());
        }

        for (int i = 0; i < numChildren; i++)
          dfs.push_front(node->children_m[i]);

      }
    }
}


#if 0
//--------------------------------------------- printTreeBFS -----
void NewickTree_t::printTreeBFS()
{
  if (root_m == NULL)
  {
    fprintf(stderr, "ERROR: Root is NULL!\n");
  }

  //Breath first search
  queue<NewickNode_t *> bfs;
  bfs.push(root_m);
  NewickNode_t *node;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    int numChildren = node->children_m.size();
    if ( numChildren==0 ) // leaf
    {
      printf("%s%s\n",indent[node->depth_m],node->label.c_str());
    }
    else
    {
      printf("%s*\n",indent[node->depth_m]);
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }
}
#endif

//--------------------------------------------- indexNewickNodes -----
void NewickTree_t::indexNewickNodes( map<int, NewickNode_t*> & idx2node)
/*
   populates a hash table of <node index> => <pointer to the node>
   returns number of leaves of the tree
*/
{
  if ( idx2node_m.size() == 0 )
  {
    //Breath first search
    queue<NewickNode_t *> bfs;
    bfs.push(root_m);
    NewickNode_t *node;

    while ( !bfs.empty() )
    {
      node = bfs.front();
      bfs.pop();

      idx2node_m[node->idx] = node;

      int numChildren = node->children_m.size();
      if ( numChildren != 0 ) // leaf
      {
	for (int i = 0; i < numChildren; i++)
	{
	  bfs.push(node->children_m[i]);
	}
      }
    }
  }

  idx2node = idx2node_m;
}

//--------------------------------------------- leafLabels -----
char ** NewickTree_t::leafLabels()
{
  char **leafLabel = (char**)malloc((nLeaves_m) * sizeof(char*));

  queue<NewickNode_t *> bfs;
  bfs.push(root_m);
  NewickNode_t *node;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    int numChildren = node->children_m.size();
    if ( numChildren==0 ) // leaf
    {
      leafLabel[node->idx] = strdup((node->label).c_str());
    }
    else
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }

  return leafLabel;
}

// ------------------------- assignIntNodeNames ------------------------
void NewickTree_t::assignIntNodeNames()
/*
  Asssign i1, i2, i3 etc names to internal nodes, where iK label corresponds to index -K
*/
{
  queue<NewickNode_t *> bfs;
  bfs.push(root_m);
  NewickNode_t *node;
  //int negIdx;   // holds negative index of a nod
  char str[16]; // holds negIdx

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    int numChildren = node->children_m.size();
    if ( numChildren ) // not leaf
    {
      //negIdx = -node->idx;
      sprintf(str,"%d", -node->idx);
      node->label = string("i") + string(str);

      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }
}




// ------------------------------ saveCltrMemb ------------------------
void NewickTree_t::saveCltrMemb( const char *outFile,
				 vector<int> &nodeCut,
				 int *annIdx,
				 map<int, string> &idxToAnn)
/*
  save clustering membership data to outFile
  output file format:

  <leaf ID> <cluster ID> <annotation str or NA if absent>

  Parameters:

  outFile  - output file
  nodeCut  - min-node-cut
  annIdx   - annIdx[i] = <0-based index of the annotation string assigned to the i-th
                          element of data table>
	     annIdx[i] = -1 if the i-th element has no annotation
  idxToAnn - idxToAnn[annIdx] = <annotation string associated with annIdx>
*/
{
  map<int, NewickNode_t*> idx2node;
  indexNewickNodes(idx2node);

  FILE *file = fOpen(outFile,"w");

  fprintf(file,"readId\tclstr\tannot\n");

  int n = nodeCut.size();
  for ( int i = 0; i < n; ++i )
  {
    if ( nodeCut[i] >= 0 )
    {
      fprintf(file,"%s\t%d\t%s\n",idx2node[nodeCut[i]]->label.c_str(),i,
	      idxToAnn[ annIdx[nodeCut[i]] ].c_str() );
    }
    else
    {
      queue<NewickNode_t *> bfs;
      bfs.push(idx2node[nodeCut[i]]);
      NewickNode_t *node;

      while ( !bfs.empty() )
      {
	node = bfs.front();
	bfs.pop();

	int numChildren = node->children_m.size();
	if ( numChildren==0 ) // leaf
	{
	  fprintf(file,"%s\t%d\t%s\n",node->label.c_str(),i,
		  idxToAnn[ annIdx[node->idx] ].c_str() );
	}
	else
	{
	  for (int i = 0; i < numChildren; i++)
	  {
	    bfs.push(node->children_m[i]);
	  }
	}
      }
    }
  }

  fclose(file);
}


// ------------------------------ saveCltrMemb2 ------------------------
void NewickTree_t::saveCltrMemb2( const char *outFile,
				  vector<int> &nodeCut,
				  int *annIdx,
				  map<int, string> &idxToAnn)
/*
  save clustering membership data to outFile
  output file format:

  <cluster_1>:
      <leaf1>  <annotation of leaf1>
      <leaf2>  <annotation of leaf2>
       ...
       ...

  Parameters:

  outFile  - output file
  nodeCut  - min-node-cut
  annIdx   - annIdx[i] = <0-based index of the annotation string assigned to the i-th
                          element of data table>
	     annIdx[i] = -1 if the i-th element has no annotation
  idxToAnn - idxToAnn[annIdx] = <annotation string associated with annIdx>
*/
{
  map<int, NewickNode_t*> idx2node;
  indexNewickNodes(idx2node);

  FILE *file = fOpen(outFile,"w");

  int n = nodeCut.size();
  for ( int i = 0; i < n; ++i )
  {
    fprintf(file,"Cluster %d:\n",i);

    if ( nodeCut[i] >= 0 )
    {
      fprintf(file,"\t%s\t%s\n",idx2node[nodeCut[i]]->label.c_str(),
	      idxToAnn[ annIdx[nodeCut[i]] ].c_str() );
    }
    else
    {
      queue<NewickNode_t *> bfs;
      bfs.push(idx2node[nodeCut[i]]);
      NewickNode_t *node;

      while ( !bfs.empty() )
      {
	node = bfs.front();
	bfs.pop();

	int numChildren = node->children_m.size();
	if ( numChildren==0 ) // leaf
	{
	  fprintf(file,"\t%s\t%s\n",node->label.c_str(),
		  idxToAnn[ annIdx[node->idx] ].c_str() );
	}
	else
	{
	  for (int i = 0; i < numChildren; i++)
	  {
	    bfs.push(node->children_m[i]);
	  }
	}
      }
    }
  }

  fclose(file);
}


// this is going to be used to sort map<string,int> in saveCltrAnnStats()
class sort_map
{
  public:
	string key;
	int val;
};

bool Sort_by(const sort_map& a ,const sort_map& b)
{
	return a.val > b.val;
}

// ------------------------------ saveCltrAnnStats --------------------------------
void NewickTree_t::saveCltrAnnStats( const char *outFile,
				     vector<int> &nodeCut,
				     int *annIdx,
				     map<int, string> &idxToAnn)
/*
  save annotation summary statistics of the clustering induced by nodeCut to outFile
  output file format:

  <cluster_1>:
      <annotation 1,1> <count1,1> <perc1,1>
      <annotation 1,2> <count1,2> <perc1,2>
       ...
       ...
  <cluster_2>:
      <annotation 2,1> <count2,1> <perc2,1>
      <annotation 2,2> <count2,2> <perc2,2>
       ...
       ...

  countX,Y is the count of leaves in cluster_X and Y-th annotation (sorted by count)
  percX,Y is the percentage of these leaves in the given cluster

  Parameters:

  outFile  - output file
  nodeCut  - min-node-cut
  annIdx   - annIdx[i] = <0-based index of the annotation string assigned to the i-th
                          element of data table>
	     annIdx[i] = -1 if the i-th element has no annotation
  idxToAnn - idxToAnn[annIdx] = <annotation string associated with annIdx>
*/
{
  map<int, NewickNode_t*> idx2node;
  indexNewickNodes(idx2node);

  FILE *file = fOpen(outFile,"w");

  int n = nodeCut.size();
  for ( int i = 0; i < n; ++i )
  {
    fprintf(file,"Cluster %d:\n",i);

    if ( nodeCut[i] >= 0 )
    {
      fprintf(file,"\t%s\t1\t100%%\n", idxToAnn[ annIdx[nodeCut[i]] ].c_str() );
    }
    else
    {
      map<string,int> annCount; // count of a given annotation string in the subtree of a given cut-node

      queue<NewickNode_t *> bfs;
      bfs.push(idx2node[nodeCut[i]]);
      NewickNode_t *node;

      while ( !bfs.empty() )
      {
	node = bfs.front();
	bfs.pop();

	int numChildren = node->children_m.size();
	if ( numChildren==0 ) // leaf
	{
	  annCount[ idxToAnn[ annIdx[node->idx] ] ]++;
	}
	else
	{
	  for (int i = 0; i < numChildren; i++)
	  {
	    bfs.push(node->children_m[i]);
	  }
	}
      }

      // printing annotation counts in the order of their counts
      map<string,int>::iterator it;
      vector< sort_map > v;
      vector< sort_map >::iterator itv;
      sort_map sm;
      double sum = 0;

      for (it = annCount.begin(); it != annCount.end(); ++it)
      {
	sm.key = (*it).first;
	sm.val = (*it).second;
	v.push_back(sm);
	sum += sm.val;
      }

      sort(v.begin(),v.end(),Sort_by);

      for (itv = v.begin(); itv != v.end(); ++itv)
	fprintf(file,"\t%s\t%d\t%.1f%%\n", ((*itv).key).c_str(), (*itv).val, 100.0 * (*itv).val / sum);
      if ( v.size() > 1 )
	fprintf(file,"\tTOTAL\t%d\t100%%\n", (int)sum);
    }
  }
  fclose(file);
}



// ------------------------------ saveNAcltrAnnStats ------------------------
vector<string> NewickTree_t::saveNAcltrAnnStats( const char *outFile,
						 vector<int> &nodeCut,
						 int *annIdx,
						 map<int, string> &idxToAnn)
/*
  It is a version of saveCltrAnnStats() that reports only clusters containig query (NA) elements.
*/
{
  map<int, NewickNode_t*> idx2node;
  indexNewickNodes(idx2node);

  FILE *file = fOpen(outFile,"w");

  vector<string> selAnn; // vector of taxons which are memembers of clusters that contain query seqeunces
  int n = nodeCut.size();
  for ( int i = 0; i < n; ++i )
  {
    if ( nodeCut[i] < 0 )
    {
      map<string,int> annCount; // count of a given annotation string in the subtree of a given cut-node

      queue<NewickNode_t *> bfs;
      bfs.push(idx2node[nodeCut[i]]);
      NewickNode_t *node;
      bool foundNA = false;

      while ( !bfs.empty() )
      {
	node = bfs.front();
	bfs.pop();

	int numChildren = node->children_m.size();
	if ( numChildren==0 ) // leaf
	{
	  annCount[ idxToAnn[ annIdx[node->idx] ] ]++;
	  if ( annIdx[node->idx] == -2) foundNA = true;
	}
	else
	{
	  for (int i = 0; i < numChildren; i++)
	  {
	    bfs.push(node->children_m[i]);
	  }
	}
      }

      // printing annotation counts in the order of their counts
      if ( foundNA )
      {
	map<string,int>::iterator it;
	vector< sort_map > v;
	vector< sort_map >::iterator itv;
	sort_map sm;
	double sum = 0;

	for (it = annCount.begin(); it != annCount.end(); ++it)
	{
	  sm.key = (*it).first;
	  sm.val = (*it).second;
	  selAnn.push_back(sm.key);
	  v.push_back(sm);
	  sum += sm.val;
	}

	sort(v.begin(),v.end(),Sort_by);

	fprintf(file,"Cluster %d:\n",i);
	for (itv = v.begin(); itv != v.end(); ++itv)
	  fprintf(file,"\t%s\t%d\t%.1f%%\n", ((*itv).key).c_str(), (*itv).val, 100.0 * (*itv).val / sum);
	if ( v.size() > 1 )
	  fprintf(file,"\tTOTAL\t%d\t100%%\n", (int)sum);
      }
    }
  }
  fclose(file);

  return selAnn;
}


// ------------------------------ saveNAcltrAnnStats ------------------------
void NewickTree_t::saveNAcltrAnnStats( const char *outFileA, // file to with cluster stats of clusters with at least minNA query sequences
				       const char *outFileB, // file to with cluster stats of clusters with less than minNA query sequences
				       const char *outFileC, // file with reference IDs of sequences for clusters in outFileA
				       const char *outFileD, // file with taxons, for clusters in outFileA, with the number of sequences >= minAnn
				       vector<int> &nodeCut,
				       int *annIdx,
				       map<int, string> &idxToAnn,
				       int minNA,
				       int minAnn)
/*
  It is a version of saveCltrAnnStats() that reports only clusters containig at
  least minNA query sequences to outFileA and less than minNA sequences to
  outFileB.

  output vector contains taxons with at least minAnn elements (in a single cluster).
*/
{
  map<int, NewickNode_t*> idx2node;
  indexNewickNodes(idx2node);

  FILE *fileA = fOpen(outFileA,"w");
  FILE *fileB = fOpen(outFileB,"w");
  FILE *fileC = fOpen(outFileC,"w");

  int nCltrs = 0; // number of clusters with the number of query seq's >= minNA
  int nTx    = 0; // number of selected taxons from the above clusters with the number of reads >= minAnn
  int nRef   = 0; // number of reference sequence from the above clusters with taxons with the number of reads >= minAnn

  map<string,bool> selTx; // hash table of selected taxons as in nTx; nTx = number of elements in selTx

  int n = nodeCut.size();
  for ( int i = 0; i < n; ++i )
  {
    if ( nodeCut[i] < 0 ) // internal node
    {
      map<string,int> annCount; // count of a given annotation string in the subtree of a given cut-node

      queue<NewickNode_t *> bfs;
      bfs.push(idx2node[nodeCut[i]]);
      NewickNode_t *node;
      int naCount = 0;

      while ( !bfs.empty() )
      {
	node = bfs.front();
	bfs.pop();

	int numChildren = node->children_m.size();
	if ( numChildren==0 ) // leaf
	{
	  annCount[ idxToAnn[ annIdx[node->idx] ] ]++;
	  if ( annIdx[node->idx] == -2) naCount++;
	}
	else
	{
	  for (int i = 0; i < numChildren; i++)
	  {
	    bfs.push(node->children_m[i]);
	  }
	}
      }

      // printing annotation counts in the order of their counts
      if ( naCount >= minNA )
      {
	nCltrs++;

	map<string,int>::iterator it;
	vector< sort_map > v;
	vector< sort_map >::iterator itv;
	sort_map sm;
	double sum = 0;

	for (it = annCount.begin(); it != annCount.end(); ++it)
	{
	  sm.key = (*it).first;
	  sm.val = (*it).second;
	  v.push_back(sm);
	  sum += sm.val;
	}

	sort(v.begin(),v.end(),Sort_by);

	fprintf(fileA,"Cluster %d:\n",i);
	for (itv = v.begin(); itv != v.end(); ++itv)
	  fprintf(fileA,"\t%s\t%d\t%.1f%%\n", ((*itv).key).c_str(), (*itv).val, 100.0 * (*itv).val / sum);
	if ( v.size() > 1 )
	  fprintf(fileA,"\tTOTAL\t%d\t100%%\n", (int)sum);


	// traverse the subtree again selecting sequences with annCount of the
	// corresponding taxons that have at least minAnn elements to be printed
	// in fileC

	queue<NewickNode_t *> bfs2;
	bfs2.push(idx2node[nodeCut[i]]);
	NewickNode_t *node;

	while ( !bfs2.empty() )
	{
	  node = bfs2.front();
	  bfs2.pop();

	  int numChildren = node->children_m.size();
	  if ( numChildren==0 ) // leaf
	  {
	    if ( annCount[ idxToAnn[ annIdx[node->idx] ] ] > minAnn && idxToAnn[ annIdx[node->idx] ] != "NA" )
	    {
	      selTx[ idxToAnn[ annIdx[node->idx] ] ] = true;
	      nRef++;
	      fprintf(fileC, "%s\t%s\n",
		      node->label.c_str(),
		      idxToAnn[ annIdx[node->idx] ].c_str() );
	    }
	  }
	  else
	  {
	    for (int i = 0; i < numChildren; i++)
	    {
	      bfs2.push(node->children_m[i]);
	    }
	  }
	}
      }
      else
      {
	map<string,int>::iterator it;
	vector< sort_map > v;
	vector< sort_map >::iterator itv;
	sort_map sm;
	double sum = 0;

	for (it = annCount.begin(); it != annCount.end(); ++it)
	{
	  sm.key = (*it).first;
	  sm.val = (*it).second;
	  v.push_back(sm);
	  sum += sm.val;
	}

	sort(v.begin(),v.end(),Sort_by);

	fprintf(fileB,"Cluster %d:\n",i);
	for (itv = v.begin(); itv != v.end(); ++itv)
	  fprintf(fileB,"\t%s\t%d\t%.1f%%\n", ((*itv).key).c_str(), (*itv).val, 100.0 * (*itv).val / sum);
	if ( v.size() > 1 )
	  fprintf(fileB,"\tTOTAL\t%d\t100%%\n", (int)sum);
      }
    }
  }
  fclose(fileA);
  fclose(fileB);
  fclose(fileC);


  // extracting only keys of selTx;
  FILE *fileD = fOpen(outFileD,"w");
  //vector<string> selTxV;  // vector of taxons which are memembers of clusters that contain query seqeunces
  for(map<string,bool>::iterator it = selTx.begin(); it != selTx.end(); ++it)
  {
    //selTxV.push_back(it->first);
    fprintf(fileD,"%s\n",(it->first).c_str());
  }
  fclose(fileD);

  nTx = selTx.size();

  fprintf(stderr,"\nNumber of clusters with the number of query seq's >= %d: %d\n",minNA, nCltrs);
  fprintf(stderr,"Number of selected taxons from the above clusters with the number of reads >= %d: %d\n",minAnn, nTx);
  fprintf(stderr,"Number of selected reference sequence from the above clusters with taxons with the number of reads >= %d: %d\n",minAnn, nRef);
}

// comparison between integers producing descending order
bool gr (int i,int j) { return (i>j); }

// ------------------------------ saveNAtaxonomy0 ------------------------
void NewickTree_t::saveNAtaxonomy0( const char *outFile,
				    const char *logFile, // recording sizes of two most abundant taxa. If there is only one (besides query sequences), recording the number of elements of the dominant taxon and 0.
				    vector<int> &nodeCut,
				    int *annIdx,
				    map<int, string> &idxToAnn,
				    int minAnn,
				    int &nNAs,
				    int &nNAs_with_tx,
				    int &tx_changed,
				    int &nClades_modified)
/*
  Assigning taxonomy to query sequences and modifying taxonomy of sequences with
  non-monolityc clades using majority vote.

  This is a special case of saveNAtaxonomy() that does taxonomy assignments of
  query sequences from clusters containing at least minAnn reference
  sequences. The taxonomy is based on the most abundant reference taxon.

  I just extended the algorithm so that taxonomy of all elements of the given
  cluster is set to the taxonomy of the cluster with the largest number of
  sequences. If there are ties, that is two taxonomies are most abundant, then
  there is no change.
*/
{
  map<int, NewickNode_t*> idx2node;
  indexNewickNodes(idx2node);

  map<string,int> otuFreq; // <taxonomic rank name> => number of OTUs already assigned to it

  FILE *file = fOpen(outFile,"w");
  FILE *logfile = fOpen(logFile,"w");

  nNAs = 0;
  nNAs_with_tx = 0;
  tx_changed = 0;
  nClades_modified = 0;

  int n = nodeCut.size();
  for ( int i = 0; i < n; ++i )
  {
    map<string,int> annCount; // count of a given annotation string in the subtree of a given node
    int naCount = 0;

    queue<NewickNode_t *> bfs;
    bfs.push(idx2node[nodeCut[i]]);
    NewickNode_t *node;

    while ( !bfs.empty() )
    {
      node = bfs.front();
      bfs.pop();

      int numChildren = node->children_m.size();
      if ( numChildren==0 ) // leaf
      {
	//if ( i==262 )printf("\n%s\t%s", node->label.c_str(), idxToAnn[ annIdx[node->idx] ].c_str());

	annCount[ idxToAnn[ annIdx[node->idx] ] ]++;
	if ( annIdx[node->idx] == -2) naCount++;
	if ( idxToAnn[ annIdx[node->idx] ] == "NA" ) nNAs++;
      }
      else
      {
	for (int j = 0; j < numChildren; j++)
	{
	  bfs.push(node->children_m[j]);
	}
      }
    }
    //if ( i==262 )printf("\n\n");

    //if ( naCount >= minNA && annCount.size() >= minAnn )
    if ( annCount.size() >= (unsigned)minAnn ) // Assigning taxonomy to query sequences and modifying taxonomy of sequences with non-monolityc clades using majority vote
    {
      string cltrTx;
      int maxCount = 0; // maximal count among annotated sequences in the cluster
      map<string,int>::iterator it;

      // Make sure there are not ties for max count

      // Copying counts into a vector, sorting it and checking that the first two
      // values (when sorting in the decreasing order) are not equal to each
      // other.
      vector<int> counts;

      for (it = annCount.begin(); it != annCount.end(); ++it)
      {
	if ( (*it).first != "NA" )
	  counts.push_back((*it).second);

	if ( (*it).first != "NA" && (*it).second > maxCount )
	{
	  maxCount = (*it).second;
	  cltrTx   = (*it).first;
	}
      }

      sort( counts.begin(), counts.end(), gr );

      #if 0
      if ( i==231 )
      {
	fprintf(stderr,"\ncltrTx: %s\n", cltrTx.c_str());
	fprintf(stderr,"\nCluster %d\nannCount\n",i);
	for (it = annCount.begin(); it != annCount.end(); ++it)
	{
	  fprintf(stderr,"\t%s\t%d\n", ((*it).first).c_str(), (*it).second);
	}
	fprintf(stderr,"\nsorted counts: ");
	for (int j = 0; j < counts.size(); j++)
	  fprintf(stderr,"%d, ", counts[j]);
	fprintf(stderr,"\n\n");
      }
      #endif

      // Traverse the subtree again assigning taxonomy to all query sequences
      int counts_size = counts.size();
      if ( counts_size==1 || ( counts_size>1 && counts[0] > counts[1] ) )
      {
	nClades_modified++;

	if ( counts_size > 1 )
	  fprintf(logfile, "%d\t%d\n", counts[0], counts[1]);
	else
	  fprintf(logfile, "%d\t0\n", counts[0]);

	queue<NewickNode_t *> bfs2;
	bfs2.push(idx2node[nodeCut[i]]);
	NewickNode_t *node;

	while ( !bfs2.empty() )
	{
	  node = bfs2.front();
	  bfs2.pop();

	  int numChildren = node->children_m.size();
	  if ( numChildren==0 ) // leaf
	  {
	    //if ( idxToAnn[ annIdx[node->idx] ] == "NA" && !cltrTx.empty() )
	    if ( idxToAnn[ annIdx[node->idx] ] !=cltrTx  && !cltrTx.empty() )
	    {
	      //idxToAnn[ annIdx[node->idx] ] = cltrTx;
	      fprintf(file, "%s\t%s\n", node->label.c_str(), cltrTx.c_str());
	      if ( idxToAnn[ annIdx[node->idx] ] == "NA" )
		nNAs_with_tx++;
	      else
		tx_changed++;
	    }
	  }
	  else
	  {
	    for (int j = 0; j < numChildren; j++)
	      bfs2.push(node->children_m[j]);
	  }
	} // end while ( !bfs2.empty() )
      } // end  if ( counts.size()==1 || ( counts.size()>1 && counts[0] > counts[1] ) )

    } // end if ( naCount >= minNA && annCount.size() >= minAnn )

  } // end   for ( int i = 0; i < n; ++i )

  fclose(file);
  fclose(logfile);
}

// ------------------------------ saveNAtaxonomy ------------------------
void NewickTree_t::saveNAtaxonomy( const char *outFile,
				   const char *logFile,
				   const char *genusOTUsFile,
				   vector<int> &nodeCut,
				   int *annIdx,
				   map<int, string> &idxToAnn,
				   int minNA,
				   map<string, string> &txParent)
/*
  Taxonomy assignments of query sequences from clusters containing either only
  query sequences or query sequences with at least minNA query sequences and
  minAnn reference sequences; In the last case the taxonomy is based on the most
  abundant reference taxon. In the previous case, the taxonomy is determined by
  the first ancestor with reference sequences. The content of reference sequences
  determines taxonomy. Here are some rules.

  If the ancestor contains reference sequences of the same species, then we
  assign this species name to the OTU.

  If the ref seq's are from the same genus, then we assign to the OTU the name of
  the most abundant ref sequence. A justification for this comes from an
  observation that often OTUs form a clade in the middle of which there are two
  subtrees with two different species of the same genus. I hypothesize that this
  is a clade of one species and the other one is an error. Therefore, the name of
  the OTU should be the name of the more abundant species, and the other species
  should be renamed to the other one.

  Alternatively, the name could be of the form
  <genus name>_<species1>.<species2>. ... .<speciesN>
  where species1, ... , speciesN are species present in the ancestral node

  If the ref seq's are from different genera, then we assign to the OTU a name
  <tx rank>_OTUi
  where tx rank is the taxonomic rank of common to all ref seq's of the ancestral subtree.

  In the extreme case of Bacteria_OTU I am merging
*/
{
  FILE *file = fOpen(outFile,"w");
  FILE *fh = fOpen(logFile,"w");
  FILE *genusFH = fOpen(genusOTUsFile,"w");

  //map<string,bool> selTx; // hash table of selected taxons as in nTx; nTx = number of elements in selTx
  map<int, NewickNode_t*> idx2node;
  indexNewickNodes(idx2node);

  map<string,int> otuFreq; // <taxonomic rank name> => number of OTUs already assigned to it

  int n = nodeCut.size();
  for ( int i = 0; i < n; ++i )
  {
    //if ( nodeCut[i] < 0 ) // internal node
    {
      map<string,int> annCount; // count of a given annotation string in the subtree of a given cut-node

      queue<NewickNode_t *> bfs;
      bfs.push(idx2node[nodeCut[i]]);
      NewickNode_t *node;
      int naCount = 0;

      while ( !bfs.empty() )
      {
	node = bfs.front();
	bfs.pop();

	int numChildren = node->children_m.size();
	if ( numChildren==0 ) // leaf
	{
	  annCount[ idxToAnn[ annIdx[node->idx] ] ]++;
	  if ( annIdx[node->idx] == -2) naCount++;
	}
	else
	{
	  for (int i = 0; i < numChildren; i++)
	  {
	    bfs.push(node->children_m[i]);
	  }
	}
      }

      if ( naCount >= minNA )
      {
	string cltrTx; // = "Unclassified"; // taxonomy for all query reads of the cluster with no annotated sequences

	if ( annCount.size() > 1 ) // there is at least one annotated element; setting cltrTx to annotation with the max frequency
	{
	  int maxCount = 0; // maximal count among annotated sequences in the cluster
	  map<string,int>::iterator it;

	  for (it = annCount.begin(); it != annCount.end(); ++it)
	  {
	    if ( (*it).first != "NA" && (*it).second > maxCount )
	    {
	      maxCount = (*it).second;
	      cltrTx   = (*it).first;
	    }
	  }
	}
	else // there are no annotated sequences in the cluster; we are looking for the nearest ancestral node with known taxonomy
	{
	  node = idx2node[nodeCut[i]]; // current node
	  node = node->parent_m; // looking at the parent node

	  map<string,int> sppFreq;
	  int nSpp = cltrSpp(node, annIdx, idxToAnn, sppFreq); // table of species frequencies in the subtree of 'node'
	  while ( nSpp==0 ) // if the parent node does not have any annotation sequences
	  {
	    node = node->parent_m; // move up to its parent node
	    nSpp = cltrSpp(node, annIdx, idxToAnn, sppFreq);
	  }

	  // Debug
	  // fprintf(stderr, "\n--- i: %d\tCluster %d\n", i, node->idx);
	  // fprintf(stderr,"Number of species detected in the Unassigned node: %d\n",nSpp);

	  // besides writing frequencies I am finding taxonomic parents of each species
	  map<string,int> genus;  // this is a map not vector so that I avoid duplication and find out frequence of a genus
	  map<string,int>::iterator it;
	  for ( it = sppFreq.begin(); it != sppFreq.end(); it++ )
	  {
	    //if ( (*it).first != "NA" )
	    //fprintf(fh, "%s\t%d\n", ((*it).first).c_str(), (*it).second);
	    if ( !txParent[(*it).first].empty() )
	      genus[ txParent[(*it).first] ] += (*it).second;
	  }

	  #if 0
	  // Checking of there is more than one genus present
	  // If so, we are going higher until we get only one taxonomic rank.
	  int nRks = genus.size();
	  map<string,int> txRk = genus; // taxonomic rank
	  map<string,int> pTxRk;        // parent taxonomic rank
	  while ( nRks>1 )
	  {
	    for ( it = txRk.begin(); it != txRk.end(); it++ )
	      if ( !txParent[(*it).first].empty() )
		pTxRk[ txParent[(*it).first] ]++;

	    // fprintf(fh, "\nTx Rank Frequencies\n");
	    // for ( it = pTxRk.begin(); it != pTxRk.end(); it++ )
	    //   fprintf(fh, "%s\t%d\n", ((*it).first).c_str(), (*it).second);

	    nRks = pTxRk.size();
	    txRk = pTxRk;
	    pTxRk.clear();

	    // Debug
	    // fprintf(stderr,"nRks: %d\n",nRks);
	  }
	  it = txRk.begin();
	  string otuRank = (*it).first;
	  #endif

	  // assign OTU ID
	  if ( genus.size() == 1 ) // assign species name of the most frequent OTU
	  {
	    int maxCount = 0; // maximal count among annotated sequences in the cluster
	    map<string,int>::iterator it;
	    for (it = sppFreq.begin(); it != sppFreq.end(); ++it)
	    {
	      if ( (*it).first != "NA" && (*it).second > maxCount )
	      {
		maxCount = (*it).second;
		cltrTx   = (*it).first;
	      }
	    }
	  }
	  else
	  {
	    // assigning the same OTU id to all OTUs with the same sppFreq profile
	    string sppProf;
	    map<string,int>::iterator it;
	    for (it = sppFreq.begin(); it != sppFreq.end(); ++it)
	    {
	      if ( (*it).first != "NA" )
	      {
		sppProf += (*it).first;
	      }
	    }

	    int maxCount = 0; // maximal count among genera
	    string otuGenus;
	    for (it = genus.begin(); it != genus.end(); ++it)
	    {
	      if ( (*it).first != "NA" && (*it).second > maxCount )
	      {
		maxCount = (*it).second;
		otuGenus = (*it).first;
	      }
	    }

	    it = otuFreq.find(sppProf);
	    int k = 1;
	    if ( it != otuFreq.end() )
	      k = it->second + 1;
	    otuFreq[sppProf] = k;

	    char kStr[8];
	    sprintf(kStr,"%d",k);
	    cltrTx = otuGenus + string("_OTU") + string(kStr);
	  }


	  // Writing sppFreq to log file
	  fprintf(fh,"\n==== %s\tAncestral node ID: %d ====\n\n",cltrTx.c_str(), node->idx);


	  fprintf(fh, "-- Species Frequencies\n");
	  int totalSppLen = 0;
	  for ( it = sppFreq.begin(); it != sppFreq.end(); it++ )
	  {
	    fprintf(fh, "%s\t%d\n", ((*it).first).c_str(), (*it).second);
	    totalSppLen += ((*it).first).length();
	  }

	  // For each species of sppFreq, determine its taxonomy parent and write them to the log file
	  fprintf(fh, "\n-- Genus Frequencies\n");
	  for ( it = genus.begin(); it != genus.end(); it++ )
	    fprintf(fh, "%s\t%d\n", ((*it).first).c_str(), (*it).second);

	  if ( genus.size() == 1 ) // for genus OTUs print OTU name, ancestral node ID and species names into genusOTUsFile
	  {
	    //fprintf(stderr, "%s\t genus.size(): %d\n", cltrTx.c_str(), (int)genus.size());

	    // number of quary sequences in the OTU
	    it = sppFreq.find("NA");
	    if ( it == sppFreq.end() )
	      fprintf(stderr,"ERROR: OTU %s does not contain any quary sequences\n",cltrTx.c_str());

	    // OTU_ID, number of query sequences in the OTU, ID of the ancestral node that contains ref sequences
	    fprintf(genusFH,"%s\t%d\t%d\t",cltrTx.c_str(), (*it).second, node->idx);

	    int n = sppFreq.size();
	    char *sppStr = (char*)malloc((totalSppLen+n+1) * sizeof(char));
	    sppStr[0] = '\0';
	    int foundSpp = 0;
	    for ( it = sppFreq.begin(); it != sppFreq.end(); it++ )
	    {
	      if ( (*it).first != "NA" && (*it).first != "Unclassified" )
	      {
		if ( foundSpp )
		  strncat(sppStr, ":",1);
		strncat(sppStr, ((*it).first).c_str(), ((*it).first).length());
		foundSpp = 1;
	      }
	    }
	    fprintf(genusFH,"%s\n", sppStr);
	    free(sppStr);
	  }

	  #if 0
	  // Higher taxonomic order frequencies
	  nRks = genus.size();
	  txRk = genus;  // taxonomic rank
	  pTxRk.clear(); // parent taxonomic rank
	  while ( nRks>1 )
	  {
	    for ( it = txRk.begin(); it != txRk.end(); it++ )
	      if ( !txParent[(*it).first].empty() )
		pTxRk[ txParent[(*it).first] ]++;

	    fprintf(fh, "\n-- Tx Rank Frequencies\n");
	    for ( it = pTxRk.begin(); it != pTxRk.end(); it++ )
	      fprintf(fh, "%s\t%d\n", ((*it).first).c_str(), (*it).second);

	    nRks = pTxRk.size();
	    txRk = pTxRk;
	    pTxRk.clear();
	  }
	  #endif

	  //Debugging
	  //fprintf(stderr,"OTU ID: %s\n",cltrTx.c_str());

	} // end if ( annCount.size() > 1 )

	// traverse the subtree again assigning taxonomy to all query sequences
	queue<NewickNode_t *> bfs2;
	bfs2.push(idx2node[nodeCut[i]]);
	NewickNode_t *node;

	while ( !bfs2.empty() )
	{
	  node = bfs2.front();
	  bfs2.pop();

	  int numChildren = node->children_m.size();
	  if ( numChildren==0 ) // leaf
	  {
	    if ( idxToAnn[ annIdx[node->idx] ] == "NA" )
	    {
	      fprintf(file, "%s\t%s\n", node->label.c_str(), cltrTx.c_str());
	    }
	  }
	  else
	  {
	    for (int i = 0; i < numChildren; i++)
	      bfs2.push(node->children_m[i]);
	  }
	} // end while ( !bfs2.empty() )

      } // if ( naCount >= minNA )

    } // end if ( nodeCut[i] < 0 ) // internal node

  } // end   for ( int i = 0; i < n; ++i )


  fclose(fh);
  fclose(genusFH);
  fclose(file);
}


// ------------------------------ cltrSpp ------------------------
// Traverses the subtree of the current tree rooted at '_node' and finds
// annotation strings (taxonomies) of reference sequences present in the subtree
// (if any). sppFreq is the frequency table of found species.
//
// Returns the size of sppFreq
//
int NewickTree_t::cltrSpp(NewickNode_t *_node,
			  int *annIdx,
			  map<int, string> &idxToAnn,
			  map<string,int> &sppFreq)
{
  queue<NewickNode_t *> bfs;
  bfs.push(_node);
  NewickNode_t *node;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    int numChildren = node->children_m.size();
    if ( numChildren==0 ) // leaf
    {
      sppFreq[ idxToAnn[ annIdx[node->idx] ] ]++;
    }
    else
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }

  return int(sppFreq.size());
}




// ------------------------------ maxAnn ------------------------
// traverse the subtree of the current tree rooted at '_node' and find annotation
// string with max frequency Return empty string if there are no annotation
// sequences in the subtree.
string NewickTree_t::maxAnn( NewickNode_t *_node,
			     int *annIdx,
			     map<int, string> &idxToAnn)

{
  map<string,int> annCount; // count of a given annotation string in the subtree of the given node
  queue<NewickNode_t *> bfs;
  bfs.push(_node);
  NewickNode_t *node;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    int numChildren = node->children_m.size();
    if ( numChildren==0 ) // leaf
    {
      annCount[ idxToAnn[ annIdx[node->idx] ] ]++;
    }
    else
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }

  string cltrTx; // = "Unclassified"; // taxonomy for all query reads of the cluster with no annotated sequences

  if ( annCount.size() > 1 ) // there is at least one annotated element; setting cltrTx to annotation with the max frequency
  {
    int maxCount = 0; // maximal count among annotated sequences in the cluster
    map<string,int>::iterator it;

    for (it = annCount.begin(); it != annCount.end(); ++it)
    {
      if ( (*it).first != "NA" && (*it).second > maxCount )
      {
	maxCount = (*it).second;
	cltrTx   = (*it).first;
      }
    }
  }

  return cltrTx;
}


// ------------------------------ findAnnNode ------------------------
// Find a node in lable 'label'
NewickNode_t* NewickTree_t::findAnnNode(std::string label)
{
  queue<NewickNode_t *> bfs;
  bfs.push(root_m);
  NewickNode_t *node = 0;
  NewickNode_t *searchedNode = 0;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    if ( node->label == label )
    {
      searchedNode = node;
      break;
    }

    int numChildren = node->children_m.size();
    if ( numChildren ) // not leaf
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }

  return searchedNode;
}


// ------------------------------ printGenusTrees ------------------------
// For each genus print a subtree rooted at the root nodes of the genus
// * to each genus assign a vector of species/OTU's cut-nodes generated by vicut
// * find a node with the smallest depth (from the root)
// * walk ancestors of this node until all OTU nodes are in the subtree of that node
// * print that tree to a file
void NewickTree_t::printGenusTrees( const char *outDir,
				    vector<int> &nodeCut,
				    int *annIdx,
				    int minNA,
				    map<int, string> &idxToAnn,
				    map<string, string> &txParent)
{
  //fprintf(stderr,"Entering printGenusTrees()\n");

  char s[1024];
  sprintf(s,"mkdir -p %s",outDir);
  system(s);

  // === Assign to each genus a vector of species/OTU nodes
  map<string, vector<NewickNode_t*> > genusSpp; // <genus> => vector of OTU/species node

  map<int, NewickNode_t*> idx2node;
  indexNewickNodes(idx2node);

  int n = nodeCut.size();
  for ( int i = 0; i < n; ++i )
  {
    //fprintf(stderr,"i: %d\r", i);

    if ( nodeCut[i] < 0 ) // internal node
    {
      map<string,int> annCount; // count of a given annotation string in the subtree of a given cut-node

      queue<NewickNode_t *> bfs;
      bfs.push(idx2node[nodeCut[i]]);
      NewickNode_t *node;
      int naCount = 0;

      while ( !bfs.empty() )
      {
	node = bfs.front();
	bfs.pop();

	int numChildren = node->children_m.size();
	if ( numChildren==0 ) // leaf
	{
	  annCount[ idxToAnn[ annIdx[node->idx] ] ]++;
	  if ( annIdx[node->idx] == -2) naCount++;
	}
	else
	{
	  for (int i = 0; i < numChildren; i++)
	  {
	    bfs.push(node->children_m[i]);
	  }
	}
      }

      if ( naCount >= minNA )
      {
	string cltrTx; // taxonomy for all query reads of the cluster with no annotated sequences

	if ( annCount.size() > 1 ) // there is at least one annotated element; setting cltrTx to annotation with the max frequency
	{
	  int maxCount = 0; // maximal count among annotated sequences in the cluster
	  map<string,int>::iterator it;

	  for (it = annCount.begin(); it != annCount.end(); ++it)
	  {
	    if ( (*it).first != "NA" && (*it).second > maxCount )
	    {
	      maxCount = (*it).second;
	      cltrTx   = (*it).first;
	    }
	  }

	  genusSpp[ txParent[ cltrTx ] ].push_back( idx2node[nodeCut[i]] );
	}
	else // there are no annotated sequences in the cluster; we are looking for the nearest ancestral node with known taxonomy
	{

	  //fprintf(stderr, "In if ( naCount >= minNA ) else\n");

	  node = idx2node[nodeCut[i]]; // current node
	  node = node->parent_m; // looking at the parent node

	  map<string,int> sppFreq;
	  int nSpp = cltrSpp(node, annIdx, idxToAnn, sppFreq); // table of species frequencies in the subtree of 'node'
	  while ( nSpp==0 ) // if the parent node does not have any annotation sequences
	  {
	    node = node->parent_m; // move up to its parent node
	    nSpp = cltrSpp(node, annIdx, idxToAnn, sppFreq);
	  }

	  //fprintf(stderr, "nSpp: %d\n", nSpp);

	  map<string,int> genus;  // this is a map not vector so that I avoid duplication and find out frequence of a genus
	  map<string,int>::iterator it;
	  for ( it = sppFreq.begin(); it != sppFreq.end(); it++ )
	  {
	    if ( !txParent[(*it).first].empty() )
	      genus[ txParent[(*it).first] ]++;
	  }

	  //fprintf(stderr, "genus.size(): %d\n", (int)genus.size());

	  if ( genus.size() == 1 ) //
	  {
	    it = genus.begin();
	    //fprintf(stderr,"genus: %s\n", (*it).first.c_str());
	    genusSpp[ (*it).first ].push_back( node );
	  }

	} // end if ( annCount.size() > 1 )

      } // if ( naCount >= minNA )

    } // end if ( nodeCut[i] < 0 ) // internal node

  } // end   for ( int i = 0; i < n; ++i )



  // === For each genus, find a node with the smallest depth (from the root)
  // === Walk this node until all OTU nodes are in the subtree of that node
  // === Print that tree to a file
  map<string, vector<NewickNode_t*> >::iterator itr;
  for ( itr = genusSpp.begin(); itr != genusSpp.end(); itr++ )
  {
    vector<NewickNode_t*> v = (*itr).second;

    int vn = (int)v.size();
    //Debugging
    // fprintf(stderr,"%s: ", (*itr).first.c_str());
    // for ( int j = 0; j < vn; j++ )
    //   fprintf(stderr,"(%d, %d)  ", v[j]->idx, v[j]->depth_m);
    // fprintf(stderr,"\n");

    NewickNode_t* minDepthNode = v[0];
    int minDepth = v[0]->depth_m;
    for ( int j = 1; j < vn; j++ )
    {
      if ( v[j]->depth_m < minDepth )
      {
	minDepth = v[j]->depth_m;
	minDepthNode = v[j];
      }
    }

    NewickNode_t* node = minDepthNode->parent_m;

    bool foundParentNode = hasSelSpp(node, v); // looking for a parent node of all elements of v
    while ( foundParentNode==false )
    {
      node = node->parent_m; // move up to its parent node
      foundParentNode = hasSelSpp(node, v);
    }

    // write subtree rooted in node to a file
    int n1 = strlen(outDir);
    int n2 =  ((*itr).first).length();
    char *outFile = (char*)malloc((n1+n2+1+5)*sizeof(char));
    strcpy(outFile,outDir);
    strcat(outFile,((*itr).first).c_str());
    strcat(outFile,".tree");
    FILE *fh = fOpen(outFile,"w");
    writeTree(fh, node);
    fclose(fh);

    // write leaf labels to a file
    vector<string> leaves;
    leafLabels(node, leaves);

    char *idsFile = (char*)malloc((n1+n2+1+5)*sizeof(char));
    strcpy(idsFile,outDir);
    strcat(idsFile,((*itr).first).c_str());
    strcat(idsFile,".ids");
    writeStrVector(idsFile, leaves);
  }
}


// ------------------------------ hasSelSpp ------------------------
// Traverses the subtree of the current tree rooted at '_node'
// and checks if that tree contains all elements of v
//
// Returns true if all elements of v are in the subtree of _node
//
bool NewickTree_t::hasSelSpp(NewickNode_t *_node,
			     vector<NewickNode_t*> &v)
{
  queue<NewickNode_t *> bfs;
  bfs.push(_node);
  NewickNode_t *node;
  int countFoundNodes = 0; // count of nodes of v found in the subtree of _node
  // table of indixes of nodes of v; return values of the table are not used
  // the table is used only to find is a node is in v
  map<int,bool> V;
  int n = (int)v.size();
  for ( int i = 0; i < n; i++ )
    V[v[i]->idx] = true;

  map<int,bool>::iterator it;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    it = V.find(node->idx);
    if ( it != V.end() )
      countFoundNodes++;

    int numChildren = node->children_m.size();
    if ( numChildren!=0 ) // internal node
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }

  return countFoundNodes == (int)v.size();
}



//--------------------------------------------- leafLabels -----
// Given a node in the tree, the routine gives leaf labels
// of the subtree rooted at 'node'.
void NewickTree_t::leafLabels(NewickNode_t *_node, vector<string> &leaves)
{
  leaves.clear();

  queue<NewickNode_t *> bfs;
  bfs.push(_node);
  NewickNode_t *node;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    int numChildren = node->children_m.size();
    if ( numChildren==0 ) // leaf
    {
      leaves.push_back( node->label );
    }
    else
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }
}


//--------------------------------------------- leafLabels -----
// Given a node in the tree, the routine returns a vector of pointes to leaf nodes
// of the subtree rooted at 'node'.
void NewickTree_t::leafLabels(NewickNode_t *_node, vector<NewickNode_t *> &leaves)
{
  queue<NewickNode_t *> bfs;
  bfs.push(_node);
  NewickNode_t *node;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    int numChildren = node->children_m.size();
    if ( numChildren==0 ) // leaf
    {
      leaves.push_back( node );
    }
    else
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }
}


//--------------------------------------------- leafLabels -----
// Given a node in the tree, the routine returns a vector of pointes to leaf nodes
// (excluding a selected node) of the subtree rooted at 'node'.
void NewickTree_t::leafLabels(NewickNode_t *_node, vector<NewickNode_t *> &leaves, NewickNode_t *selNode)
{
  queue<NewickNode_t *> bfs;
  bfs.push(_node);
  NewickNode_t *node;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    int numChildren = node->children_m.size();
    if ( numChildren==0 ) // leaf
    {
      if ( node != selNode )
	leaves.push_back( node );
    }
    else
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }
}


//--------------------------------------------- leafLabels -----
// Given a node in the tree, the routine gives leaf labels
// of the subtree rooted at 'node'.
void NewickTree_t::leafLabels(NewickNode_t *_node, set<string> &leaves)
{
  queue<NewickNode_t *> bfs;
  bfs.push(_node);
  NewickNode_t *node;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    int numChildren = node->children_m.size();
    if ( numChildren==0 ) // leaf
    {
      leaves.insert( node->label );
    }
    else
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);
      }
    }
  }
}


//--------------------------------------------- getSppProfs -----
void NewickTree_t::getSppProfs( vector<sppProf_t*> &sppProfs,
				vector<int> &nodeCut,
				int *annIdx,
				int minNA,
				map<int, string> &idxToAnn,
				map<string, string> &txParent)
/*
  Given nodeCut + some other parameters
  populate a vector of pointers to sppProf_t's
  with i-th element of sppProfs corresponding to the i-th
  node cut in nodeCut.
*/
{
  map<int, NewickNode_t*> idx2node;
  indexNewickNodes(idx2node);

  int n = nodeCut.size();
  for ( int i = 0; i < n; ++i )
  {
    map<string,int> annCount; // count of a given annotation string in the subtree of a given cut-node

    queue<NewickNode_t *> bfs;
    bfs.push(idx2node[nodeCut[i]]);
    NewickNode_t *node;
    int naCount = 0;

    while ( !bfs.empty() )
    {
      node = bfs.front();
      bfs.pop();

      int numChildren = node->children_m.size();
      if ( numChildren==0 ) // leaf
      {
	annCount[ idxToAnn[ annIdx[node->idx] ] ]++;
	if ( annIdx[node->idx] == -2) naCount++;
      }
      else
      {
	for (int i = 0; i < numChildren; i++)
	{
	  bfs.push(node->children_m[i]);
	}
      }
    }

    if ( naCount >= minNA )
    {
      map<string,int> sppFreq;
      node = idx2node[nodeCut[i]]; // current node
      NewickNode_t * ancNode;

      if ( annCount.size() > 1 ) // there is at least one annotated element; setting cltrTx to annotation with the max frequency
      {
	cltrSpp(node, annIdx, idxToAnn, sppFreq); // table of species frequencies in the subtree of 'node'
	ancNode = node;
      }
      else // there are no annotated sequences in the cluster; we are looking for the nearest ancestral node with known taxonomy
      {
	ancNode = node->parent_m; // looking at the parent node

	int nSpp = cltrSpp(ancNode, annIdx, idxToAnn, sppFreq); // table of species frequencies in the subtree of 'node'
	while ( nSpp==0 ) // if the parent node does not have any annotation sequences
	{
	  ancNode = ancNode->parent_m; // move up to its parent node
	  nSpp = cltrSpp(ancNode, annIdx, idxToAnn, sppFreq);
	}
	//ancNode = ancNode->parent_m;
      }

      sppProf_t *sppProf = new sppProf_t();
      sppProf->node    = node;
      sppProf->ancNode = ancNode;
      sppProf->sppFreq = sppFreq;

      sppProfs.push_back( sppProf );
    }
  }
}

#if 0
//--------------------------------------------- strMap2Set -----
set<string> NewickTree_t::strMap2Set(map<string,int> &sppFreq)
/*
  Extracts keys from map<string,int> and puts them in to a set<string>
*/
{
  set<string> sppSet;
  map<string,int>::iterator it;
  for (it = sppFreq.begin(); it != sppFreq.end(); ++it)
  {
    if ( (*it).first != "NA" )
    {
      sppSet.insert( (*it).first );
    }
  }

  return sppSet;
}
#endif

//--------------------------------------------- strMap2Str -----
string NewickTree_t::strMap2Str(map<string,int> &sppFreq)
/*
  Extracts keys from map<string,int>, sorts them and then concatenates them in to a string
*/
{
  vector<string> spp;
  map<string,int>::iterator it;
  for (it = sppFreq.begin(); it != sppFreq.end(); ++it)
  {
    if ( (*it).first != "NA" )
    {
      spp.push_back( (*it).first );
    }
  }

  sort (spp.begin(), spp.end());

  string sppStr = spp[0];
  int n = spp.size();
  for ( int i = 1; i < n; i++ )
    sppStr += string(".") + spp[i];

  return sppStr;
}


//--------------------------------------------- getAncRoots -----
void NewickTree_t::getAncRoots(vector<sppProf_t*> &sppProfs,
			       vector<sppProf_t*> &ancRootProfs) // ancestral root nodes
/*
  Finding root ancestral nodes, which are ancestral nodes of
  sppProfs, such that on the path from the ancestral node to the root of the tree
  there a no other ancestral nodes with the same species profile (excluding frequencies)
*/
{
  // creating a table of indices of all ancestral nodes except those which are root_m or whose parent is root_m
  // for the ease of identifying ancestral nodes on the path to the root_m
  map<int, int> ancNodes;
  int n = sppProfs.size();
  for ( int i = 0; i < n; ++i )
  {
    NewickNode_t *node = sppProfs[i]->ancNode;
    if ( node != root_m && node->parent_m != root_m )
    {
      //node = node->parent_m;
      //fprintf(stderr,"i: %d\tidx: %d\r",i,node->idx);
      ancNodes[ node->idx ] = i;
    }
  }

  // iterating over all ancestral nodes and checking which of them have no other
  // ancestral nodes on the path to root_m
  for ( int i = 0; i < n; ++i )
  {
    if ( sppProfs[i]->ancNode == sppProfs[i]->node // cut node with some reference sequences in its subtree
	 || sppProfs[i]->ancNode == root_m
	 || (sppProfs[i]->ancNode)->parent_m == root_m )
    {
      ancRootProfs.push_back( sppProfs[i] );
    }
    else
    {
      NewickNode_t *node = sppProfs[i]->ancNode;
      map<int,int>::iterator it;

      while ( node != root_m )
      {
	if ( (it = ancNodes.find( node->idx )) != ancNodes.end() && // walking to the root_m, another ancestral node was found
	     strMap2Str(sppProfs[i]->sppFreq) == strMap2Str( sppProfs[it->second]->sppFreq ) ) // it->second is the index of the ancestral node in sppProfs that was found on the path to root_m
	  break;
	node = node->parent_m;
      }

      if ( node == root_m ) // no ancestral nodes were found on the path to the root_m
	ancRootProfs.push_back( sppProfs[i] );
    }
  }
}

//--------------------------------------------- rmTxErrors -----
void NewickTree_t::rmTxErrors( vector<sppProf_t*> &ancRootProfs,
			       int depth,
			       int *annIdx,
			       map<int, string> &idxToAnn)

/*
  Changing taxonomic assignment in a situation when one species is contained in a
  clade of another species of the same genus.

  depth parameter controls how many ancestors we want to visit before making a decision.

  If at the depth ancestor (going up towards the root) species frequencies are
  1 (for the current node) and (depth-1) for other species of the same genus, then
  sppFreq of the current node is replaced by the one of the other species nodes.

  The modification are only applied to ancestral nodes with ancNode=node, because
  others will have homogeneous clade structure by construction.
*/
{
  int n = ancRootProfs.size();
  for ( int i = 0; i < n; ++i )
  {
    if ( ancRootProfs[i]->ancNode == ancRootProfs[i]->node )
    {
      NewickNode_t *node = ancRootProfs[i]->ancNode;
      for ( int j = 0; j < depth; j++ )
	node = node->parent_m;

      // compute sppFreq for node
      map<string,int> sppFreq;
      cltrSpp(node, annIdx, idxToAnn, sppFreq);

      // finding max count species in the original node
      // and comparing its frequency in sppFreq

      //if ( sppFreq[ (ancRootProfs[i]->node)->label ] = sppFreq[
    }
  }
}


//--------------------------------------------- printTx -----
void NewickTree_t::printTx( const char *outFile,
			    vector<sppProf_t*> &ancRootProfs)
/*
  Generating taxonomic assignment to all leaves of the tree.

  The leaves of the subtree of a root ancestral node all have the same taxonomy
  based on the majority vote at the species or genus level.
*/
{
  FILE *fh = fOpen(outFile,"w");

  int n = ancRootProfs.size();
  for ( int i = 0; i < n; ++i )
  {
    NewickNode_t *node = ancRootProfs[i]->ancNode;
    map<string,int> sppFreq = ancRootProfs[i]->sppFreq;
    string cltrTx; // taxonomy for all query reads of the cluster with no annotated sequences

    // determining taxonomy based on the majority vote
    int maxCount = 0; // maximal species count
    map<string,int>::iterator it;
    for (it = sppFreq.begin(); it != sppFreq.end(); ++it)
    {
      if ( (*it).first != "NA" && (*it).second > maxCount )
      {
	maxCount = (*it).second;
	cltrTx   = (*it).first;
      }
    }

    // propagating the taxonomy to all leaves of the subtree of the ancestral root node
    vector<string> leaves;
    leafLabels(node, leaves);
    int m = leaves.size();
    for ( int j = 0; j < m; j++ )
      fprintf(fh, "%s\t%s\n", leaves[j].c_str(), cltrTx.c_str());
  }

  fclose(fh);
}

//---------------------------------------------------------- updateLabels ----
void NewickTree_t::updateLabels( strSet_t &tx2seqIDs )
/*
  To test txSet2txTree() this routine updates tree labels so they
  include the number of sequences associated with each node
*/
{
  queue<NewickNode_t *> bfs;
  bfs.push(root_m);
  NewickNode_t *node;
  strSet_t::iterator it;
  char nStr[16];

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    if ( node != root_m )
    {
      it=tx2seqIDs.find( node->label );
      if ( it != tx2seqIDs.end() )
      {
	int n = (it->second).size();
	sprintf(nStr,"%d",n);
	node->label += string("{") + string(nStr) + string("}");
      }
      else
      {
	fprintf(stderr,"ERROR in %s at line %d: node with label %s not found in tx2seqIDs map\n",
		__FILE__,__LINE__, node->label.c_str());
	exit(1);
      }
    }

    int numChildren = node->children_m.size();
    if ( numChildren )
    {
      for (int i = 0; i < numChildren; i++)
	bfs.push(node->children_m[i]);
    }
  }
}


//---------------------------------------------------------- txSet2txTree ----
void NewickTree_t::txSet2txTree( strSet_t &tx2seqIDs )
/*
  tx2seqIDs maps tx (assumed to be labels of the leaves of the tree) into set of
  seqIDs corresponding to the given taxon, tx.

  The routine extends tx2seqIDs to map also internal nodes to seqIDs associated
  with the leaves of the corresponding subtree.
*/
{
  queue<NewickNode_t *> bfs;
  bfs.push(root_m);
  NewickNode_t *node;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    int numChildren = node->children_m.size();
    if ( numChildren )
    {
      for (int i = 0; i < numChildren; i++)
      {
	bfs.push(node->children_m[i]);

	set<string> seqIDs; // sequence IDs from all leaves of the subtree of node->children_m[i]
	set<string> leaves;
	leafLabels(node->children_m[i], leaves);
	set<string>::iterator it;
	for (it=leaves.begin(); it!=leaves.end(); it++)
	{
	  set<string> ids = tx2seqIDs[ *it ];
	  set_union( seqIDs.begin(), seqIDs.end(),
		     ids.begin(), ids.end(),
		     inserter(seqIDs, seqIDs.begin()) );
	}

	if ( (node->children_m[i])->label.empty() )
	{
	  cerr << "ERROR in " << __FILE__ << " at line " << __LINE__ << ": internal node is missing lable\n" << endl;
	}

	tx2seqIDs[ (node->children_m[i])->label ] = seqIDs;
      }
    }
  }
}


//---------------------------------------------------------- inodeTx ----
void NewickTree_t::inodeTx( const char *fullTxFile, map<string, string> &inodeTx )
/*
  Creating <interna node> => <taxonomy> table inodeTx
*/
{
  // parse fullTxFile that has the following structure

  // BVAB1	g_Shuttleworthia	f_Lachnospiraceae	o_Clostridiales	c_Clostridia	p_Firmicutes	d_Bacteria
  // BVAB2	g_Acetivibrio	f_Ruminococcaceae	o_Clostridiales	c_Clostridia	p_Firmicutes	d_Bacteria
  // BVAB3	g_Acetivibrio	f_Ruminococcaceae	o_Clostridiales	c_Clostridia	p_Firmicutes	d_Bacteria
  // Dialister_sp._type_1	g_Dialister	f_Veillonellaceae	o_Clostridiales	c_Clostridia	p_Firmicutes	d_Bacteria

  const int NUM_TX = 7;
  char ***txTbl;
  int nRows, nCols;
  readCharTbl( fullTxFile, &txTbl, &nRows, &nCols );

  #if 0
  fprintf(stderr,"fullTxFile txTbl:\n");
  printCharTbl(txTbl, 10, nCols); // test
  #endif

  map<string, vector<string> > fullTx; // fullTx[speciesName] = vector of the corresponding higher rank taxonomies
  // for example
  // BVAB1	g_Shuttleworthia	f_Lachnospiraceae	o_Clostridiales	c_Clostridia	p_Firmicutes	d_Bacteria
  // corresponds to
  // fullTx[BVAB1] = (g_Shuttleworthia, f_Lachnospiraceae, o_Clostridiales, c_Clostridia, p_Firmicutes, d_Bacteria)
  charTbl2strVect( txTbl, nRows, nCols, fullTx);

  map<string, vector<string> >::iterator it1;

  #if 0
  map<string, vector<string> >::iterator it1;
  for ( it1 = fullTx.begin(); it1 != fullTx.end(); it1++ )
  {
    fprintf(stderr, "%s: ", (it1->first).c_str());
    vector<string> v = it1->second;
    for ( int i = 0; i < (int)v.size(); i++ )
      fprintf(stderr, "%s ", v[i].c_str());
    fprintf(stderr, "\n");
  }

  for ( it1 = fullTx.begin(); it1 != fullTx.end(); it1++ )
  {
    vector<string> v = it1->second;
    fprintf(stderr, "%s: %d\n", (it1->first).c_str(), (int)v.size());
  }
  fprintf(stderr, "after\n");
  #endif

  // traverse the tree and for each internal node find consensus taxonomic rank
  // of all leaves of the corresponding subtree
  queue<NewickNode_t *> bfs;
  bfs.push(root_m);
  NewickNode_t *node;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    int numChildren = node->children_m.size();
    if ( numChildren ) // internal node
    {
      for (int i = 0; i < numChildren; i++)
	bfs.push(node->children_m[i]);

      vector<string> spp;
      leafLabels(node, spp);

      int txIdx;
      set<string> tx;
      set<string>::iterator it;
      for ( txIdx = 0; txIdx < NUM_TX; txIdx++ )
      {
	tx.erase( tx.begin(), tx.end() ); // erasing all elements of tx

	int n = spp.size();
	for ( int j = 0; j < n; j++ )
	{
	  it1 = fullTx.find( spp[j] );
	  if ( it1==fullTx.end() )
	    fprintf(stderr, "%s not found in fullTx\n", spp[j].c_str());
	  tx.insert( fullTx[ spp[j] ][txIdx] );
	}

	if ( tx.size() == 1 )
	  break;
      }

      it = tx.begin(); // now it points to the first and unique element of tx
      inodeTx[ node->label ] = *it;
    }
  }
}


//---------------------------------------------------------- modelIdx ----
// populate model_idx fields of all nodes of the tree
//
// when node labels correspond to model labels model_idx is the index of the node
// label in MarkovChains2_t's modelIds_m
void NewickTree_t::modelIdx( vector<string> &modelIds )
{
  map<string, int> modelIdx;
  int n = modelIds.size();
  for ( int i = 0; i < n; i++ )
    modelIdx[ modelIds[i] ] = i;

  queue<NewickNode_t *> bfs;
  bfs.push(root_m);
  NewickNode_t *node;

  map<string, int>::iterator it;

  while ( !bfs.empty() )
  {
    node = bfs.front();
    bfs.pop();

    if ( node != root_m )
    {
      it = modelIdx.find( node->label );
      if ( it != modelIdx.end() )
	node->model_idx = modelIdx[ node->label ];
      #if 0
      else
	fprintf(stderr, "ERROR in %s at line %d: Node label %s not found in modelIds\n",
		__FILE__, __LINE__, (node->label).c_str());
      #endif
    }

    int numChildren = node->children_m.size();
    if ( numChildren ) // internal node
    {
      for (int i = 0; i < numChildren; i++)
	bfs.push(node->children_m[i]);
    }
  }
}


// ----------------------------------------- readNewickNode() ---------------------------------
//
// This code was lifted/adopted from the spimap-1.1 package
// ~/devel/MCclassifier/packages/spimap-1.1/src
// http://compbio.mit.edu/spimap/
// https://academic.oup.com/mbe/article-lookup/doi/10.1093/molbev/msq189
//

float readDist( FILE *infile )
{
    float dist = 0;
    fscanf(infile, "%f", &dist);
    return dist;
}


char readChar(FILE *stream, int &depth)
{
    char chr;
    do {
        if (fread(&chr, sizeof(char), 1, stream) != 1) {
            // indicate EOF
            return '\0';
        }
    } while (chr == ' ' || chr == '\n');

    // keep track of paren depth
    if (chr == '(') depth++;
    if (chr == ')') depth--;

    return chr;
}


char readUntil(FILE *stream, string &token, const char *stops, int &depth)
{
    char chr;
    token = "";
    while (true) {
        chr = readChar(stream, depth);
        if (!chr)
            return chr;

        // compare char to stop characters
        for (const char *i=stops; *i; i++) {
            if (chr == *i)
                return chr;
        }
        token += chr;
    }
}

string getIdent( int depth )
{
   string identStr = "";
   for ( int i = 0; i < depth; i++ )
   {
	  identStr += string("  ");
   }

   return identStr;
}

NewickNode_t *readNewickNode(FILE *infile, NewickTree_t *tree, NewickNode_t *parent, int &depth)
{
  #define DEBUG_RNN 0
  #if DEBUG_RNN
   string ident = getIdent( depth );
   fprintf(stderr, "%sIn readNewickNode() depth: %d\n", ident.c_str(), depth);
  #endif

   char chr, char1;
   NewickNode_t *node;
   string token;

   // read first character
   if ( !(char1  = readChar(infile, depth)) )
   {
	  fprintf(stderr, "\n\n\tERROR: Unexpected end of file\n\n");
	  exit(1);
   }

   if ( char1 == '(' ) // read internal node
   {
	  int depth2 = depth;
	  if ( parent )
	  {
		 node = parent->addChild();
		 #if DEBUG_RNN
		 if ( depth )
		 {
			string ident = getIdent( depth );
			fprintf( stderr, "%sCurrently at depth: %d - Adding a child to a parent\n", ident.c_str(), depth );
		 }
		 else
		 {
			fprintf( stderr, "Currently at depth: %d - Adding a child to a parent\n", depth );
		 }
		 #endif
	  }
	  else
	  {
		 node = new NewickNode_t();
		 tree->setRoot( node );
		 #if DEBUG_RNN
		 fprintf( stderr, "Currently at depth: %d - Creating a root\n", depth );
         #endif
	  }


	  // updating node's depth
	  node->depth_m = depth - 1;

	  // read all child nodes at this depth
	  while ( depth == depth2 )
	  {
		 NewickNode_t *child = readNewickNode(infile, tree, node, depth);
		 if (!child)
		   return NULL;
	  }

	  // Assigning a numeric index to the node
	  node->idx = tree->getMinIdx();
	  tree->decrementMinIdx();

	  // read branch_length for this node
	  // this fragment assumes that the internal nodes do not have labels
	  char chr = readUntil(infile, token, "):,;", depth);

	  if ( !token.empty() ) // Reading node's label
	  {
		 node->label = trim(token.c_str());

         #if DEBUG_RNN
		 string ident = getIdent( depth );
		 fprintf( stderr, "%sNode name: %s\n", ident.c_str(), token.c_str() );
         #endif
	  }

	  if (chr == ':')
	  {
		 node->branch_length = readDist( infile );

		 #if DEBUG_RNN
		 string ident = getIdent( depth );
		 fprintf( stderr, "%snode->branch_length: %f\n", ident.c_str(), node->branch_length );
		 #endif

		 if ( !(chr = readUntil(infile, token, "):,;", depth)) )
		   return NULL;
	  }
	  //if (chr == ';' && depth == 0)
	  //    return node;

	  #if DEBUG_RNN
	  string ident = getIdent( node->depth_m );
	  fprintf( stderr, "%sNode's depth: %d\n", ident.c_str(), node->depth_m );
	  #endif

	  return node;
   }
   else
   {
	  // Parsing leaf

	  if (parent)
	  {
		node = parent->addChild();
		#if DEBUG_RNN
		string ident = getIdent(depth);
		fprintf( stderr, "%sCurrently at depth: %d - Found leaf: adding it as a child to a parent\n", ident.c_str(), depth );
		#endif
	  }
	  else
	  {
		 fprintf( stderr, "\n\n\tERROR: Found leaf without a parent\n\n" );
		 exit(1);
	  }

	  // Assigning to the leaf a numeric index
	  node->idx = tree->leafCount();
	  tree->incrementLeafCount();

	  // updating leaf's depth
	  node->depth_m = depth;

	  // Reading leaf label
	  if ( !(chr = readUntil(infile, token, ":),", depth)) )
		return NULL;
	  token = char1 + trim(token.c_str());
	  node->label = token;

	  #if DEBUG_RNN
	  string ident = getIdent( depth );
	  fprintf( stderr, "%sLeaf name: %s\tdepth: %d\tidx: %d\n", ident.c_str(), token.c_str(), node->depth_m, node->idx );
      #endif

	  // read distance for this node
	  if (chr == ':')
	  {
		 node->branch_length = readDist( infile );

		 #if DEBUG_RNN
		 fprintf( stderr, "%sLeaf's branch length: %f\n", ident.c_str(), node->branch_length );
         #endif

		 if ( !(chr = readUntil( infile, token, ":),", depth )) )
		   return NULL;
	  }

	  return node;
   }
}


NewickTree_t *readNewickTree( const char *filename )
{
   FILE *infile = fOpen(filename, "r");

   int depth = 0;
   NewickTree_t *tree = new NewickTree_t();
   NewickNode_t * node = readNewickNode(infile, tree, NULL, depth);
   tree->setRoot( node );

   fclose(infile);

   return tree;
}
