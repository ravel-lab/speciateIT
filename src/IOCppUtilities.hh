#ifndef IOUTILITIES_HH
#define IOUTILITIES_HH

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


#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>

//#include "inputPar.hh"

using namespace std;

// -- vicut stuff

int * parseAnnotation(const char *inFile, char **rowIds, int nIds,
                      map<string,int> &annToIdx, bool verbose);

void parseQueryIDs(const char *inFile, map<string, bool> &queryIDs);

void queryIDsToAnnIdx ( map<string, bool> &queryIDs,
			char **leafIds, int nIds,
			map<string,int> &annToIdx,
			int *annIdx);

void read2cols(const char *inFile, map<string,string> &tbl);

// ---

string trim( const char *word );

void printStringSet(set<string> &S);

void childrenMap( char ***tbl, int nRows, int nCols, map<string, set<string> > &children);
void txRankTbl( char ***tbl, int nRows, int nCols, vector< set<string> > &txRank);

typedef map<string, set<string> > strSet_t;
void txTbl2txSet( char ***tbl, int nRows, strSet_t &tx2seqIDs);
void charTbl2strVect( char ***tbl, int nRows, int nCols, map<string, vector<string> > &id2rest);

void writeStrVector(const char *outFile, vector<string> &v, const char *sep="\n");

/// returns number of sequences in a fasta file
int numRecordsInFasta( const char *file );

/// read 2 column of strings table
void read_2s_tbl( const char *inFile, map<string, string> &sMap );

/// read matrix of numerical data, set header to 1 if input file has a header
void readMatrix( const char *inputFile, double ***matrix, int *nrow, int *ncol, int header );

/// read lines of the input file
void readLines( const char *inFile, vector<char*> &lines );
void readLines( const char *inFile, vector<string> &lines );
void readLines( const char *inFile, char ***list, int *nRows, int header );
//void readLines( inPar_t *p );

/// read a fasta file
void readFasta( const char *file, map<string,string> &seqTbl);

/// read one record of fasta file at a time
bool getNextFastaRecord( FILE *fp, char *&id, char *&seq, int &seqLen);
bool getNextFastaRecord( FILE *fp, char *&id, char *data, size_t alloc, char *seq, int &seqLen);

/// read first sequence record of fasta file
bool readFastaRecord(const char *file, char *&id,
		     char *&header, char *&seq,
		     int &seqLen);

/// parses comma separated list of ints
void parseCommaList( char *listStr, vector<int> &list );

/// prints to stderr list of strings with header 'name'
void printCharList( char **list, int len, const char *name );

/// write fragStartPos\tprobs to out file handle
void writeProbs( FILE *out,
		 int fragStartPos,
		 double *probs, int nProbs);

void writeProbs( FILE *out,
		 const char *seqId,
		 double *probs, int nProbs);

void writeClassification( FILE *out,
			  const char *seqId,
			  vector< pair<string, double> > &score,
			  map<string, vector<string> > &fullTx);

void writeALogOdds( FILE *out,
		    const char *seqId,
		    double *x, int n,
		    int m,
		    vector<string> &path,
		    string &nodeLabel);

/// write header of the fragment probability table
void writeHeader( FILE *out, vector<char *> &chromoIds );

/// write vector of doubles to file
void writeDVector( vector<double> &m, const char *yFile );

/// write a matrix to a file
void writeMatrix( double **m, int nRows, int nCols, const char *file );

// print table of doubles to stdout
void printDblTbl(double **tbl, int nRows, int nCols);

/// print vector of class T elements to stdout
template<typename T>
void printVector(vector<T> &v, const char *sep="\t", int width = 0)
{
  cout.precision(10);
  if ( width > 0 )
    cout.width(width);

  int n = v.size();
  for ( int i = 0; i < n; ++i )
    cout << v[i] << sep;
  cout << endl;

  if ( width > 0 )
    cout.width(1);
}

template<typename T>
void printVector(const char *outFile, vector<T> &v, const char *header=NULL)
{
  ofstream fout( outFile, ios::out );
  fout.precision(15);

  if ( header )
    fout << header << endl;

  int n = v.size();
  for ( int i = 0; i < n; ++i )
    fout << v[i] << endl;
  fout << endl;
  fout.close();
}

/// print vector of class T elements to stdout
template<typename T>
void printArray(T *a, int n)
{
  cout.precision(10);
  for ( int i = 0; i < n; ++i )
    cout << a[i] << " ";
  cout << endl;
}

/// check if a class, T, instant is a member of a vector<T>
template<typename T>
int exists_in_vector( vector<T> &v, T a )
{
  typename vector<T>::iterator it = find ( v.begin(), v.end(), a );

  int ret = 1;
  if ( it == v.end() )
    ret = 0;

  return ret;
}

/// check if a class, T, instant is a member of a map<T, S>
template<typename T, typename S>
int exists_in_map( map<T, S> &m, T a )
{
  typename map<T, S>::iterator it = m.find( a );

  int ret = 1;
  if ( it == m.end() )
    ret = 0;

  return ret;
}


#endif
