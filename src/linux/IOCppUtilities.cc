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

#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "IOCppUtilities.hh"
#include "IOCUtilities.h"
#include "CUtilities.h"
#include "CppUtilities.hh"
#include "strings.hh"

#define LINE_LEN  10000

// typedef map<string, set<string> > strSet_t;

//---------------------------------------------------------- printStringSet ----
void printStringSet(set<string> &S)
{
  set<string>::iterator it;
  for (it=S.begin(); it!=S.end(); it++)
    cout << *it << " ";
  cout << endl;
}


//---------------------------------------------------------- txTbl2txSet ----
void txTbl2txSet( char ***tbl, int nRows, strSet_t &tx2seqIDs)
/*
  It is assumed that tbl has two columns: seqID and taxon.

  tx2seqIDs maps tx into set of seqIDs corresponding to the given taxon, tx.
*/
{
  for ( int i = 0; i < nRows; ++i )
  {
    tx2seqIDs[ string(tbl[i][1]) ].insert( string(tbl[i][0]) );
  }
}

//---------------------------------------------------------- txRankTbl ----
// creating taxonomic rank index (species = 0, genus = 1, etc ) => set of
// corresponding taxonomic ranks in fullTx
void txRankTbl( char ***tbl, int nRows, int nCols, vector< set<string> > &txRank)
{
  int nCols1 = nCols - 1;
  for ( int i = 0; i < nCols1; ++i )
  {
    set<string> s;
    for ( int j = 0; j < nRows; ++j )
    {
      //fprintf(stderr, "i=%d j=%d inserting tbl[j][i]=%s\n", i, j, tbl[j][i]);
      s.insert( string(tbl[j][i]) );
    }
    txRank.push_back(s);
  }
}

//---------------------------------------------------------- childrenMap ----
// creating children mapping
// children[ d_Bacteria ] = set( all phyla in fullTx table that are children of d_Bacteria )
void childrenMap( char ***tbl, int nRows, int nCols, map<string, set<string> > &children)
{
  // order of taxonomic ranks in tbl
  // BVAB1	g_Shuttleworthia	f_Lachnospiraceae	o_Clostridiales	c_Clostridia	p_Firmicutes	d_Bacteria

  int nCols1 = nCols - 1;
  for ( int i = 0; i < nRows; ++i )
    for ( int j = nCols1; j > 0; j-- )
    {
      children[ string(tbl[i][j]) ].insert( string(tbl[i][j-1]) );

      if ( string(tbl[i][j]) == "_" )
      {
	fprintf(stderr, "parent=_\ni=%d\n\n",i);
	for ( int j = 0; j < nCols; j++ )
	  fprintf(stderr,"%s  ",tbl[i][j]);
	  fprintf(stderr,"\n\n");
      }
    }
}


//---------------------------------------------------------- charTbl2strVect ----
void charTbl2strVect( char ***tbl, int nRows, int nCols, map<string, vector<string> > &id2rest)
/*
  tbl has at least two columns

  id2rest maps element of the first column to the vector of the remaining ones
*/
{
  for ( int i = 0; i < nRows; ++i )
    for ( int j = 1; j < nCols; ++j )
      id2rest[ string(tbl[i][0]) ].push_back( string(tbl[i][j]) );
}

// ------------------------------ writeStrVector ------------------------
void writeStrVector(const char *outFile, vector<string> &v, const char *sep)
{
  FILE *file = fOpen(outFile,"w");
  int n = v.size();
  for ( int i = 0; i < n; ++i )
    fprintf(file, "%s%s",v[i].c_str(), sep);
  fclose(file);
}

//--------------------------------------------------- numRecordsInFasta ----
/// number of sequences in a fasta file
int numRecordsInFasta( const char *file )
{
  int nRecs;
  char cmd[1024];

  char *tFile = (char*)malloc(2*strlen(file)*sizeof(char));
  tFile[0] = '\0';
  strcat(tFile, file);
  strcat(tFile, ".size");

  sprintf(cmd,"grep '>' %s | wc -l > %s",file,tFile);
  system(cmd);

  FILE *in = fOpen(tFile, "r");
  fscanf(in,"%d",&nRecs);
  fclose(in);

  //cerr << "cmd=" << cmd << endl;
  sprintf(cmd,"rm -f %s",tFile);
  system(cmd);
  free(tFile);

  return nRecs;
}


//-------------------------------------------------- getNextFastaRecord ----
/// Memory allocation is moved ouside of the routine
///
/// Reads a record from a fasta file. Returns true on
/// success, otherwise false.
///
/// id     - sequence identifier
/// seq    - sequence string
/// seqLen - sequence length
/// header - fasta header
///
/// Addapted from Adam Phyllipy's Sequence_t class
///
bool getNextFastaRecord( FILE *fp, char *&id, char *data, size_t alloc, char *seq, int &seqLen)
{
  int ch;

  //-- Find the beginning of the fasta record
  while ( (ch = getc(fp)) != '>' )
    if ( ch == EOF )
      return false;

  //-- Get the fasta header
  if ( !ReadLine(fp, data, alloc) )
    return false;

  //STRDUP(id, strlopspace(data));
  id = strlopspace(data);

  //-- Get the sequence data
  seqLen = 0;
  while ( (ch = getc(fp)) != '>' && ch != EOF )
  {
    if ( isspace(ch) ) continue;

    if ( size_t(seqLen + 1) >= alloc )
    {
      cerr << "Error: " << __FILE__
	   << " in getNextFastaRecord() at line " << __LINE__
	   << "\tseqLen+1>alloc" << endl;
      return false;
      exit(1);
    }

    seq[seqLen++] = toupper(ch);
  }

  seq[seqLen] = '\0';

  if ( ch == '>' ) ungetc(ch, fp);

  return true;
}

//-------------------------------------------------- getNextFastaRecord ----
/// Reads a record from a fasta file. Returns true on
/// success, otherwise false.
///
/// id     - sequence identifier
/// seq    - sequence string
/// seqLen - sequence length
/// header - fasta header
///
/// Addapted from Adam Phyllipy's Sequence_t class
///
bool getNextFastaRecord( FILE *fp, char *&id, char *&seq, int &seqLen)
{
  int ch;

  //-- Find the beginning of the fasta record
  while ( (ch = getc(fp)) != '>' )
    if ( ch == EOF )
      return false;

  //-- Reserve default amount
  size_t alloc = 1*1024*1024;
  char *data;
  MALLOC(data, char*, alloc * sizeof(char));
  MALLOC(seq, char*, alloc * sizeof(char));

  //-- Get the fasta header
  if ( !ReadLine(fp, data, alloc) )
    return false;

  STRDUP(id, strlopspace(data));

  //-- Get the sequence data
  seqLen = 0;
  while ( (ch = getc(fp)) != '>' && ch != EOF )
  {
    if ( isspace(ch) ) continue;

    if ( size_t(seqLen + 1) >= alloc )
    {
      alloc = (seqLen + 1) * 2;
      REALLOC(seq, char*, alloc * sizeof(char));
    }

    seq[seqLen++] = toupper(ch);
  }

  seq[seqLen] = '\0';
  REALLOC(seq, char*, (seqLen+1) * sizeof(char));

  if ( ch == '>' ) ungetc(ch, fp);

  free(data);

  return true;
}


//----------------------------------------------- readFasta ----------
void readFasta( const char *file, map<string,string> &seqTbl)
/*
  for each record of a fasta file with sequence ID, id, and sequence, seq, create
  a pair seqTbl.insert(make_pair(string(id),string(seq)))
*/
{
  char *id, *seq;
  int seqLen;

  FILE *in = fOpen(file, "r");

  while ( getNextFastaRecord( in, id, seq, seqLen) )
  {
    seqTbl.insert(make_pair(string(id),string(seq)));
    free(id);
    free(seq);
  }

  fclose(in);
}



//----------------------------------------------- readFastaRecord ----
/// Reads a single sequence record from a fasta file. Returns true on
/// success, otherwise false.
///
/// id     - sequence identifier
/// seq    - sequence string
/// seqLen - sequence length
/// header - fasta header
///
/// Addapted from Adam Phyllipy's Sequence_t class readFastaRecord()
///
bool readFastaRecord( const char *file, char *&id,
	        char *&header, char *&seq, int &seqLen)
{
  FILE *fp = fOpen(file, "r");
  int ch;

  //-- Find the beginning of the fasta record
  while ( (ch = getc(fp)) != '>' )
    if ( ch == EOF )
      return false;

  //-- Reserve default amount
  size_t alloc = 8*1024*1024;
  char *data;
  MALLOC(data, char*, alloc * sizeof(char));
  MALLOC(seq, char*, alloc * sizeof(char));

  //-- Get the fasta header
  if ( !ReadLine(fp, data, alloc) )
    return false;

  STRDUP(header, chomp(data));
  STRDUP(id, strlopspace(data));

  //-- Get the sequence data
  seqLen = 0;
  while ( (ch = getc(fp)) != '>' && ch != EOF )
  {
    if ( isspace(ch) ) continue;

    if ( size_t(seqLen + 1) >= alloc )
    {
      alloc = (seqLen + 1) * 2;
      REALLOC(seq, char*, alloc * sizeof(char));
    }

    seq[seqLen++] = toupper(ch);
  }

  seq[seqLen] = '\0';

  if ( ch == '>' ) ungetc(ch, fp);

  fclose(fp);
  free(data);

  return true;
}

// -------------------------------- writeProbs ------------------------
// write fragStartPos\tprobs to out file handle
void writeProbs( FILE *out,
		 int fragStartPos,
		 double *probs, int nProbs)
{
  fprintf(out,"%d\t",fragStartPos);

  int nProbs1 = nProbs - 1;
  int i = 0;
  for ( ; i < nProbs1; ++i )
    fprintf(out,"%f\t", probs[i]);
  fprintf(out,"%f\n", probs[i]);
}


// -------------------------------- writeProbs ------------------------
// write seqId\tprobs to out file handle
void writeProbs( FILE *out,
		 const char *seqId,
		 double *probs, int nProbs)
{
  fprintf(out,"%s\t",seqId);

  int nProbs1 = nProbs - 1;
  int i = 0;
  for ( ; i < nProbs1; ++i )
    fprintf(out,"%f\t", probs[i]);
  fprintf(out,"%f\n", probs[i]);
}

// -------------------------------- writeClassification ------------------------
// print scores of a given seqId
//
// If the reference tree is missing some taxonomic ranks (because it would have
// only one child), the missing rank inherits the score of the higher rank.
void writeClassification( FILE *out,
			  const char *seqId,
			  vector< pair<string, double> > &score,
			  map<string, vector<string> > &fullTx)
{
  #if 0
  fprintf(out,"%s",seqId);
  int n = score.size();
  for ( int i = 0; i < n; ++i )
    fprintf(out,"\t%s\t%.1f", score[i].first.c_str(), score[i].second);
  fprintf(out,"\n");
  #endif

  int n = score.size();
  map<string, vector<string> >::iterator it = fullTx.find( score[n-1].first );
  if ( it == fullTx.end() )
  {
    fprintf(stderr, "ERROR in %s at line %d: Could not find %s in fullTx table\n",
	    __FILE__, __LINE__, score[n-1].first.c_str());
  }

  map<string, double> scoreMap;
  for ( int i = 0; i < n; ++i )
    scoreMap[ score[i].first ] = score[i].second;

  map<string, double>::iterator it2 = scoreMap.find( score[0].first );
  if ( it2 == scoreMap.end() )
  {
    fprintf(stderr, "ERROR in %s at line %d: Could not find %s in scoreMap table\n",
	    __FILE__, __LINE__, score[0].first.c_str());
    exit(1);
  }

  vector<string> tx = fullTx[ score[n-1].first ];
  tx.pop_back(); // remove d_Bacteria;
  reverse(tx.begin(), tx.end());
  tx.push_back( score[n-1].first );

  #if 0
  fprintf(out,"score[n-1].first: %s\n\n", score[n-1].first.c_str());
  n = tx.size();
  fprintf(out,"tx: ");
  for ( int i = 0; i < n; ++i )
    fprintf(out,"  %s", tx[i].c_str());
  fprintf(out,"\n");
  #endif

  fprintf(out,"%s",seqId);
  fprintf(out,"\t%s\t%.1f", score[0].first.c_str(), score[0].second);
  n = tx.size();
  for ( int i = 1; i < n; ++i )
  {
    it2 = scoreMap.find( tx[i] );
    if ( it2 == scoreMap.end() )
      scoreMap[ tx[i] ] = scoreMap[ tx[i-1] ];
    fprintf(out,"\t%s\t%.1f", tx[i].c_str(), scoreMap[ tx[i] ]);
  }
  fprintf(out,"\n");
}


// -------------------------------- writeALogOdds ------------------------
void writeALogOdds( FILE *out,
		    const char *seqId,
		    double *x, int n,
		    int m,
		    vector<string> &path,
		    string &nodeLabel)
{
  fprintf(out,"%s\t%s\t",seqId, nodeLabel.c_str());

  int i = 0;
  for ( ; i < m; ++i )
    fprintf(out,"%.1f\t", x[i]);
  //fprintf(out,"%.1f:%s\t", x[i], path[i].c_str());
  //fprintf(out,"\n" );
  #if 1
  int n1 = n - 1;
  for ( ; i < n1; ++i )
    fprintf(out,"%.1f\t", x[i]);
  if ( i < n1 )
    fprintf(out,"%.1f", x[i] );
  fprintf(out,"\n");
  #endif
}

// -------------------------------- writeALogOdds ------------------------
void writeALogOdds( FILE *out,
		    const char *seqId,
		    double x,
		    string &nodeLabel)
{
  fprintf(out,"%s\t%s\t%.1f\n",seqId, nodeLabel.c_str(), x);
}

// -------------------------------- writeHeader ------------------------
/// write header of the fragment probability table
void writeHeader( FILE *out, vector<char *> &ids )
{
  if ( !ids.size() )
  {
    cerr << "Error: writeHeader() in " << __FILE__
	 << " at line " << __LINE__
	 << " ids.size()=0"  << endl;
    return;
  }

  fprintf(out,"seqId\t");

  int n = ids.size() - 1;
  int i = 0;

  for ( ; i < n; ++i )
    fprintf(out,"%s\t", ids[i]);
  fprintf(out,"%s\n", ids[i]);
}


//------------------------------------------------------ printCharList ----
void printCharList( char **list, int len, const char *name )
{
  cerr << name << ": ";
  for ( int i = 0; i < len; ++i )
    cerr << list[i] << " ";
  cerr << endl;
}


//------------------------------------------------------ parseCommaList ----
void parseCommaList( char *listStr, vector<int> &list )
{
  list.clear();

  char *digitStr, *brkt;

  for (digitStr = strtok_r(listStr, ":", &brkt);
       digitStr;
       digitStr = strtok_r(NULL, ":", &brkt))
    list.push_back(atoi(digitStr));
}


//------------------------------------------------------ readMatrix ----
/*
  Read the whole input file into a buffer and then process it line by line
  Input file is expected to have the same number of fields (tab delimited) in each row
  matrix[i] holds the content of the i-th row of the input file.
  set header=1 if the input file has a header, otherwise set it to 0
  nrow is the number of lines of the input file (excluding header line)
  ncol is the number of columns
*/
void readMatrix( const char *inputFile, double ***matrix, int *nrow, int *ncol, int header )
{
  char line[LINE_LEN]; // array for holding a line
  int lineLen = 0;     // length of the line read from the buffer
  int offset = 0;      // file position offset
  int count = 0;       // line counter
  char *ptr = NULL;    // auxiliary pointer used with strtod function
  char * field;
  int nread;           // number of bytes read into the buffer

  size_t bufferSize = fileSize(inputFile);
  char * buffer     = (char *) malloc( bufferSize * sizeof(char) );
  double ** data    = (double **)malloc( bufferSize * sizeof(double *) );
  int in            = open(inputFile, O_RDONLY);

  if ( (nread = read( in, buffer, bufferSize-1 )) )
  {
    lineLen = 0;
    offset = 0;
    buffer[nread] = '\0';

    // check the number of columns
    lineLen = readLine( offset, buffer, bufferSize, line, LINE_LEN );
    offset += lineLen;
    field = strtok ( line, " \t" );
    *ncol = 1;

    while( field )
    {
      field = strtok ( NULL, " \t\n" );
      (*ncol)++;
    }
    (*ncol)--;

    if ( !header )
      offset = 0;

    while ( offset < nread )
    {
      lineLen = readLine( offset, buffer, bufferSize, line, LINE_LEN );
      offset += lineLen;

      data[count] = (double *) malloc ( (*ncol) * sizeof ( double ) );
      field = strtok ( line, " \t" );
      data[count][0] = strtod ( field, &ptr );

      int col = 1;
      while ( col < *ncol )
      {
	field = strtok ( NULL, " \t\n" );
	data[count][col++] = strtod ( field, &ptr );
      }
      count++;
    }
  }

  close(in);
  free ( buffer );

  data = (double **)realloc(data, count * sizeof(double *));
  *matrix = data;
  *nrow = count;
}



/*
  Reads the rows of the input file, inFile, into lines
  lines[i] is the i-th row of the input file (it is assumed there is no header)
*/
void readLines( const char *inFile, vector<char*> &lines )
{
  lines.clear();
  size_t bufferSize = fileSize(inFile);
  int in            = open(inFile, O_RDONLY);

  char *buffer;
  MALLOC(buffer, char* , bufferSize * sizeof(char) );

  char line[LINE_LEN]; // array for holding a line
  int nread;           // number of bytes read into the buffer
  char *data;

  if ( (nread = read ( in, buffer, bufferSize-1 )) )
  {
    int lineLen   = 0;  // length of the line read from the buffer
    int offset    = 0;  // file position offset
    buffer[nread] = '\0';

    while ( offset < nread )
    {
      lineLen = readLine( offset, buffer, bufferSize, line, LINE_LEN );
      offset += lineLen;
      line[lineLen] = '\0';

      STRDUP(data,line);
      lines.push_back(data);
    }
  }

  close(in);
  free(buffer);
}

/*
  Reads the rows of the input file into a list
  list[i] is the i-th row of the input file (excluding header line)
  set header=1 if the input file has a header, otherwise set it to 0
  nRows is the number of lines of the input file (excluding header line)
*/
void readLines( const char *inputFile, char ***list, int *nRows, int header )
{
  char line[LINE_LEN]; // array for holding a line
  int lineLen = 0;     // length of the line read from the buffer
  int offset = 0;      // file position offset
  int count = 0;       // line counter
  int nread;           // number of bytes read into the buffer

  size_t bufferSize = fileSize(inputFile);
  char * buffer     = (char *) malloc( bufferSize * sizeof(char) );
  char ** data      = (char **)malloc( bufferSize * sizeof(char *) );
  int in            = open(inputFile, O_RDONLY);

  if ( (nread = read ( in, buffer, bufferSize-1 )) )
  {
    lineLen = 0;
    offset = 0;
    buffer[nread] = '\0';

    if ( header )
    {
      lineLen = readLine( offset, buffer, bufferSize, line, LINE_LEN );
      offset += lineLen;
    }

    while ( offset < nread )
    {
      lineLen = readLine( offset, buffer, bufferSize, line, LINE_LEN );
      offset += lineLen;
      line[lineLen] = '\0';
      STRDUP(data[count],line);
      count++;
    }
  }

  close(in);
  free(buffer);

  data = (char **)realloc(data, count * sizeof(char *));
  *list = data;
  *nRows = count;
}

//---------------------------------------------------------- writeMatrix ----
void writeMatrix( double **m, int nRows, int nCols, const char *file )
{
  FILE *out = fOpen(file,"w");
  int nCols1 = nCols-1;
  int i, j;

  for ( i = 0; i < nRows; ++i )
  {
    for ( j = 0; j < nCols1; ++j )
      fprintf(out,"%lf\t",m[i][j]);
    fprintf(out,"%lf\n",m[i][j]);
  }

  fclose(out);
}

//---------------------------------------------------------- printDblTbl ----
void printDblTbl(double **tbl, int nRows, int nCols)
{
  for ( int i = 0; i < nRows; ++i )
  {
    for ( int j = 0; j < nCols; ++j )
      cout << tbl[i][j] << " ";

    cout << endl;
  }
  cout << endl;
}
