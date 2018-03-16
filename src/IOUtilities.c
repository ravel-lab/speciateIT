/*
 * C Library of Input/Output routines
 * Copyright (C) 2010 Pawel Gajer pgajer@gmail.com
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 *
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

#include "IOUtilities.h"
#include "CUtilities.h"

//----------------------------------------------------------------- fOpen ----
//! open a file handle and throw an error message if it cannot be opened
FILE *_fOpen ( const char *file, const char *format, const char * cppfile, int line )
{
    FILE *f;

    errno = 0;

    if ( file == NULL )
        file = "\0"; // fopen() may choke on a NULL file

    if ( ( f = fopen ( file, format ) ) == NULL )
    {
        if ( line == -1 )
            fprintf ( stderr, "Cannot open %s for %s: %s\n\n", file, format, strerror ( errno ) );
        else
            fprintf ( stderr, "%s line:%d  Cannot open %s for %s: %s\n\n",
                      cppfile, line, file, format, strerror ( errno ) );

        exit ( 1 );
    }
    return f;
}

char* GetLine(FILE* inputfile)
/*
  The function GetLine reads one line from the inputfile, and returns it as a
  null-terminated string. If inputfile is at EOF, a null pointer is returned.
  The calling routine should free the char* returned by GetLine.
*/
{ int c;
  int n = 0;
  int size = 1023;
  char* line = malloc((size+1)*sizeof(char));
  while ((c = getc(inputfile))!=EOF && c!='\r' && c!='\n')
  { if (n == size)
    { size *= 2;
      line = realloc(line,(size+1)*sizeof(char));
    }
    line[n] = (char)c;
    n++;
  }
  if (c=='\r')
  { c = getc(inputfile);
    if (c!='\n' && c!=EOF) ungetc(c,inputfile);
  }
  if (n==0 && c==EOF)
  { free(line);
    return 0;
  }
  line[n] = '\0';
  line = realloc(line,(n+1)*sizeof(char));
  return line;
}

static char* tokenize(char* s)
{
  char* p = s;
  while (1)
  {
    if (*p=='\0') return NULL;
    if (*p=='\t')
    {
      *p = '\0';
      return p+1;
    }
    p++;
  }
  /* Never get here */
  return NULL;
}


char * readTable( const char *inFile, double ***matrix, int *nRows, int *nCols,
		  char ***rowNames, char ***colNames )
/*
  Read a table with row and column names
  Input format
  PtId  \tcol1name\tcol2name\t ... \tcolNname
  row1Id\tx_1     \tx_2     \t ... \tx_dim

  Alternatively header may not have RowId lable
  col1name\tcol2name\t ... \tcolNname
  row1Id\tx_1     \tx_2     \t ... \tx_dim

  Reads the whole input file into a buffer and then process it line by line

  readTable() returns 0 on success and error string on error.

  All error messages are allocated with malloc, even if not strictly
  necessary. The reason is that any error message then can (and should) be
  safely freed.
*/
{
  int row, column;           /* Counters for data matrix */
  int fileRow, fileColumn;   /* Counters for rows and columns in the file */
  int n;
  int nFileColumns;
  char* line;
  char* s;
  char* token;

  FILE *file = fOpen(inFile,"r");

  // Parse header line (first line) to find out what the columns are

  line=GetLine(file);
  if(!line)
  {
    const char text[] = "Error: Attempt to read empty file";
    const int m = strlen(text) + 1;
    char* error = malloc(m*sizeof(char));
    strcpy(error,text);
    return error;
  }

  while(line[0]=='\0') /* Ignore completely empty lines */
  {
    free(line);
    line = GetLine(file);
    if(!line)
    {
      const char text[] = "Error: Failed to find first line in file";
      const int m = strlen(text) + 1;
      char* error = malloc(m*sizeof(char));
      strcpy(error,text);
      return error;
    }
  }

  s = tokenize(line);
  fileColumn = 1;

  while (s)
  {
    s = tokenize(s);
    fileColumn++;
  }
  free(line);

  nFileColumns = fileColumn;


  // check the number of elements in the second row
  // it can be either = nFileColumns
  // or
  // = nFileColumns+1
  line=GetLine(file);
  s = tokenize(line);
  fileColumn = 1;
  fileRow = 1;

  while (s)
  {
    s = tokenize(s);
    fileColumn++;
  }
  free(line);

  if (fileColumn < nFileColumns)
  {
    int n = 1024;
    char* text = malloc(n*sizeof(char));
    sprintf (text,
	     "Error reading line %d: only %d columns available (%d needed)",
	     fileRow, fileColumn, nFileColumns);
    n = strlen(text) + 1;
    text = realloc(text,n*sizeof(char));
    return text;
  }

  if (fileColumn > nFileColumns+1)
  {
    int n = 1024;
    char* text = malloc(n*sizeof(char));
    sprintf (text,
	     "Error reading line %d: %d columns given (%d needed)",
	     fileRow, fileColumn, nFileColumns);
    n = strlen(text) + 1;
    text = realloc(text,n*sizeof(char));
    return text;
  }

  int hasRowIdlabel = 0;
  if ( fileColumn == nFileColumns )
    hasRowIdlabel = 1;

  nFileColumns = fileColumn;

  if (nFileColumns < 2)
  { const char text[] = "Error: less than two columns found in the file";
    const int m = strlen(text) + 1;
    char* error = malloc(m*sizeof(char));
    strcpy(error,text);
    return error;
  }

  /* Check if the other rows in the file have the same number of columns */
  fileRow = 2;
  while ((line = GetLine(file)))
  { if (line[0]=='\0') free(line); /* Ignore completely empty lines */
    else
      /* Parse the first column to find out what the rows contain */
    {
      fileColumn = 1; /* One more columns than tabs */
      for (s=line; (*s)!='\0'; s++) if(*s=='\t') fileColumn++;
      free(line);
      fileRow++;

      if (s==NULL)
      {
	int n = 1024;
        char* text = malloc(n*sizeof(char));
        sprintf (text, "Error reading line %d: Gene name is missing", fileRow);
	n = strlen(text) + 1;
        text = realloc(text,n*sizeof(char));
 	return text;
      }

      if (fileColumn < nFileColumns)
      {
	int n = 1024;
        char* text = malloc(n*sizeof(char));
        sprintf (text,
                 "Error reading line %d: only %d columns available (%d needed)",
                 fileRow, fileColumn, nFileColumns);
	n = strlen(text) + 1;
        text = realloc(text,n*sizeof(char));
	return text;
      }

      if (fileColumn > nFileColumns)
      {
	int n = 1024;
        char* text = malloc(n*sizeof(char));
        sprintf (text,
                 "Error reading line %d: %d columns given (%d needed)",
                 fileRow, fileColumn, nFileColumns);
	n = strlen(text) + 1;
        text = realloc(text,n*sizeof(char));
	return text;
      }
    }
  }

  int nNumRows = fileRow-1;
  *nRows = nNumRows;

  int nNumCols = nFileColumns-1;
  *nCols = nNumCols;


  /* Read the first line into a string */
  fseek (file, 0, SEEK_SET);
  line = GetLine(file);

  if(!line)
  {
    char text[] = "Error finding UniqID keyword";
    const int m = strlen(text) + 1;
    char* error = malloc(m*sizeof(char));
    strcpy(error,text);
    return error;
  }

  /* Allocate space for column labels and save them */

  char **colLabel = (char **)malloc(nNumCols*sizeof(char*)); // don't save rowIdlabel

  column = 0;
  s = line;

  if ( !hasRowIdlabel )
  {
    token = s;
    s = tokenize(s);
    n = strlen(token);
    colLabel[column] = malloc((n+1)*sizeof(char));
    strcpy(colLabel[column],token);
    column++;
  }

  while (column < nNumCols)
  {
    token = s;
    s = tokenize(s);
    n = strlen(token);
    colLabel[column] = malloc((n+1)*sizeof(char));
    strcpy(colLabel[column],token);
    column++;
  }
  free(line);

  /* Allocate space for data */
  double **data = (double **)malloc(nNumRows*sizeof(double*));
  for (row = 0; row < nNumRows; row++)
    data[row] = malloc(nNumCols*sizeof(double));

  /* Allocate space for row labels */
  char **rowLabel = (char **)malloc(nNumRows*sizeof(char*));

  // read numerical data and row labels
  row = 0;
  while ((line=GetLine(file)))
  {
    if (strlen(line) > 1) /* Ignore completely empty lines */
    {
      column = 0;
      fileColumn = 0;
      s = line;

      while (s)
      {
	token = s;
	s = tokenize(s);

	if (fileColumn==0)
	{
	  const int n = strlen(token) + 1;
	  rowLabel[row] = malloc(n*sizeof(char));
	  strcpy (rowLabel[row],token);
	}
	else
	{
	  char* error = NULL;
	  data[row][column] = 0;

	  if (token[0]!='\0') /* Otherwise it is a missing value */
	  {
	    double number = strtod(token, &error);

	    if (!(*error))
	      data[row][column] = number;
	  }
	  column++;
	}

	fileColumn++;
      }
      row++;
      free(line);
    }
  }

  *rowNames = rowLabel;
  *colNames = colLabel;
  *matrix   = data;
  fclose(file);

  return 0;
}


// linked list of table rows
typedef struct tblRow
{
  //int nCols; the user needs to know the number of elements in row array
  char **row;
  struct tblRow *next;
} tblRow_t;

//--------------------------------- readCharTbl ---------------------------
int readCharTbl( const char *inFile, char ****tbl, int *nRows, int *nCols)
/*
  Reads a tab delimited character table with no header.

  (*tbl)[i] is the i-th row of a table (hence char ** array)
  nRows and nCols is the number of rows, columns, respectively.

  Return value is 0 on success and 1 on error.
*/
{
  char* line;
  char* s;
  char* token;

  FILE *file = fOpen(inFile,"r");

  // Parse the first line to figure out the number of columns of the table
  if( !( line=GetLine(file) ) )
  {
    fprintf(stderr,"Error in %s at line %d: Attempt to read empty file %s\n",__FILE__,__LINE__, inFile);
    return 1;
  }

  while(line[0]=='\0') // Ignore completely empty lines
  {
    free(line);
    line = GetLine(file);
    if(!line)
    {
      fprintf(stderr,"Error in %s at line %d: Failed to find the first line in %s\n",__FILE__,__LINE__, inFile);
      return 1;
    }
  }

  s = tokenize(line);
  *nCols = 1;

  while (s)
  {
    s = tokenize(s);
    (*nCols)++;
  }
  free(line);

  //fprintf(stderr,"in readCharTbl() nCols=%d\n",*nCols);

  // rewind
  fseek (file, 0, SEEK_SET);

  tblRow_t *rowRootNode;
  NEW(rowRootNode, tblRow_t);

  tblRow_t *rowNode  = 0;
  tblRow_t *oldRowNode = rowRootNode;
  char *str;

  int rowIdx = 0;
  int colIdx = 0;

  while ( ( line=GetLine(file) ) )
  {
    if (strlen(line) > 1) // Ignore completely empty lines
    {

      NEW(rowNode, tblRow_t);
      oldRowNode->next = rowNode;
      MALLOC(rowNode->row, char**, (*nCols)*sizeof(char*));

      s = line;
      colIdx = 0;
      while (s)
      {
	token = s;
	s = tokenize(s);

	if (token[0]!='\0') /* Otherwise it is a missing value */
	{
	  STRDUP(str, token);
	  rowNode->row[colIdx++] = str;
	}

	if ( colIdx > *nCols )
	{
	  fprintf(stderr,"ERROR in %s at line %d: Number of columns in row %d exceeds nCols=%d\n",__FILE__,__LINE__,rowIdx,*nCols);
	  return 1;
	}
      }
    }
    rowIdx++;
    oldRowNode = rowNode;
    free(line);
  }

  rowNode = rowRootNode->next;

  //fprintf(stderr,"in readCharTbl() rowIdx=%d\n", rowIdx);

  MALLOC((*tbl), char***, rowIdx*sizeof(char**));
  int i;
  for ( i = 0; i < rowIdx; i++, rowNode=rowNode->next )
  {
    (*tbl)[i] = rowNode->row;
  }

  *nRows = rowIdx;
  fclose(file);

  return 0;
}

//---------------------------------------------------------- printCharTbl ----
void printCharTbl(char ***tbl, int nRows, int nCols)
// prints to stdout character table
{
  int i, j;
  for ( i = 0; i < nRows; ++i )
  {
    for ( j = 0; j < nCols; ++j )
      printf("%s ", tbl[i][j]);
    printf("\n");
  }
  printf("\n");
}


//---------------------------------------------------------- writeCharTbl ----
void writeCharTbl(char *inFile, const char ***tbl, int nRows, int nCols )
// writes character table to a file
{
  FILE *out = fOpen(inFile, "w");
  int i, j;
  for ( i = 0; i < nRows; ++i )
  {
    for ( j = 0; j < nCols; ++j )
      fprintf(out, "%s ", tbl[i][j]);
    fprintf(out, "\n");
  }
  fprintf(out, "\n");

  fclose(out);
}
