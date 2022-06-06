//
// Reads a numeric (tab or comma delimited numerical table with header and row names) and then prints it to stdout
//
// Usage: readTable_test -f <file>
//
// clang -std=c++11 -lc++abi -lstdc++ -I/Users/pgajer/devel/speciateIT/src readTable_test.cc -o readTable_test
//

// Example
//
// readTable_test ~/devel/speciateIT/data/test_data/sample_3x4_num_tbl_with_rownames_and_no_header.csv


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "IOUtilities.h"
#include "CUtilities.h"

using namespace std;

//----------------------------------------------------------------- comma_tokenize ----
/*!
  Returns the pointer to a comma separated token
*/
static char* comma_tokenize(char* s)
{
  char* p = s;
  while (1)
  {
    if (*p=='\0') return NULL;
    if (*p==',')
    {
      *p = '\0';
      return p+1;
    }
    p++;
  }
  /* Never get here */
  return NULL;
}

/*!
  Finds tokens separated by the commans. The first token is assumed to be a string
  and the remaining tokens are double's.
*/
void comma_tokenize_line(char *line)
{
    if (strlen(line) == 0 )
      DIEM("The line is empty!");

    char* s = line;

    // extracting the string
    char* token = s;
    s = comma_tokenize(s);
    string rowID = string(token);
    fprintf(stderr, "rowID=%s", rowID.c_str());

    char* error = NULL;
    while (s)
    {
      token = s;
      s = comma_tokenize(s);

      if ( token[0] != '\0') /* Otherwise it is a missing value */
      {
        double f = strtod(token, &error);
        if ( !(*error) )
          fprintf(stderr, ",%lf", f);
        else
          perror("ERROR");
      }
    }
    fprintf(stderr, "\n");
}

/*!
  Counts the number of commas in a string
*/
int count_commas(char* s)
{
  char* p = s;
  if ( *p=='\0' ) return 0;

  int counter = 0;
  while (1)
  {
    if ( *p==',' )
      counter++;
    else if ( *p=='\0' )
      break;
    p++;
  }

  return counter;
}


/*!
  Reads a numeric CSV file with row IDs and with or without header.

  \param inFile - an input file name.
  \param tbl    - output table
  \param _nRows - a reference to the number of rows of tbl.
  \param _nCols - a reference to the number of columns of tbl.
  \param rowNames - row names.
  \param colNames - column names.

*/
void read_csv(const char *inFile,
              double ***tbl,
              int *_nRows,
              int *_nCols,
              vector<string> &rowNames,
              vector<string> &colNames )
{
    FILE *file = fOpen(inFile,"r");


    // Check the number of columns of the first line: ncol1
    // Counting commans of the first line

    char *line = GetLine(file);
    int ncol1 = count_commas(line);
    free(line);

    // Check the number of columns of the second line: ncol2
    line = GetLine(file);
    int ncol2 = count_commas(line);
    free(line);

    // Either ncol1=ncol2 or ncol1=ncol2-1
    // If this is not the case, throw an error
    if ( !(ncol1==ncol2 || ncol1==ncol2-1) )
    {
      fprintf(stderr, "ncol1=%d  ncol2=%d\n", ncol1, ncol2);
      DIEM("ncol1!=ncol2 AND ncol1!=ncol2-1");
    }

    // If ncol1=ncol2-1, we assume that the first row is the header
    // Otherwise, we assume that the file does not have the header.

    // Check that the number of columns in all other rows is ncol2
    int ncol;
    int rowCounter = 3;
    while ( (line = GetLine(file)) )
    {
      ncol = count_commas(line);
      free(line);
      if ( ncol != ncol2 )
      {
        fprintf(stderr, "ncol=%d  ncol2=%d\n", ncol, ncol2);
        DIEM("ncol != ncol2");
      }
      rowCounter++;
    }

    fprintf(stderr, "ncol1=%d  ncol2=%d\n", ncol1, ncol2);
    fprintf(stderr, "All rows starting at row 2 have the same number of comma separated %d fields\n", ncol);
    fprintf(stderr, "Found %d rows\n",rowCounter);

    // Go back to the first line
    fseek(file, 0, SEEK_SET);

    // if the file has the header, read it
    // For now assuming that the file does not have the header

    // Read all lines of the file storing the first token into rowNames
    // and all remaining ones in the table tbl

    // Parsing the first line
    line = GetLine(file);

    comma_tokenize_line(line);


    fclose(file);
}

int main(int argc, char *argv[])
{
    if (argc != 2) {
      fprintf(stderr, "Usage: %s <file>, where 'file' is the name of a tab or comma delimited file holding a numeric table with a header and row names.\n", argv[0]);
      exit(EXIT_FAILURE);
    }

    char *inFile = argv[1];

    vector<string> rowNames;
    vector<string> colNames;
    int nRows;
    int nCols;
    double **tbl;
    read_csv(inFile,
             &tbl,
             &nRows,
             &nCols,
             rowNames,
             colNames);

    exit(EXIT_SUCCESS);
}
