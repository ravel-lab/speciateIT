//
// Given a CSV file of a numeric table, tbl, (for example a rlpp table) and the
// number of bins, EMD's between the first row, x, of tbl and the i-th (i=1, ...
// , nr) row y are computed, where nr is the number of rows of tbl.
//
// Pawel Gajer
// Tuesday, May 31 2022
//
// Compile with
// make -f Makefile_emd_test

#include <iostream>
#include <getopt.h>
#include "csv.hpp"
#include "emd.h"
#include "CppStatUtilities.hh"
#include "CStatUtilities.h"
//#include "IOUtilities.hh"
#include "CUtilities.h"

using namespace csv;
using namespace std;


//----------------------------------------------------------- printUsage ----
void printUsage( const char *s )
{
    cout << endl

         << "USAGE " << endl << endl
         << s << "-f <file > [Options]" << endl << endl

         << "\twhere"
         << endl
         << "\tOptions:\n"
         << "\t-f|-in-file <iFile>   - A CSV file holding a numeric table 'tbl'.\n"
         << "\t-o|--out-file <oFile> - The name of an output file. If not supplied, the EMD's will be sent to stdout.\n"
         << "\t-b|--num-bins <nbins> - The number of bins for histograms of the rows of tbl. Default value: 100.\n"
         << "\t-v                    - A verbose mode.\n\n"
         << "\t-h|--help             - This message\n\n"

         << "\tThe format of the outout data is a one column table with EMD's of the first row of 'tbl' and all other rows.\n"

         << "\n\tExample: \n"
         << s << " -f sIT_models " << endl << endl;
}


//----------------------------------------------- printHelp ----
void printHelp( const char *s )
{
    cout << endl
         << "emd_test estimates the EMD's between the first row of 'tbl' and all other rows." << endl;

    printUsage(s);
}


//================================================= inPar2_t ====
//! holds input parameters
class inPar_t
{
public:
  inPar_t();
  ~inPar_t();

  char *inFile;      /// The name of the input CSV file holding a numerical table.
  char *outFile;     /// The name of the output file.
  int nbins;         /// The number of bins in the histograms of each pairs of rows of the input table.
  bool verbose;      /// Set to true for verbose mode.

  void print();
};


//============================== local sub-routines =========================

void parseArgs( int argc, char ** argv, inPar_t *p );

double* row_to_darray(int i, vector<CSVRow> &rows, int ncols);
double dist(feature_t *F1, feature_t *F2);
double hist_emd(double *x,
                int     nx,
                double *y,
                int     ny,
                int     nbins);

//============================== main ======================================
int main(int argc, char **argv)
{
    //-- setting up init parameters
    inPar_t *inPar = new inPar_t();

    //-- parsing input parameters
    parseArgs(argc, argv, inPar);

    if ( inPar->verbose )
      inPar->print();

    #if 0
    if ( !inPar->inFile )
    {
      fprintf(stderr, "ERROR: Missing input file. Please specify it using the -f flag.");
      exit(EXIT_FAILURE);
    }
    #endif

    //
    // Read input table
    //
    inPar->inFile = strdup("/Users/pgajer/devel/speciateIT/data/vag_sIT_models_V3V4/rLpps/Abiotrophia_defectiva.csv");

    CSVFormat format;
    format.trim({ ' ' })
      .variable_columns(VariableColumnPolicy::THROW) // .variable_columns(false)
      .no_header()
      .delimiter(',');

    CSVReader reader(inPar->inFile, format);
    vector<CSVRow> rows(reader.begin(), reader.end());
    int nrows = rows.size();
    auto info = get_file_info(inPar->inFile);
    int ncols = info.n_cols;
    cout << "Found " << nrows << " rows and " << ncols << " columns" << endl << endl;

    // Extracting the first row (only numeric values) of the table
    double *x = row_to_darray(0, rows, ncols);
    #if 0
    for ( int i = 0; i < ncols-1; i++ )
    {
      if ( i < ncols-2 )
        cout << x[i] << ",";
      else if ( i == ncols-2 )
        cout << x[i] << endl;
    }
    #endif


    // Remaining rows
    int nnum = ncols - 1;
    for( int j = 1; j < nrows; j++ )
    {
      double *y = row_to_darray(j, rows, ncols);
      double d =  hist_emd(x, nnum, y, nnum, inPar->nbins);
      cout << "emd(x,y[" << j << "])=" << d << endl;

      #if 0
      for ( int i = 0; i < ncols-1; i++ )
      {
        if ( i < ncols-2 )
          cout << x[i] << ",";
        else if ( i == ncols-2 )
          cout << x[i] << endl;
      }
      #endif
      free(y);
    }

    free(x);

    return EXIT_SUCCESS;
}




//----------------------------------------------------------- parseArgs ----
//! parse command line arguments
void parseArgs( int argc, char ** argv, inPar_t *p )
{
    int c, errflg = 0;
    optarg = NULL;

    static struct option longOptions[] = {
      {"in-file"   ,required_argument, 0, 'f'},
      {"out-file"  ,required_argument, 0, 'o'},
      {"num-bins"  ,required_argument, 0, 'b'},
      {"help"            ,no_argument, 0,   0},
      {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv,"b:f:o:vh",longOptions, NULL)) != -1)
      switch (c)
      {
        case 'b':
          p->nbins = atoi(optarg);
          break;

        case 'f':
          p->inFile = strdup(optarg);
          break;

        case 'o':
          p->outFile = strdup(optarg);
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
}

//------------------------------------------------- constructor ----
inPar_t::inPar_t()
{
    inFile   = NULL;
    outFile  = NULL;
    nbins    = 100;
}

//------------------------------------------------- constructor ----
inPar_t::~inPar_t()
{
  if ( inFile )
    free(inFile);

  if ( outFile )
    free(outFile);
}

//------------------------------------------------------- print ----
void inPar_t::print()
{
    if ( inFile )
      cerr << inFile << endl;
    else
      cerr << "MISSING" << endl;

    if ( outFile )
      cerr << outFile << endl;
    else
      cerr << "MISSING" << endl;

    cerr << "nbins: " << nbins << endl;
}


//----------------------------------------------------------- dist ----
double dist(feature_t *F1, feature_t *F2)
{
    double dX = *F1 - *F2;
    return sqrt( dX * dX );
}

// ------------------------------------------------------- hist_emd ----
/*!
  Computes the Earth Movers Distance (EMD) between histograms of two arrays

  \param x   - A double array of size xSize
  \param nx  - The number of elements in x.
  \param y   - A double array of size ySize
  \param ny  - The number of elements in y.
  \param nbins - The number of bins for the histograms of x and y.

  Returns EMD between the histograms of x and y.
*/
double hist_emd(double *x,
                int     nx,
                double *y,
                int     ny,
                int     nbins)
{

    double *hmid;
    MALLOC(hmid, double*, nbins * sizeof(double));

    double *hx;
    MALLOC(hx, double*, nbins * sizeof(double));

    double *hy;
    MALLOC(hy, double*, nbins * sizeof(double));

    hist(x, nx, y, ny, nbins, hmid, hx, hy);

    signature_t s1 = { nbins, (feature_t *)hmid, hx };
    signature_t s2 = { nbins, (feature_t *)hmid, hy };

    double e = emd(&s1, &s2, dist, 0, 0);

    free(hmid);
    free(hx);
    free(hy);

    return e;
}

//----------------------------------------------------------- row_to_darray ----
/*!
  Creates a double array from the i-th element of vector<CSVRow>.

  The first entry of the vector is the row name, which is skipped and only numeric
  values are returned.

  \param i     - The i-th element of 'rows' array.
  \param rows  - The rows of the input table.
  \param ncols - The number of columns in the table.

  \return An array of numeric values of the corresponding row.
*/
double* row_to_darray(int i, vector<CSVRow> &rows, int ncols)
{
    int nnum = ncols - 1; // The row constains a row ID and then numeric values; This is the number of numeric values in the row.

    double *darray;
    MALLOC(darray, double*, nnum * sizeof(double));

    CSVRow row = rows[i];
    i = 0; // field counter
    for ( CSVField& field: row )
    {
      if ( i > 0 )
        darray[i-1] = field.get<double>();
      i++;
    }

    return darray;
}
