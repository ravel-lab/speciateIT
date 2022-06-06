//
// First example from
// https://github.com/vincentlaucsb/csv-parser#c-version
//

// compile with
// make -f Makefile_CSVReader_test



// Examples:

// CSVReader_test ../data/test_data/sample_3x4_num_tbl_with_rownames_and_no_header.csv

// [~/devel/speciateIT/src_tests]% time CSVReader_test ../data/vag_sIT_models_V3V4/MC6.log10cProb
// Found 3042 rows and 16385 columns
// CSVReader_test ../data/vag_sIT_models_V3V4/MC6.log10cProb  2.28s user 0.69s system 59% cpu 5.018 total

// [~/devel/speciateIT/src_tests]% time CSVReader_test ../data/vag_sIT_models_V3V4/MC7.log10cProb
// Found 3042 rows and 65537 columns
// CSVReader_test ../data/vag_sIT_models_V3V4/MC7.log10cProb  8.96s user 3.28s system 54% cpu 22.316 total

#include <iostream>
#include "csv.hpp"
#include "CUtilities.h"

using namespace csv;
using namespace std;


double* row_to_darray(int i, vector<CSVRow> &rows, int ncols);

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " [file]" << std::endl;
        exit(1);
    }

    std::string file = argv[1];

    CSVFormat format;
    format.trim({ ' ' })
      .variable_columns(VariableColumnPolicy::THROW) // .variable_columns(false)
      .no_header()
      .delimiter(',');

    //CSVReader reader("/Users/pgajer/devel/speciateIT/data/vag_sIT_models_V3V4/rLpps/Abiotrophia_defectiva.csv", format);
    //std::string file = argv[1];
    //string file = "/Users/pgajer/devel/speciateIT/data/test_data/sample_3x4_num_tbl_with_rownames_and_no_header.csv";

    #if 0
    CSVStat stats(file);
    auto col_names = stats.get_col_names();
    auto min = stats.get_mins(), max = stats.get_maxes(),
      means = stats.get_mean(), vars = stats.get_variance();

    for (size_t i = 0; i < col_names.size(); i++) {
      cout << col_names[i] << endl
                << "Min: " << min[i] << endl
                << "Max: " << max[i] << endl
                << "Mean: " << means[i] << endl
                << "Var: " << vars[i] << endl;
    }

    auto info = get_file_info(file);
    cout << file << endl
              << "Columns: " << internals::format_row(info.col_names, ", ")
              << "Dimensions: " << info.n_rows << " rows x " << info.n_cols << " columns" << endl
              << "Delimiter: " << info.delim << endl;
    #endif

    CSVReader reader(file, format);

    //int nrows = reader.n_rows();
    vector<CSVRow> rows(reader.begin(), reader.end());
    int nrows = rows.size();
    auto info = get_file_info(file);
    int ncols = info.n_cols;
    cout << "Found " << nrows << " rows and " << ncols << " columns" << endl << endl;


    #if 1
    //for (CSVRow& row: reader) // Input iterator
    for (CSVRow& row: rows) // Input iterator
    {
      int i = 0; // field counter
      for (CSVField& field: row)
      {
        // By default, get<>() produces a std::string.
        // A more efficient get<string_view>() is also available, where the resulting
        // string_view is valid as long as the parent CSVRow is alive
        if ( i > 0 )
          cout << field.get<double>() << ",";
        i++;
      }
      cout << endl;
    }
    #endif


    #if 0
    cout << "\n\n========================================\n"
         << "The first row: ";
    CSVRow row = rows[0];
    int i = 0; // field counter
    for (CSVField& field: row)
    {
      if ( i > 0 && i < ncols-1 )
        cout << field.get<double>() << ",";
      else if ( i == ncols-1 )
        cout << field.get<double>();
      i++;
    }
    cout << endl;

    cout << "\n\n========================================\n"
         << "Remaining rows" << endl;
    for( int j = 1; j < nrows; j++ )
    {
      row = rows[j];
      int i = 0; // field counter
      for (CSVField& field: row)
      {
        if ( i > 0 && i < ncols-1 )
          cout << field.get<double>() << ",";
        else if ( i == ncols-1 )
          cout << field.get<double>();
        i++;
      }
      cout << endl;
    }
    #endif

    cout << "\n\n========================================\n"
         << "Utilizing row_to_darray() to extract numeric values from the rows of the input table" << endl
         << "The first row: ";
    double *x = row_to_darray(0, rows, ncols);
    for ( int i = 0; i < ncols-1; i++ )
    {
      if ( i < ncols-2 )
        cout << x[i] << ",";
      else if ( i == ncols-2 )
        cout << x[i] << endl;
    }
    free(x);

    cout << "\n\n========================================\n"
         << "Remaining rows" << endl;
    for( int j = 1; j < nrows; j++ )
    {
      x = row_to_darray(j, rows, ncols);
      for ( int i = 0; i < ncols-1; i++ )
      {
        if ( i < ncols-2 )
          cout << x[i] << ",";
        else if ( i == ncols-2 )
          cout << x[i] << endl;
      }
      free(x);
    }

    return 0;
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
