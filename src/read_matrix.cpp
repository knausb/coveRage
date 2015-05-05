#include <Rcpp.h>
#include <zlib.h>
#include "vcfRCommon.h"

using namespace Rcpp;

// Number of records to report progress at.
const int nreport = 1000;

// Size of the block of memory to use for reading.
#define LENGTH 0x1000


// Modified from:
// http://stackoverflow.com/a/236803
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
//    return elems;
}



//' @title File input/output
//' @name File input/output
//' 
//' @description Fast but featureless input of tabular data in either *.txt or *.gz format.
//' 
//' 
//' @rdname read_matrix
//' @aliases file_stats
//' 
//' @param filename name of a file
//' @param sep character which delimits columns
//' @param nrows number of rows to read
//' @param skip number of rows to skip
//' @param verbose should verbose output be generated
//' 
//' @return \strong{file_stats} returns a three element vector.
//' 'Total_rows' reports the total number of rows read.
//'  This is either the number of rows in the file or the number of skipped rows and the number of rows read in.
//'  'Rows' is the number of rows read in.
//'  This is either the same as nrows or however many rows were read in after skip and before the end of the file (when less than nrows).
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector file_stats( std::string filename,
                                char sep = '\t',
                                int nrows = -1,
                                int skip = 0,
                                int verbose = 1) {

  // Initialize return datastructure.
//  Rcpp::NumericVector stats(3);
  Rcpp::IntegerVector stats(3);
  stats.names() = Rcpp::StringVector::create("Total_rows", "Rows", "Columns");
  
  // Open the file stream.
  gzFile file;
  file = gzopen (filename.c_str(), "r");
  if (! file) {
    Rcpp::Rcerr << "gzopen of " << filename << " failed: " << strerror (errno) << ".\n";
    return stats;
  }

  // Scroll through buffer.
  std::string lastline = "";
  while (1) {
    Rcpp::checkUserInterrupt();
    int err;
    int bytes_read;
    char buffer[LENGTH];

    bytes_read = gzread (file, buffer, LENGTH - 1);
    buffer[bytes_read] = '\0';
    
    // Recast buffer from char to string.
    std::string mystring(reinterpret_cast<char*>(buffer));
    
    // Add last line to begining of new input.
    mystring = lastline + mystring;
    
    // Initialize vector of strings for parsed buffer.
    // Delimit buffer with newline characters.
    std::vector < std::string > line_vec;
    char line_split = '\n'; // Must be single quotes!
    split(mystring, line_split, line_vec);

    // Scroll through lines derived from the buffer.
    for(int i=0; i < line_vec.size() - 1; i++){
      // Increment line counter
      stats[0]++;

      // Count columns in first row past skipped rows.
      if( stats[0] == skip + 1 ){
        std::vector < std::string > column_vec;  // Initialize vector of strings for parsed buffer.
        split(line_vec[i], sep, column_vec);
        stats[2] = column_vec.size();
      }
      
      if( stats[0] > skip ){
        stats[1]++;
      }

      if( nrows > 0 & stats[1] > nrows - 1 ){
        gzclose (file);
        if( verbose == 1){
          Rcpp::Rcout << "\nCompleted: " << stats[0] << " lines.\n";
        }
        return stats;
      }

      if( stats[0] % nreport == 0 & verbose == 1){
          Rcpp::Rcout << "\rProcessed line: " << stats[0];
      }
    }
    // Manage the last line.
    lastline = line_vec[line_vec.size() - 1];

    // Check for EOF or errors.
    if (bytes_read < LENGTH - 1) {
      if (gzeof (file)) {
        break;
      }
      else {
        const char * error_string;
        error_string = gzerror (file, & err);
        if (err) {
          Rcpp::Rcerr << "Error: " << error_string << ".\n";
          return stats;
        }
      }
    } 
  }

  gzclose (file);
  
  if( verbose == 1){
    Rcpp::Rcout << "\nCompleted: " << stats[0] << " lines.\n";
  }

  return stats;
}



//' @rdname read_matrix
//' @aliases read_matrix
//' 
//' @param ncols number of columns for the matrix
//' 
//' @return \strong{read_matrix} returns a matrix of strings of dimension specified by nrows and ncols.
//' 
//' @seealso
//' \href{http://cran.r-project.org/web/packages/readr/index.html}{readr}
//' \href{http://cran.r-project.org/web/packages/data.table/index.html}{data.table::fread}
//'
//' @export
// [[Rcpp::export]]
Rcpp::StringMatrix read_matrix( std::string filename,
                                char sep = '\t',
                                int nrows = 1,
                                int ncols = 1,
                                int skip = 0,
                                int verbose = 1) {

  // Initialize return data structure.
  Rcpp::StringMatrix mymatrix(nrows, ncols);
  Rcpp::IntegerVector stats(3); // Helper.
  
  // Open the file stream.
  gzFile file;
  file = gzopen (filename.c_str(), "r");
  if (! file) {
    Rcpp::Rcerr << "gzopen of " << filename << " failed: " << strerror (errno) << ".\n";
    return mymatrix;
  }
 
 
   // Scroll through buffer.
  std::string lastline = "";
//  int nline = 0;
//  int rownum = 0;
  while (1) {
    Rcpp::checkUserInterrupt();
    int err;
    int bytes_read;
    char buffer[LENGTH];
    
    bytes_read = gzread (file, buffer, LENGTH - 1);
    buffer[bytes_read] = '\0';

    std::string mystring(reinterpret_cast<char*>(buffer));  // Recast buffer from char to string.
    mystring = lastline + mystring;

    std::vector < std::string > line_vec;  // Initialize vector of strings for parsed buffer.
    char line_split = '\n'; // Must be single quotes!
//    vcfRCommon::strsplit(mystring, svec, line_split);
    split(mystring, line_split, line_vec);
  
    // Scroll through lines derived from the buffer.
      for(int i=0; i < line_vec.size() - 1; i++){
        // Increment line counter
        stats[0]++;
        
        if( stats[0] % nreport == 0 & verbose == 1){
          Rcpp::Rcout << "\rProcessed line: " << stats[0];
        }

        if( stats[0] >= skip + 1 ){
          // Load line into matrix.
          std::vector < std::string > column_vec;  // Initialize vector of strings for parsed buffer.
          char col_split = sep; // Must be single quotes!
          vcfRCommon::strsplit(line_vec[i], column_vec, col_split);
          split(line_vec[i], sep, column_vec);
          
          if(mymatrix.ncol() > column_vec.size()){
            Rcerr << "Warning: more matrix rows than input elements on line: " << stats[0] << "\n";
            Rcerr << "Using as many input elements that fit.\n";
            for(int j = 0; j < column_vec.size(); j++){
              mymatrix(stats[1], j) = column_vec[j];
            }
          } else if (mymatrix.ncol() < column_vec.size()){
            Rcerr << "Warning: more input elements than matrix rows on line: " << stats[0] << "\n";
            Rcerr << "Using as many input elements that fit.\n";
            for(int j = 0; j < mymatrix.ncol(); j++){
              mymatrix(stats[1], j) = column_vec[j];
            }
          } else if (mymatrix.ncol() == column_vec.size()){
            for(int j = 0; j < mymatrix.ncol(); j++){
              mymatrix(stats[1], j) = column_vec[j];
            }          
          }
          stats[1]++;
          if( nrows > 0 & stats[1] > nrows - 1 ){
            gzclose (file);
            if( verbose == 1){
              Rcpp::Rcout << "\nCompleted: " << stats[0] << " lines.\n";
            }
            return mymatrix;
          }
        }
      }
    // Manage the last line.
    lastline = line_vec[line_vec.size() - 1];

    // Check for EOF or errors.
    if (bytes_read < LENGTH - 1) {
      if (gzeof (file)) {
        break;
      }
      else {
        const char * error_string;
        error_string = gzerror (file, & err);
        if (err) {
          Rcpp::Rcerr << "Error: " << error_string << ".\n";
          return mymatrix;
        }
      }
    }
  }

  gzclose (file);
  
  if( verbose == 1){
    Rcpp::Rcout << "\nCompleted: " << stats[0] << " lines.\n";
  }

  return mymatrix;
}


