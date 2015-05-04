#include <Rcpp.h>
#include <zlib.h>
#include "vcfRCommon.h"

using namespace Rcpp;

// Number of records to report progress at.
const int nreport = 1000;

// Size of the block of memory to use for reading.
#define LENGTH 0x1000


// [[Rcpp::export]]
Rcpp::NumericVector file_stats( std::string filename,
                                char sep = '\t',
                                int nrows = -1,
                                int skip = 0,
                                int verbose = 1) {

  // Initialize return datastructure.
  Rcpp::NumericVector stats(3);
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

    std::string mystring(reinterpret_cast<char*>(buffer));  // Recast buffer from char to string.
    mystring = lastline + mystring;

    std::vector < std::string > svec;  // Initialize vector of strings for parsed buffer.
    char line_split = '\n'; // Must be single quotes!
    vcfRCommon::strsplit(mystring, svec, line_split);
  
    // Scroll through lines derived from the buffer.
      for(int i=0; i < svec.size() - 1; i++){
      stats[0]++;
      
      int nrec = stats[0];
      if( nrec % nreport == 0 & verbose == 1){
        Rcpp::Rcout << "\rProcessed line: " << stats[0];
      }
      
      if( stats[0] > skip){
        if(nrows >= 0 & stats[1] < nrows){
          stats[1]++;
        } else if (nrows < 0){
          stats[1]++;
        }
      }
    }
    // Manage the last line.
    lastline = svec[svec.size() - 1];

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

  
  // Count columns in last line.
  std::vector < std::string > column_vec;  // Initialize vector of strings for parsed buffer.
  char col_split = sep; // Must be single quotes!
  vcfRCommon::strsplit(svec[svec.size() - 1], column_vec, col_split);
  stats[2] = column_vec.size();
  }

  gzclose (file);
  
  if( verbose == 1){
    Rcpp::Rcout << "\nCompleted: " << stats[0] << " lines.\n";
  }

  return stats;
}



// [[Rcpp::export]]
Rcpp::StringMatrix read_matrix( std::string filename,
                                char sep = '\t',
                                int nrows = 1,
                                int ncols = 1,
                                int skip = 0,
                                int verbose = 1) {

  Rcpp::StringMatrix matrix(nrows, ncols);
  
  
  return matrix;
}