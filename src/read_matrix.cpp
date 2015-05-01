#include <Rcpp.h>
#include <fstream>
#include <zlib.h>
//#include <errno.h>
#include "strsplit.h"

//using namespace Rcpp;

// Number of records to report progress at.
const int nreport = 1000;

// Size of the block of memory to use for reading.
#define LENGTH 0x1000


// [[Rcpp::export]]
Rcpp::NumericVector vcf_stats(std::string x) {
  // Scroll through file and collect the number of
  // rows and columns for the returned matrix.

  Rcpp::NumericVector stats(2);
  stats.names() = Rcpp::StringVector::create("rows", "columns");

  gzFile file;
  file = gzopen (x.c_str(), "r");
  if (! file) {
    Rcpp::Rcerr << "gzopen of " << x << " failed: " << strerror (errno) << ".\n";
    return stats;
  }

  // Scroll through buffers.
  std::string lastline = "";
  while (1) {
    Rcpp::checkUserInterrupt();
    int err;
    int bytes_read;
    char buffer[LENGTH];
    bytes_read = gzread (file, buffer, LENGTH - 1);
    buffer[bytes_read] = '\0';

    std::string mystring(reinterpret_cast<char*>(buffer));  // Recast buffer as a string.
    mystring = lastline + mystring;
    std::vector < std::string > svec;  // Initialize vector of strings for parsed buffer.
    
    char line_split = '\n'; // Must be single quotes!
//    strsplit a1;
//    strsplit::strsplit(mystring, svec, line_split);

    // Scroll through lines derived from the buffer.
    for(int i=0; i < svec.size() - 1; i++){
//      stat_line(stats, svec[i]);
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
  }
  gzclose (file);

  return stats;
}



