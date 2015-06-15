#include <Rcpp.h>
#include <zlib.h>
// #include "vcfRCommon.h"

using namespace Rcpp;

// Number of records to report progress at.
const int nreport = 1000;

// Size of the block of memory to use for reading.
#define LENGTH 0x1000


std::vector< std::string > strsplit(std::string mystring, std::string delimiter){
    // Modified from:
    // http://stackoverflow.com/a/14266139
    size_t pos = 0;
    std::string token;
    std::vector< std::string > newvector;
    
    while ((pos = mystring.find(delimiter)) != std::string::npos) {
      token = mystring.substr(0, pos);
      newvector.push_back(token);
      mystring.erase(0, pos + delimiter.length());
    }
    newvector.push_back(mystring);
    return newvector;
}


void proc_line(Rcpp::StringVector mystring,
               std::string line,
               std::string sep,
               Rcpp::IntegerVector cols){

  std::vector < std::string > temp_vec = strsplit(line, sep);

  int i = 0;
  while( i < cols.size() ){
//    Rcpp::Rcout << "cols[ cols.size() - 1 ]: " << cols[ cols.size() - 1 ] << "\n";
//    Rcpp::Rcout << "temp_vec.size: " << temp_vec.size() << "\n";
    if( cols[ cols.size() - 1 ] > temp_vec.size() - 1 ){
//      Rcpp::Rcerr << "More column elements selected than exist!\n";
      mystring[0] = "ERROR";
      return;
    } else {    
      mystring[i] = temp_vec[ cols[i] ];
      i++;
    }
  }
  
}




//' @title File input
//' @name File input
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
//' @details \strong{file_stats} returns a three element vector containing a summary of a file's contents.
//' 'Total_rows' reports the total number of rows read.
//'  This is either the number of rows in the file or the number of skipped rows and the number of rows read in.
//'  'Rows' is the number of rows read in.
//'  This is either the same as nrows or however many rows were read in after skip and before the end of the file (when less than nrows).
//'  'Columns' is the number of columns resulting after delimiting with sep.
//'  This information is intended to be used with read_matrix.
//'  
//' @return \strong{file_stats} returns a three element vector.
//'  
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector file_stats( std::string filename,
                                std::string sep = "\t",
                                int nrows = -1,
                                int skip = 0,
                                int verbose = 1) {

  // Initialize return datastructure.
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
    
    // Delimit buffer with newline characters.
    std::vector < std::string > line_vec = strsplit(mystring, "\n");

    // Scroll through lines derived from the buffer.
    for(int i=0; i < line_vec.size() - 1; i++){
      Rcpp::checkUserInterrupt();
      // Increment line counter
      stats[0]++;

      // Count columns in first row past skipped rows.
      if( stats[0] == skip + 1 ){
        std::vector < std::string > column_vec = strsplit(line_vec[i], sep);
        stats[2] = column_vec.size();
//        Rcpp::Rcout << "Row for column counting:\n";
//        Rcpp::Rcout << line_vec[i];
//        Rcpp::Rcout << "\n";
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
//' @param cols vector of column numbers to include in the matrix
//' 
//' @details \strong{read_matrix} returns a matrix of strings of dimension specified by nrows and cols.
//' The cols parameter is a vector of integers specifying which columns to read in.
//'
//' @return \strong{read_matrix} returns a matrix of strings
//'
//' 
//' @seealso
//' \href{http://cran.r-project.org/web/packages/readr/index.html}{readr}
//' \href{http://cran.r-project.org/web/packages/data.table/index.html}{data.table::fread}
//'
//' @export
// [[Rcpp::export]]
Rcpp::StringMatrix read_matrix( std::string filename,
                                std::string sep = "\t",
                                int nrows = 1,
                                Rcpp::IntegerVector cols = 0,
                                int skip = 0,
                                int verbose = 1) {

  cols.sort();
  if(cols[0] == 0){
    Rcerr << "User must specify which (positive integer) columns to extract from the file.\n";
    Rcpp::StringMatrix mymatrix(1,1);
    mymatrix(0,0) = NA_STRING;
    return mymatrix;
  }
  cols = cols - 1; // R is 1-based
  
  // Initialize return data structure.
//  Rcpp::StringMatrix mymatrix(nrows, ncols);
  Rcpp::StringMatrix mymatrix(nrows, cols.size());
  
  Rcpp::IntegerVector stats(3); 
  // stats is a helper where:
  // the first element is line number,
  // the second element counts rows included
  // and the third element isn't actually used here.
  
  // Open the file stream.
  gzFile file;
  file = gzopen (filename.c_str(), "r");
  if (! file) {
    Rcpp::Rcerr << "gzopen of " << filename << " failed: " << strerror (errno) << ".\n";
    return mymatrix;
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

    // Parse buffer on newline.
    std::vector < std::string > line_vec = strsplit(mystring, "\n");

    // Scroll through lines derived from the buffer.
    for(int i=0; i < line_vec.size() - 1; i++){
      Rcpp::checkUserInterrupt();
      // Increment line counter
      stats[0]++;
      
      if( stats[0] % nreport == 0 & verbose == 1){
        Rcpp::Rcout << "\rProcessed line: " << stats[0];
      }

      if( stats[0] >= skip + 1 ){
        // Parse line base on sep.
        Rcpp::StringVector xx_row = mymatrix(stats[1], Rcpp::_);
        proc_line(xx_row, line_vec[i], sep, cols);
        if(xx_row[0] == "ERROR"){
          Rcpp::Rcerr << "More columns selected than exist!\n";
          Rcpp::Rcerr << "File row: " << stats[0] << "\n";
          Rcpp::Rcerr << "Returning matrix up to bad row.\n";
          return mymatrix;
        }        
        mymatrix(stats[1], Rcpp::_) = xx_row;

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


        /*
        std::vector < std::string > column_vec = strsplit(line_vec[i], sep);
          
        if(mymatrix.ncol() > column_vec.size()){
          Rcerr << "Warning: more matrix columns than input elements on line: " << stats[0] << "\n";
          Rcerr << "Using as many input elements that fit.\n";
          for(int j = 0; j < column_vec.size(); j++){
            mymatrix(stats[1], j) = column_vec[j];
          }
        } else if (mymatrix.ncol() < column_vec.size()){
          Rcerr << "Warning: more input elements than matrix columns on line: " << stats[0] << "\n";
          Rcerr << "Using as many input elements that fit.\n";
          for(int j = 0; j < mymatrix.ncol(); j++){
            mymatrix(stats[1], j) = column_vec[j];
          }
        } else if (mymatrix.ncol() == column_vec.size()){
          for(int j = 0; j < mymatrix.ncol(); j++){
            mymatrix(stats[1], j) = column_vec[j];
          }          
        }
        */
        
        
        