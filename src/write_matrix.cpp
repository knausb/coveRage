#include <Rcpp.h>
#include <fstream>
#include <zlib.h>

// using namespace Rcpp;

// Number of records to report progress at.
const int nreport = 1000;

//' @name write_matrix
//' @title write_matrix
//' @rdname write_matrix
//' 
//' @description Write matrix data to a gzipped file.
//' 
//' @param filename filename for output
//' @param mymatrix matrix to be written to file
//' @param sep delimiting character
//' @param verbose should verbose output be generated?
//'
//' @details
//' Writes matrix data to a gzipped file delimited by sep.
//' Data is appended to the file.
//' This is intended to allow header information to be included by other calls.
//' It should also allow files to grow with incremental processes.
//'
//' @export
// [[Rcpp::export]]
void write_matrix( std::string filename,
                   Rcpp::StringMatrix mymatrix,
                   std::string sep = "\t",
                   int verbose = 1) {

  gzFile fi = gzopen( filename.c_str(), "ab" );
//  gzFile *fi = (gzFile *)gzopen(filename.c_str(),"ab");

  int i = 0;
  for( i = 0; i<mymatrix.nrow(); i++){
    Rcpp::checkUserInterrupt();
    
    if( i % nreport == 0 & verbose == 1){
      Rcpp::Rcout << "\rProcessed line: " << i;
    }
    
    std::string myelement = Rcpp::as< std::string >(mymatrix(i, 0));
    gzwrite(fi, (char *)myelement.c_str(), myelement.size());

    for(int j=1; j<mymatrix.ncol(); j++){
      Rcpp::checkUserInterrupt();
      std::string myelement = Rcpp::as< std::string >(mymatrix(i, j));
      gzwrite(fi, (char *)sep.c_str(), sep.size());
      gzwrite(fi, (char *)myelement.c_str(), myelement.size());
    }
    gzwrite(fi, "\n", strlen("\n"));
  }

  gzclose(fi);

  if( verbose == 1 ){
    Rcpp::Rcout << "\nWriting complete.\n";
    Rcpp::Rcout << "Processed " << i << "lines.\n";
  }

  return;
}
