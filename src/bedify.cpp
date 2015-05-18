#include <Rcpp.h>
//using namespace Rcpp;


//' @title Parse data by a bed file
//' @rdname bedify
//' 
//' @description seperate a data matrix by using bed format data.
//' 
//' @param myBed matrix of bed format data
//' @param myData StringMatrix to be sorted
//' 
//' @details
//' 
//' Bed format data contain at least four columns.
//' The first column indicates the chromosome (i.e., supercontig, scaffold, contig, etc.).
//' The second cotains the starting positions.
//' The third the ending positions.
//' The fourth are the names of the features.
//' All subsequent columns are ignored here.
//' In an attempt to optimize performance the data are expected to be formatted as a character matrix.
//' The starting and end positions are converted to numerics internally.
//' 
//' 
//' \href{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}{Bed format} at UCSC
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List bedify( Rcpp::StringMatrix myBed, Rcpp::IntegerMatrix myData ) {
  
  Rcpp::List myList(myBed.nrow());
  Rcpp::StringVector myNames(myBed.nrow());
  
  // Scroll over bed rows (features).
  for(int i=0; i<myBed.nrow(); i++){
    Rcpp::checkUserInterrupt();
//    Rcpp::Rcout << "Feature: " << i << "\n";
    
    // Count rows in feature
    std::string temp = Rcpp::as< std::string >(myBed(i,1));
    int start = stoi(temp);
    temp = Rcpp::as< std::string >(myBed(i,2));
    int end = stoi(temp);
    
    // Manage if on reverse strand
    if(end < start){
      int tmp = start;
      start = end;
      end = tmp;
    }
    
    int nrows = end  - start + 1;
//    Rcpp::Rcout << "nrows: " << nrows << "\n";
    
    // Initialize return matrix.
    Rcpp::IntegerMatrix myMatrix(nrows, myData.ncol());
    std::fill(myMatrix.begin(), myMatrix.end(), NA_INTEGER);
    for(int j=0; j<myMatrix.nrow(); j++){
      myMatrix(j,0) = start + j;
    }

    int j=0;  // matrix to bedify row counter
//    int k = 0; // Out matrix row counter.
    
    // Increment to feature beginning.
    while(myData(j,0) < start){
      Rcpp::checkUserInterrupt();
      j++;
    }
    
    // Process to the end.
    while( myData(j,0) <= end & j < myData.nrow() ){
      Rcpp::checkUserInterrupt();
//      Rcpp::Rcout << "myData: " << myData(j,0) << ", ";
      int k = myData(j,0) - start;
//      Rcpp::Rcout << "  k: " << k << "\n";
      myMatrix( k, Rcpp::_ ) = myData( j - 0, Rcpp::_ );
      j++;
    }
    
//    Rcpp::Rcout << "  Found feature.\n";
//    Rcpp::Rcout << "  myData(j, 0): " << myData(j, 0) << ".\n";
//    Rcpp::Rcout << "  start: " << start << ".\n";
//    Rcpp::Rcout << "  end: " << end << ".\n";

    // Process feature.
//    while(myData(j,0) <= end + 0){
//      Rcpp::checkUserInterrupt();
//      Rcpp::Rcout << "j is: " << j;
//      Rcpp::Rcout << ", k is: " << k;
//      Rcpp::Rcout << ", start is " << start;
//      Rcpp::Rcout << "\n";
      
      // Handle when missing data at begining.
//      while( myData(j,0) > start + k ){
//        myMatrix(k, 0) = start + k;
//        k++;
//      }


//      j++;
//      k++;
//    }
    
    colnames(myMatrix) = Rcpp::CharacterVector::create("POS", "A", "C", "G", "T", "N", "*", "a", "c", "g", "t", "n");

    myList(i) = myMatrix;
    myNames(i) = myBed(i,3);
  }

  myList.attr("names") = myNames;

  return myList;
}


