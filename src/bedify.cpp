#include <Rcpp.h>
//using namespace Rcpp;


//' @title Parse data by a bed file
//' @rdname bedify
//' 
//' @description seperate a data matrix by using bed format data.
//' 
//' @param bed data.frame
//' @param mydata StringMatrix to be sorted
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List bedify( Rcpp::DataFrame bed, Rcpp::StringMatrix mydata ) {
  
  Rcpp::StringMatrix data1(3,2);
  Rcpp::StringMatrix data2(3,2);
  
  Rcpp::List myList = Rcpp::List::create(Rcpp::Named("origDataFrame")=data1, Rcpp::Named("newDataFrame")=data2);
  
  return myList;
}


