

#include <Rcpp.h>

//using namespace Rcpp;
//using namespace RcppParallel;

//' @name baf_stats
//' @title baf_stats
//' @rdname baf_stats
//' 
//' @param calls vector of pileup calls
//' @param quals vector of pileup calls
//' @param ref vector of reference alleles
//' @param minq minimum quality for call to be retained
//'
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix baf_stats(Rcpp::CharacterVector calls,
                              Rcpp::CharacterVector quals,
                              Rcpp::CharacterVector ref,
                              int minq) {
   Rcpp::IntegerMatrix stats(calls.size(), 11);
//   Rcpp::IntegerMatrix stats(3, 3);
//   Rcpp::IntegerMatrix stats(3, 11);
   
//   colnames(stats) = Rcpp::CharacterVector::create("A", "B", "C");
//   Rcpp::colnames(stats) = Rcpp::CharacterVector::create("A", "B", "C");
//   Rcpp::StringVector colnames(11) ;
//   stats.attr("colnames") = Rcpp::CharacterVector::create('A', 'C', 'G', 'T', 'N', '*', 'a', 'c', 'g', 't', 'n');
   colnames(stats) = Rcpp::CharacterVector::create("A", "C", "G", "T", "N", "*", "a", "c", "g", "t", "n");

   return stats;
}
