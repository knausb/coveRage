#include <Rcpp.h>
//using namespace Rcpp;


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
   
   Rcpp::StringVector colnames(11) ;
//  colnames(0) = "Primero";
//  colnames(1) = "Segundo";
//  colnames(2) = "Tercero";
//  colnames(3) = "Cuarto";
//   stats.attr("colnames") = Rcpp::CharacterVector::create('A', 'C', 'G', 'T', 'N', '*', 'a', 'c', 'g', 't', 'n');
   colnames(stats) = Rcpp::CharacterVector::create("A", "C", "G", "T", "N", "D", "a", "c", "g", "t", "n");

   return stats;
}
