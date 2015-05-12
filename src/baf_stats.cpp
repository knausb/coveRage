// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rcpp.h>

using namespace RcppParallel;
//using namespace Rcpp;


// http://gallery.rcpp.org/articles/parallel-distance-matrix/
// http://gallery.rcpp.org/articles/parallel-vector-sum/


struct bafstats_p : public Worker {
  
    // Input vector
    const RVector< std::string > invect;
//    const RVector< Rcpp::String > invect;
//    const RVector< double > invect;
    
    // Output matrix
    RMatrix< int > outmat;
    
    // initialize from Rcpp input and output matrixes (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
//    bafstats_p(const Rcpp::CharacterVector invect, Rcpp::IntegerMatrix outmat) : invect(invect), outmat(outmat) {}
    bafstats_p(const Rcpp::StringVector invect, Rcpp::IntegerMatrix outmat) : invect(invect), outmat(outmat) {}
//    bafstats_p(const Rcpp::CharacterVector invect) : invect(invect) {}
//    bafstats_p(Rcpp::IntegerMatrix outmat) : outmat(outmat) {}



};





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
