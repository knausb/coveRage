// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rcpp.h>

//using namespace RcppParallel;
//using namespace Rcpp;

// http://gallery.rcpp.org/articles/parallel-distance-matrix/
// http://gallery.rcpp.org/articles/parallel-vector-sum/

// http://www.cplusplus.com/doc/tutorial/classes/

struct bafstats_p : public RcppParallel::Worker {
  
  // Input matrix
  const RcppParallel::RMatrix< Rcpp::String > inmat1;
//  const RcppParallel::RMatrix< std::string > inmat1;
  
  // Output matrix
  RcppParallel::RMatrix< int > outmat;
  
  // Threshold
  int qmin;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  bafstats_p(const Rcpp::StringMatrix inmat2, Rcpp::IntegerMatrix outmat, int qmin)
    : inmat1(inmat1), outmat(outmat), qmin(qmin) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      
      
//      RcppParallel::RMatrix< Rcpp::String >::Row row2 = inmat1.row(i);
//      Rcpp::Rcout << row2[0] << "\n";
//      outmat(i,0)++;
      Rcpp::String calls = &inmat1(i,1);
//      for(int j=0; j < calls.length(); j++){

//      }
//      std::string calls = inmat1(i,1);
//      for(int j=0; j < inmat1(i,1).length(); j++){
        
//      }

      // Input data
//      Rcpp::Rcout << inmat1(i,0) << "\n";
//      Rcpp::String ref   = inmat1(i,0);
//      std::string ref   = inmat1(i,0);
//      std::string calls = inmat1(i,1);
//      std::string quals = inmat1(i,2);
      
      // row we will operate on
//      RcppParallel::RMatrix< int >::Row row2 = outmat.row(i);
      
      
      

    }
  }
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
//' @details
//' The reference alleles must be in all upper case.
//' See \code{toupper} if they are not.
//'
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix baf_stats(Rcpp::StringVector calls,
                              Rcpp::StringVector quals,
                              Rcpp::StringVector ref,
                              int minq = 0) {
//   Rcpp::IntegerMatrix stats(calls.size(), 11);
//   Rcpp::IntegerMatrix stats(3, 3);
//   Rcpp::IntegerMatrix stats(3, 11);


    // Agregate the input matrix
    Rcpp::StringMatrix inmat(calls.size(), 3);
    inmat(Rcpp::_, 0) = ref;
    inmat(Rcpp::_, 1) = calls;
    inmat(Rcpp::_, 2) = quals;
    
    // allocate the matrix we will return
    Rcpp::IntegerMatrix outmat(calls.size(), 11);
    colnames(outmat) = Rcpp::CharacterVector::create("A", "C", "G", "T", "N", "*", "a", "c", "g", "t", "n");

   // create the worker
//   bafstats_p bafstats_p(inmat, outmat);
   bafstats_p bafstats_p(inmat, outmat, minq);
     
   // call it with parallelFor
   parallelFor(0, inmat.nrow(), bafstats_p);


//   return stats;
   return outmat;
}


//' @name baf_stats_st
//' @title baf_stats_st
//' @rdname baf_stats
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix baf_stats_st(Rcpp::StringVector calls,
                              Rcpp::StringVector quals,
                              Rcpp::StringVector ref,
                              int minq = 0) {

  // allocate the matrix we will return
  Rcpp::IntegerMatrix outmat(calls.size(), 11);
  colnames(outmat) = Rcpp::CharacterVector::create("A", "C", "G", "T", "N", "*", "a", "c", "g", "t", "n");

  // Agregate the input matrix
  Rcpp::StringMatrix inmat(calls.size(), 3);
  inmat(Rcpp::_, 0) = ref;
  inmat(Rcpp::_, 1) = calls;
  inmat(Rcpp::_, 2) = quals;
  
  for(int i=0; i<inmat.nrow(); i++){
    std::string ref(inmat(i,0));
    std::string calls(inmat(i,1));
    std::string quals(inmat(i,2));
//    Rcpp::Rcout << calls << "\n";

//    Rcpp::Rcout << calls.size() << ": " <<calls << " " << quals.size() << ": " << quals << "\n";
    for(int j=0; j<calls.size(); j++){
      
      if(calls[j] == '$'){
        // Read ends
        calls.erase(j, 1);
      }
      if(calls[j] == '^'){
        // Read beginnings include ascii mapping quality
        calls.erase(j, 2);
      }
      
      int qual = static_cast<int>(quals[j]);
      qual = qual - 33; // Sanger format.
      if( qual >= minq ){

          if(calls[i] == '.'){
            if(ref == "A"){
              outmat(i,0)++;
            } else if(ref == "C"){
              outmat(i,1)++;
            } else if(ref == "G"){
              outmat(i,2)++;
            } else if(ref == "T"){
              outmat(i,3)++;
            } else if(ref == "N"){
              outmat(i,4)++;
            }
          } else if (calls[i] == ','){
            if(ref == "A"){
              outmat(i,6)++;
            } else if(ref == "C"){
              outmat(i,7)++;
            } else if(ref == "G"){
              outmat(i,8)++;
            } else if(ref == "T"){
              outmat(i,9)++;
            } else if(ref == "N"){
              outmat(i,10)++;
            }
          } else if (calls[i] == 'A'){
            outmat(i,0)++;
          } else if (calls[i] == 'C'){
            outmat(i,1)++;
          } else if (calls[i] == 'G'){
            outmat(i,2)++;
          } else if (calls[i] == 'T'){
            outmat(i,3)++;
          } else if(calls[i] == 'N'){
            outmat(i,4)++;
          } else if(calls[i] == '*'){
            outmat(i,5)++;
          } else if(calls[i] == 'a'){
            outmat(i,6)++;
          } else if(calls[i] == 'c'){
            outmat(i,7)++;
          } else if(calls[i] == 'g'){
            outmat(i,8)++;
          } else if(calls[i] == 't'){
            outmat(i,9)++;
          } else if(calls[i] == 'n'){
            outmat(i,10)++;
          }
      }
    }


  }

  return outmat;
}