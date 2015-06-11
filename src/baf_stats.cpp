// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rcpp.h>

//using namespace RcppParallel;
//using namespace Rcpp;


void call_pu(Rcpp::IntegerVector counts, char call, std::string ref){
  
  if(call == '.'){
    if(ref == "A"){
      counts(0)++;
    } else if(ref == "C"){
      counts(1)++;
    } else if(ref == "G"){
      counts(2)++;
    } else if(ref == "T"){
      counts(3)++;
    } else if(ref == "N"){
      counts(4)++;
    }

  } else if (call == ','){
    if(ref == "A"){
      counts(6)++;
    } else if(ref == "C"){
      counts(7)++;
    } else if(ref == "G"){
      counts(8)++;
    } else if(ref == "T"){
      counts(9)++;
    } else if(ref == "N"){
      counts(10)++;
    }

  } else if (call == 'A'){
    counts(0)++;
  } else if (call == 'C'){
    counts(1)++;
  } else if (call == 'G'){
    counts(2)++;
  } else if (call == 'T'){
    counts(3)++;
  } else if(call == 'N'){
    counts(4)++;
  } else if(call == '*'){
    counts(5)++;
  } else if(call == 'a'){
    counts(6)++;
  } else if(call == 'c'){
    counts(7)++;
  } else if(call == 'g'){
    counts(8)++;
  } else if(call == 't'){
    counts(9)++;
  } else if(call == 'n'){
    counts(10)++;
  }
  
}





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
//' @param inMatrix input matrix
//' 
//' @details The character matrix \strong{inMatrix} should consist of columns  5 columns. 
//' The first column is the chromosome name (and is not presently used).
//' The second column is the chromosomal position.
//' The third column is the reference allele.
//' The fourth column is a string of calls.
//' The fifth column is a string of qualities.
//' This is expected to come from mpileup output.
//' Note that while mpileup can include data for multiple samples, here we need to process each sample seperately.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix baf_stats_st(Rcpp::StringMatrix inMatrix,
                                 int minq = 0) {

  // allocate the matrix we will return
  Rcpp::IntegerMatrix outMatrix(inMatrix.nrow(), 12);
  colnames(outMatrix) = Rcpp::CharacterVector::create("POS", "A", "C", "G", "T", "N", "*", "a", "c", "g", "t", "n");
  
  // Iterate over rows (sites or POSitions)
  for(int i=0; i<inMatrix.nrow(); i++){

    std::string ref   = Rcpp::as< std::string >(inMatrix(i,2));
    std::string calls = Rcpp::as< std::string >(inMatrix(i,3));
    std::string quals = Rcpp::as< std::string >(inMatrix(i,4));
    
    Rcpp::IntegerVector out_row = outMatrix(i,Rcpp::_);
    out_row.erase(0);  // Remove POS

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
      qual = qual - 33; // Sanger encoding.
      if( qual >= minq ){
        call_pu(out_row, calls[j], ref);
      }
    }
    
    // Load counts into output matrix
    std::string strPOS = Rcpp::as< std::string >(inMatrix(i,1));
    outMatrix(i,0) = std::stoi( strPOS );
    for(int j=0; j<out_row.size(); j++){
      outMatrix(i,j+1) = out_row(j);
    }

  }
  return outMatrix;
}

