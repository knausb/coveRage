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
    
    // Initialize return matrix.
    Rcpp::IntegerMatrix myMatrix(nrows, myData.ncol());
    std::fill(myMatrix.begin(), myMatrix.end(), NA_INTEGER);
    for(int j=0; j<myMatrix.nrow(); j++){
      myMatrix(j,0) = start + j;
    }

    int j=0;  // matrix to bedify row counter
    
    // Increment to feature beginning.
    while(myData(j,0) < start){
      Rcpp::checkUserInterrupt();
      j++;
    }
    
    // Process to the end.
    while( myData(j,0) <= end & j < myData.nrow() ){
      Rcpp::checkUserInterrupt();
      int k = myData(j,0) - start;
      myMatrix( k, Rcpp::_ ) = myData( j - 0, Rcpp::_ );
      j++;
    }
    
    colnames(myMatrix) = Rcpp::CharacterVector::create("POS", "A", "C", "G", "T", "N", "*", "a", "c", "g", "t", "n");

    myList(i) = myMatrix;
    myNames(i) = myBed(i,3);
  }

  myList.attr("names") = myNames;

  return myList;
}




std::vector< int > get_pos( Rcpp::StringMatrix myData ){
  std::vector< int > POS(myData.nrow());

  for(int i=0; i<POS.size(); i++){
    std::string temp = Rcpp::as< std::string >(myData(i,1));
    POS[i] = stoi(temp);
  }
  
  return(POS);
}




//' @title Parse data by a bed file
//' @rdname bedify
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List bedify_sm( Rcpp::StringMatrix myBed, Rcpp::StringMatrix myData ) {

  // Initialize return datastructure
  Rcpp::List myList(myBed.nrow());
  
  // Convert POS to ints
  std::vector< int > POS = get_pos(myData);

//  for(int i=0; i < POS.size(); i++){Rcpp::Rcout << "POS: " << POS[i] << "\t";}
//  Rcpp::Rcout << "POS size: " << POS.size() << "\t";

  // To collect feature names and name list before return
  Rcpp::StringVector myNames(myBed.nrow());

  // Scroll over bed rows (features).
  for(int i=0; i<myBed.nrow(); i++){
    Rcpp::checkUserInterrupt();
    
    // Extract StringMatrix elements and convert to int
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

    int j=0;  // matrix to bedify row counter



    // Increment to feature beginning.
    while(POS[j] < start){
//      Rcpp::Rcout << "POS: " << POS[j] << "\n";
      Rcpp::checkUserInterrupt();
      j++;
    }
//    Rcpp::Rcout << "Feature begins at: " << j << "\n";



    // Increment past all records for the feature
    int k=j;
    while(POS[k] <= end){
//      Rcpp::Rcout << "POS: " << POS[k] << "\n";
      Rcpp::checkUserInterrupt();
      k++;
    }
//    Rcpp::Rcout << "Feature ends at: " << k << "\n";
    
//    Rcpp::Rcout << "Made it.\n";
    
    // Declare and populate a return matrix
    int nrows = k - j + 0;
    Rcpp::StringMatrix myMatrix(nrows, myData.ncol());
    for(int l=j; l<k; l++){
      Rcpp::checkUserInterrupt();
      myMatrix(l-j, Rcpp::_) = myData(l, Rcpp::_);
    }
    
    myList(i) = myMatrix;
    myNames(i) = myBed(i,3);
  }

  myList.attr("names") = myNames;
  return myList;
}


