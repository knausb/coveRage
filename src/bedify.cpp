#include <Rcpp.h>
//using namespace Rcpp;


//' @title Parse data by a bed file
//' @rdname bedify
//' @name bedify
//' 
//' @description Seperate a data matrix into list elements based on coordinates from bed format data.
//' 
//' @param myBed matrix of bed format data
//' @param myData StringMatrix or IntegerMatrix to be sorted
//' 
//' @details
//' 
//' \strong{Bed format} data contain at least four columns.
//' The first column indicates the chromosome (i.e., supercontig, scaffold, contig, etc.).
//' The second cotains the starting positions.
//' The third the ending positions.
//' The fourth are the names of the features.
//' All subsequent columns are ignored here.
//' In an attempt to optimize performance the data are expected to be formatted as a character matrix.
//' The starting and end positions are converted to numerics internally.
//' 
//' The \strong{matrix format} used here is based on vcf type data.
//' Typically these data have a chromosome as the first column.
//' Each chromosome has its own coordinate system which begins at one.
//' This means that using multiple chromosomes will necessitate some fix to the coordinate systems.
//' Here I take the perspective that you should simply work on one chromosome at a time, so the chromosome information is ignored.
//' The first column is the chromosome, which I ignore.
//' The second column is the position, which is used for sorting.
//' Subsequent columns are not treated but are brought along with the subset.
//' 
//' 
//' When the matrix is of numeric form the first column, which contains the chromosome identifier (CHROM), must also be numeric.
//' This is because matrix elements must all be of the same type.
//' 
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




// @title Parse data by a bed file
// @rdname bedify
// 
//
// @export
// [[Rcpp::export]]




Rcpp::List bedify_sm( Rcpp::StringMatrix myBed, Rcpp::StringMatrix myData ) {

  // Initialize return datastructure
  Rcpp::List myList(myBed.nrow());
  Rcpp::StringVector myColNames = Rcpp::colnames(myData);
  
  // Convert POS to ints
  std::vector< int > POS = get_pos(myData);

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
      Rcpp::checkUserInterrupt();
      j++;
    }

    // Increment past all records for the feature
    int k=j;
    while( POS[k] <= end & k < POS.size() ){
      k++;
    }

    // Declare and populate a return matrix
    int nrows = k - j + 0;
    Rcpp::StringMatrix myMatrix(nrows, myData.ncol());
    
    for(int l = 0; l < myMatrix.nrow(); l++){
      Rcpp::checkUserInterrupt();
      myMatrix(l, Rcpp::_) = myData(l+j, Rcpp::_);
    }
    
    if( myColNames.size() > 0){
      Rcpp::colnames(myMatrix) = myColNames;
    }
    
    myList(i) = myMatrix;
    myNames(i) = myBed(i,3);
  }

  myList.attr("names") = myNames;
  return myList;
}


//' @title Parse data by a bed file
//' @rdname bedify
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List bedify_nm( Rcpp::StringMatrix myBed, Rcpp::NumericMatrix myData ) {

  // Initialize return datastructure
  Rcpp::List myList(myBed.nrow());
  Rcpp::StringVector myColNames = Rcpp::colnames(myData);
  
  // Convert POS to ints
//  std::vector< int > POS = get_pos(myData);
  Rcpp::NumericVector POS = myData(Rcpp::_,1);

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
      Rcpp::checkUserInterrupt();
      j++;
    }

    // Increment past all records inclusive of the feature
    int k=j;
    while( POS[k] <= end & k < POS.size() ){
      k++;
    }

    // Declare and populate a return matrix
    int nrows = k - j + 0;
    Rcpp::NumericMatrix myMatrix(nrows, myData.ncol());
    

    
    for(int l = 0; l < myMatrix.nrow(); l++){
      Rcpp::checkUserInterrupt();
      myMatrix(l, Rcpp::_) = myData(l+j, Rcpp::_);
    }
    
    if( myColNames.size() > 0){
      Rcpp::colnames(myMatrix) = myColNames;
    }

    myList(i) = myMatrix;
    myNames(i) = myBed(i,3);
  }

  myList.attr("names") = myNames;
  return myList;
}




