#include <Rcpp.h>
//using namespace Rcpp;



std::vector< int > get_pos( Rcpp::StringMatrix myData ){
  std::vector< int > POS(myData.nrow());

  for(int i=0; i<POS.size(); i++){
    std::string temp = Rcpp::as< std::string >(myData(i,1));
    POS[i] = stoi(temp);
  }
  
  return(POS);
}



//myList(i) = proc_feature(myBed(i,0), myBed(i,1), myBed(i,2), POS, myData);
Rcpp::StringMatrix proc_feature( Rcpp::StringVector myBed,
                                 std::vector< int > POS,
                                 Rcpp::StringMatrix myData
                                 ){

  // Convert Rcpp::StringVector elements to int
  std::string temp = Rcpp::as< std::string >( myBed(1) );
  int start = stoi(temp);
  temp = Rcpp::as< std::string >( myBed(2) );
  int end = stoi(temp);

//  Rcpp::IntegerVector POSmyData( Rcpp::_, 1 )

  // Manage if feature is on reverse strand
  if(end < start){
    int tmp = start;
    start = end;
    end = tmp;
  }

  // Increment to chromosome.
  int i = 0; // Data row counter
  while ( myData(i,1) != myBed(0) ){
    i++;
  }

  // Increment to POS.
  while( POS[i] < start ){
    i++;
  }

  // Increment to the end of the feature
  int j=i;
  while( POS[j] <= end & j < POS.size() ){
    j++;
  }

  // We now have the information to declare a return matrix
  int nrows = j - i;
  Rcpp::StringMatrix myMatrix(nrows, myData.ncol());
  Rcpp::colnames(myMatrix) = Rcpp::colnames(myData);

  // Populate the return matrix
  for(int k = 0; k < myMatrix.nrow(); k++){
    Rcpp::checkUserInterrupt();
    myMatrix(k, Rcpp::_) = myData(k+i, Rcpp::_);
  }

  return myMatrix;
}




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
Rcpp::List bedify( Rcpp::StringMatrix myBed, Rcpp::StringMatrix myData ) {

  // Initialize return List
  Rcpp::List myList(myBed.nrow());
  myList.attr("names") = myBed( Rcpp::_ , 3 );
  
  // Col names for each return matrix
  Rcpp::StringVector myColNames = Rcpp::colnames(myData);

  // Convert POS to ints
  std::vector< int > POS = get_pos(myData);

  // Scroll over bed rows (features).
  for( int i=0; i<myBed.nrow(); i++ ){
    Rcpp::checkUserInterrupt();
    
    myList(i) = proc_feature(myBed(i,Rcpp::_), POS, myData);
//    myList(i) = proc_feature(myBed(i,0), myBed(i,1), myBed(i,2), POS, myData);
    
/*    
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
//    myNames(i) = myBed(i,3);
*/
  }

//  myList.attr("names") = myNames;

  return myList;
}








