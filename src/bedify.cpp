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
                                 Rcpp::StringMatrix myData,
                                 int fill_missing
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
  while ( myData(i,0) != myBed(0) & i < myData.nrow() ){
    i++;
  }
  
  // If we didn't find the chromosome return NA.
  if( i == myData.nrow() & fill_missing != 1 ){
    Rcpp::StringMatrix myMatrix(1);
    myMatrix(0,0) = NA_STRING;
    return myMatrix;
  }

  // Increment to POS.
  while( POS[i] < start & i < myData.nrow() ){
    i++;
  }
  // If we didn't find the POS return NA.
  if( i == myData.nrow() & fill_missing != 1 ){
    Rcpp::StringMatrix myMatrix(1);
    myMatrix(0,0) = NA_STRING;
    return myMatrix;
  }
  
  // Increment to the end of the feature
  int j=i;
  while( POS[j] <= end & j < myData.nrow() ){
    j++;
  }

//  for(int q=0; q<myData.nrow();q++){Rcpp::Rcout << myData(q,1)<<"\n";}

  // We now have the information to declare a return matrix
  // and populate it.
  if( fill_missing != 1 ){
    // Do not fill missubg data.
    Rcpp::StringMatrix myMatrix( j-i , myData.ncol());
    Rcpp::colnames(myMatrix) = Rcpp::colnames(myData);
    // Populate the return matrix
    for(int k = 0; k < myMatrix.nrow(); k++){
      Rcpp::checkUserInterrupt();
      myMatrix(k, Rcpp::_) = myData(k+i, Rcpp::_);
    }
    return myMatrix;
  } else {
    // Fill missing data.
    Rcpp::StringMatrix myMatrix( end - start + 1 , myData.ncol());
    Rcpp::colnames(myMatrix) = Rcpp::colnames(myData);

    // Populate the return matrix
    if( i >= myData.nrow()){
      // No data
      for(int k = 0; k < myMatrix.nrow(); k++){
        Rcpp::checkUserInterrupt();
        myMatrix(k,0) = myBed(0);
        myMatrix(k,1) = std::to_string(start + k);
        myMatrix(k,2) = NA_STRING;
//        for(int m=2; m<myMatrix.ncol(); m++){
//          myMatrix(k,m) = NA_STRING;
//        }
      }
    } else {
      // Data and possibly missing data    
      int l = 0;
      for(int k = 0; k < myMatrix.nrow(); k++){
        Rcpp::checkUserInterrupt();
        
        if( i + l < myData.nrow() ){
          // We have not overrun the file yet
          temp = Rcpp::as< std::string >( myData( i+l , 1 ) );
          int myPOS = stoi(temp);

          if( myPOS == start + k ){
//            myMatrix(k, Rcpp::_) = myData(k+i, Rcpp::_);
            myMatrix(k, Rcpp::_) = myData( i + l, Rcpp::_);

            l++;
          } else {
            myMatrix(k,0) = myBed(0);
            myMatrix(k,1) = std::to_string(myPOS);
            myMatrix(k,2) = NA_STRING;
            //myMatrix(k,1) = myBed(1) + k;
            //for(int m=2; m<myMatrix.ncol(); m++){
            //  myMatrix(k,m) = NA_STRING;
            //}
          }
        } else {
          // We've overrun the rows in the file.
          myMatrix(k,0) = myBed(0);
          myMatrix(k,1) = std::to_string( start + k );
          myMatrix(k,2) = NA_STRING;
          //myMatrix(k,1) = myBed(1) + k;
          //for(int m=2; m<myMatrix.ncol(); m++){
          //  myMatrix(k,m) = NA_STRING;
          //}
        }
      }
    }
    return myMatrix;
  }

  Rcpp::Rcerr << "You should never get here, something bad has happened!\n";
}




//' @title Parse data by a bed file
//' @rdname bedify
//' @name bedify
//' 
//' @description Seperate a data matrix into list elements based on coordinates from bed format data.
//' 
//' @param myBed matrix of bed format data
//' @param myData StringMatrix or IntegerMatrix to be sorted
//' @param fill_missing include records for when there is no data (0, 1).  By default these records are omitted.
//' @param verbose should verbose output be generated (0, 1)
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
Rcpp::List bedify( Rcpp::StringMatrix myBed,
                   Rcpp::StringMatrix myData, 
                   int fill_missing = 0, 
                   int verbose = 0 ) {

  // Start a timer
  time_t result = time(nullptr);

  // Initialize return List
  Rcpp::List myList(myBed.nrow());

  Rcpp::StringVector myListNames( myBed.nrow() );
  for(int i = 0; i < myBed.nrow(); i++){
    myListNames(i) = myBed( i, 3 );
  }
  myList.attr( "names" ) = myListNames;
  
  // Col names for each return matrix
  Rcpp::StringVector myColNames = Rcpp::colnames(myData);

  // Convert POS to ints
  std::vector< int > POS = get_pos(myData);

  // Scroll over bed rows (features).
  for( int i=0; i<myBed.nrow(); i++ ){
    Rcpp::checkUserInterrupt();

    if( verbose == 1){
      Rcpp::Rcout << "Searching for feature " << i + 1 << ": " << myBed(i,3);
      Rcpp::Rcout << " on " << myBed(i,0);
      Rcpp::Rcout << " at " << time(nullptr) - result << " seconds.\n";
    }
    
    myList(i) = proc_feature( myBed(i,Rcpp::_), POS, myData, fill_missing );
//    myList(i) = proc_feature(myBed(i,0), myBed(i,1), myBed(i,2), POS, myData);
    

  }

//  myList.attr("names") = myNames;

  return myList;
}








