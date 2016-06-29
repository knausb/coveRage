#include <Rcpp.h>
#include <time.h>
//using namespace Rcpp;



std::vector< int > get_pos( Rcpp::StringMatrix myData ){
  std::vector< int > POS(myData.nrow());
  
  for(int i=0; i<POS.size(); i++){
    std::string temp = Rcpp::as< std::string >(myData(i,1));
    POS[i] = atoi(temp.c_str());
//    POS[i] = stoi(temp);
  }
  
  return(POS);
}



//myList(i) = proc_feature(myBed(i,0), myBed(i,1), myBed(i,2), POS, myData);
Rcpp::StringMatrix proc_feature( Rcpp::StringVector myBed,
                                 std::vector< int > POS,
                                 Rcpp::StringMatrix myData,
                                 int fill_missing
                                 ){
  
  // Create an empty matrix to return in exceptions.
  Rcpp::StringMatrix MT_matrix(0, myData.ncol());

  // Rcpp::StringVector myBed includes start and stop integer coordinates.
  // Convert Rcpp::StringVector elements to int
  std::string temp = Rcpp::as< std::string >( myBed(1) );
  int start = atoi(temp.c_str());
  temp = Rcpp::as< std::string >( myBed(2) );
  int end = atoi(temp.c_str());

  // Manage if feature is on reverse strand
  if(end < start){
    int tmp = start;
    start = end;
    end = tmp;
  }

  // Increment i so that chromosome in myData
  // matches the chromosome in the single BED record.
  int i = 0; // Data row counter
  while ( myData(i,0) != myBed(0) && i < myData.nrow() ){
    i++;
  }
  
  // If we didn't find the chromosome return an empty matrix.
  if( i == myData.nrow() & fill_missing != 1 ){
//    Rcpp::StringMatrix myMatrix(0, myData.ncol());
//    Rcpp::StringMatrix myMatrix(1);
//    myMatrix(0,0) = NA_STRING;
//    myMatrix( 0, myData.ncol() );
//    return myMatrix;
    return MT_matrix;
  }

  // We should now have i at the correct chromosome in myData.
  // POS is the integer recast of POS in myData.
  // We can now increment to the correct position in the chromosome
  // by incrementing to the start of teh annotation.
  //
  // Increment to POS.
  while( myData(i,0) == myBed(0) && POS[i] < start && i < myData.nrow() ){
    i++;
  }
  // If we didn't find the POS return an empty matrix.
  if( i == myData.nrow() & fill_missing != 1 ){
//    Rcpp::StringMatrix myMatrix(0, myData.ncol());
//    myMatrix(0,0) = NA_STRING;
//    return myMatrix;
    return MT_matrix;
  }
  
  // Increment to the end of the feature
  int j=i;
  while( myData(i,0) == myBed(0) && POS[j] <= end && j < myData.nrow() ){
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
//        myMatrix(k,1) = std::to_string(start + k);
        std::ostringstream stm;
        stm << start + k;
        myMatrix(k,1) = stm.str();
        
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
          int myPOS = atoi(temp.c_str());
//          int myPOS = stoi(temp);

          if( myPOS == start + k ){
//            myMatrix(k, Rcpp::_) = myData(k+i, Rcpp::_);
            myMatrix(k, Rcpp::_) = myData( i + l, Rcpp::_);

            l++;
          } else {
            myMatrix(k,0) = myBed(0);
//            myMatrix(k,1) = std::to_string(myPOS);
            std::ostringstream stm;
            stm << myPOS;
            myMatrix(k,1) = stm.str();
            
            myMatrix(k,2) = NA_STRING;
            //myMatrix(k,1) = myBed(1) + k;
            //for(int m=2; m<myMatrix.ncol(); m++){
            //  myMatrix(k,m) = NA_STRING;
            //}
          }
        } else {
          // We've overrun the rows in the file.
          myMatrix(k,0) = myBed(0);
//          myMatrix(k,1) = std::to_string( start + k );
          std::ostringstream stm;
          stm << start + k;
          myMatrix(k,1) = stm.str();

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
//' \strong{Bed format} data contain at least three columns.
//' The first column indicates the chromosome (i.e., supercontig, scaffold, contig, etc.).
//' The second cotains the starting positions.
//' The third the ending positions.
//' Optional columns are in columns four through nine.
//' For example, the fourth column may contain the names of features.
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
//' \href{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}{Bed format} at UCSC
//' 
//' 
//' @examples
//' 
//' bed <- structure(c("chr_290", "chr_4176", "chr_126921", "chr_126921", 
//' "chr_125157", "chr_125157", "chr_125157", "chr_125157", "chr_126888", 
//' "chr_126888", "47", "400", "4344", "1", "3712", "6025", "2269", 
//' "1779", "7930", "4637", "80", "500", "4967", "9066", "6566", 
//' "6450", "2933", "2226", "11939", "7913", "gene_1", "gene_2", 
//' "gene_3", "gene_4", "gene_5", "gene_6", "gene_7", "gene_8", "gene_9", 
//' "gene_10"), .Dim = c(10L, 4L), .Dimnames = list(NULL, c("chrom", 
//' "chromStart", "chromEnd", "name")))
//' 
//' 
//' vcf.matrix <- structure(c("chr_290", "chr_290", "chr_4176", "chr_4176", "chr_50514", 
//' "chr_64513", "chr_107521", "chr_121987", "chr_122006", "chr_122006", 
//' "78", "96", "406", "425", "863", "2853", "77", "103", "243", 
//' "636", "0/1:5,4:9:99:117,0,153", "0/0:9,0:9:99:0,27,255", "0/1:10,11:21:99:255,0,255", 
//' "0/1:10,11:21:99:255,0,255", "0/1:14,14:28:99:255,0,255", "0/1:29,13:42:99:255,0,255", 
//' "0/1:26,11:37:99:255,0,255", "0/1:21,14:35:99:255,0,255", "0/0:12,1:13:67:0,4,255", 
//' "0/1:55,8:63:99:99,0,255", "0/1:10,8:18:99:234,0,255", "0/0:17,0:17:99:0,51,255", 
//' "0/1:16,13:29:99:255,0,255", "0/1:16,13:29:99:255,0,255", "0/1:26,19:45:99:255,0,255", 
//' "0/1:50,19:69:99:255,0,255", "0/1:62,17:79:99:255,0,255", "0/1:95,22:117:99:255,0,255", 
//' "0/1:32,5:37:99:68,0,255", "0/1:69,21:90:99:255,0,255"), .Dim = c(10L, 
//' 4L), .Dimnames = list(NULL, c("CHROM", "POS", "sample_1", "sample_2"
//' )))
//' 
//' 
//' class(bed)
//' is.character(bed)
//' class(vcf.matrix)
//' is.character(vcf.matrix)
//' 
//' var.list <- bedify(bed, vcf.matrix)
//' table(unlist(lapply(var.list, nrow)))
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

  // Col names for each return matrix
  Rcpp::StringVector myColNames = Rcpp::colnames(myData);
  
  Rcpp::StringVector myListNames( myBed.nrow() );
  for(int i = 0; i < myBed.nrow(); i++){
    myListNames(i) = myBed( i, 3 );
  }
  myList.attr( "names" ) = myListNames;


  // Convert POS to ints
  std::vector< int > POS = get_pos(myData);

  
  // Begin to parse the data.
  //
  // Scroll over bed rows (features).
  for( int i=0; i<myBed.nrow(); i++ ){
    Rcpp::checkUserInterrupt();

    if( verbose == 1){
      Rcpp::Rcout << "Searching for feature " << i + 1 << ": " << myBed(i,3);
      Rcpp::Rcout << " on " << myBed(i,0);
      Rcpp::Rcout << " at " << time(nullptr) - result << " seconds.\n";
    }
    
    // Send one annotation (bed row) to proc_feature.
    myList(i) = proc_feature( myBed(i,Rcpp::_), POS, myData, fill_missing );
//    myList(i) = proc_feature(myBed(i,0), myBed(i,1), myBed(i,2), POS, myData);

    
    Rcpp::colnames(myList(i)) = myColNames;

  }

//  myList.attr("names") = myNames;

  return myList;
}


// EOF.