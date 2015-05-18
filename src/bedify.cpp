#include <Rcpp.h>
//using namespace Rcpp;


//' @title Parse data by a bed file
//' @rdname bedify
//' 
//' @description seperate a data matrix by using bed format data.
//' 
//' @param bed matrix of bed format data
//' @param mydata StringMatrix to be sorted
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
  
//  Rcpp::IntegerMatrix::Column starts = Rcpp::as< Rcpp::IntegerMatrix>(myBed(Rcpp::_, 1));
//  Rcpp::StringMatrix::Column starts = myBed(Rcpp::_, 1);
//  Rcpp::StringMatrix::Column ends = myBed(Rcpp::_, 2);

//  myList.attr("names") = Rcpp::CharacterVector::create();

  for(int i=0; i<myBed.nrow(); i++){
    Rcpp::Rcout << "Feature: " << i << "\n";
    
    // Count rows
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
//    Rcpp::Rcout << "nrows: " << nrows << "\n";
    
    Rcpp::IntegerMatrix myMatrix(nrows, myData.ncol());
    std::fill(myMatrix.begin(), myMatrix.end(), NA_INTEGER);
    
    // Increment to feature beginning.
    int j=0;
    while(myData(j,0) <= start){ j++;}
    Rcpp::Rcout << "  Found feature.\n";
    Rcpp::Rcout << "  myData(j, 0): " << myData(j, 0) << ".\n";
    Rcpp::Rcout << "  start: " << start << ".\n";
    Rcpp::Rcout << "  end: " << end << ".\n";

    // Process feature.
    while(myData(j,0) <= end + 1){
      Rcpp::Rcout << "j is: " << j;
      Rcpp::Rcout << ", start is " << start;
      Rcpp::Rcout << ", the difference is: " << j - start;
      Rcpp::Rcout << "\n";
      myMatrix( j - start, Rcpp::_ ) = myData( j - 1, Rcpp::_ );
      j++;
    }
    
    colnames(myMatrix) = Rcpp::CharacterVector::create("POS", "A", "C", "G", "T", "N", "*", "a", "c", "g", "t", "n");

    myList(i) = myMatrix;
    myNames(i) = myBed(i,3);
  }

  myList.attr("names") = myNames;

  return myList;
}


