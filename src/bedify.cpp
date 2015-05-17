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
    
    // Count rows
    std::string temp = Rcpp::as< std::string >(myBed(i,1));
    int start = stoi(temp);
    temp = Rcpp::as< std::string >(myBed(i,2));
    int end = stoi(temp);
    int nrows = end  - start + 1;
    Rcpp::Rcout << "nrows: " << nrows << "\n";
    
    Rcpp::IntegerMatrix myMatrix(nrows, myData.ncol());
    
    //- Rcpp::as<int>(myBed(i,2));
//    Rcpp::Rcout << "nrows: " << nrows << "\n";
//    Rcpp::IntegerMatrix myMatrix(nrows, 8);
//    Rcpp::IntegerMatrix myMatrix(nrows, myData.nrow());
//    Rcpp::IntegerMatrix myMatrix(2,2);
//    std::fill(myMatrix.begin(), myMatrix.end(), 2);
    std::fill(myMatrix.begin(), myMatrix.end(), NA_INTEGER);
    myList(i) = myMatrix;
    myNames(i) = myBed(i,3);
  }



//  Rcpp::StringMatrix data1(3,2);
//  Rcpp::StringMatrix data2(3,2);
//  Rcpp::List myList = Rcpp::List::create(Rcpp::Named("origDataFrame")=data1, Rcpp::Named("newDataFrame")=data2);

//  Rcpp::StringMatrix::Column myNames = myBed(Rcpp::_, 3);
//  Rcpp::StringVector myNames = myBed( Rcpp::_ , 2 );
  myList.attr("names") = myNames;

  return myList;
}


