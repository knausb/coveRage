// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rcpp.h>

using namespace RcppParallel;


//http://gallery.rcpp.org/articles/parallel-vector-sum/

struct Sum : public Worker
{   
   // source vector
   const RVector<double> input;
   
   // accumulated value
   double value;
   
   // constructors
   Sum(const Rcpp::NumericVector input) : input(input), value(0) {}
   Sum(const Sum& sum, Split) : input(sum.input), value(0) {}
   
   // accumulate just the element of the range I've been asked to
   void operator()(std::size_t begin, std::size_t end) {
      value += std::accumulate(input.begin() + begin, input.begin() + end, 0.0);
   }
     
   // join my value with that of another Sum
   void join(const Sum& rhs) { 
      value += rhs.value; 
   }
};


// [[Rcpp::export]]
double parallelVectorSum(Rcpp::NumericVector x) {
   
   // declare the SumBody instance 
   Sum sum(x);
   
   // call parallel_reduce to start the work
   parallelReduce(0, x.length(), sum);
   
   // return the computed sum
   return sum.value;
}


