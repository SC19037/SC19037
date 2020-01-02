#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' @title Cpp function for Exercise 9.4
//' @description Cpp function for Exercise 9.4
//' @param sigma the variance of the normality distribution
//' @param x_0 the original value
//' @param N the size
//' @return the sample x \code{n}
//' @export
// [[Rcpp::export]]
NumericMatrix rwC(double sigma, double x_0,int N) {
  NumericMatrix x(N,1);
  NumericVector U=runif(N,0,1);
  x(0,0)=x_0; 
  double y;double a;double b;
  for (int i=1; i<N; i++){
    y = rnorm(1,x(i-1,0),sigma)[0];
    a = 0.5 * exp(-abs(y));
    b = 0.5 * exp(-abs(x(i-1,0)));
    if (U[i] <= a/b){
      x(i,0) = y;
    }
    else{
      x(i,0)=x(i-1,0);
    }
  }
  return x;
}



