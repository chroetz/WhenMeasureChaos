#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
int estimateRcppTaylor5(
    NumericMatrix u, NumericVector y0, NumericVector y1, double sd, int deltaI
) {
  double maxW = -1;
  int maxIdx;
  double w;
  int n = u.nrow(), d = u.ncol();
  double dst, dst1, dst2, dst3, dst4, dst5;
  double v;
  for (int i = 0; i < n - deltaI; ++i) {
    dst = 0;
    for (int j = 0; j < d; ++j) {
      v = u(i, j) - y0(j);
      dst += v*v;
    }
    dst1 = dst/sd/sd;
    dst2 = dst*dst;
    dst3 = dst2*dst;
    dst4 = dst3*dst;
    dst5 = dst4*dst;
    w = 1 - dst1 + dst2/2 - dst3/6 + dst4/24 - dst5/120; // approx exp(-dst/sd^2)
    
    if (w < maxW) continue;
    
    dst = 0;
    for (int j = 0; j < d; ++j) {
      v = u(i+deltaI, j) - y1(j);
      dst += v*v;
    }
    dst1 = dst/sd/sd;
    dst2 = dst*dst;
    dst3 = dst2*dst;
    dst4 = dst3*dst;
    dst5 = dst4*dst;
    w *= 1 - dst1 + dst2/2 - dst3/6 + dst4/24 - dst5/120;
    
    if (w > maxW) {
      maxIdx = i;
      maxW = w;
    }
  }
  return maxIdx;
}

// [[Rcpp::export]]
int estimateRcppExp(
    NumericMatrix u, NumericVector y0, NumericVector y1, double sd, int deltaI
) {
  double maxW = -1;
  int maxIdx;
  double w;
  int n = u.nrow(), d = u.ncol();
  double dst;
  double v;
  for (int i = 0; i < n - deltaI; ++i) {
    dst = 0;
    for (int j = 0; j < d; ++j) {
      v = u(i, j) - y0(j);
      dst += v*v;
    }
    w = exp(-dst/sd/sd);
    
    if (w < maxW) continue;
    
    dst = 0;
    for (int j = 0; j < d; ++j) {
      v = u(i+deltaI, j) - y1(j);
      dst += v*v;
    }
    w *= exp(-dst/sd/sd);
    
    if (w > maxW) {
      maxIdx = i;
      maxW = w;
    }
  }
  return maxIdx;
}

