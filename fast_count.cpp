#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector fast_count(NumericVector vals, NumericVector unique_vals, NumericVector counts) {
  IntegerVector out;
  // for (NumericVector::iterator iter = vals.begin(); iter != vals.end(); ++iter) {
  //   NumericVector::iterator it =
  //     std::find(unique_vals.begin(), unique_vals.end(), *iter);
  //   out.push_back(*(counts.begin() + (it - unique_vals.begin())));
  // }
  for (int i = 0; i < vals.size(); ++i) {
    for (int j = 0; j < unique_vals.size(); ++j) {
      if (vals[i] == unique_vals[j]) {
        out.push_back(counts[j]);
      }
    }
  }
  return out;
}

/*** R
# (tmp <- sample(1:5, 20, replace = TRUE))
# (tab <- table(tmp))
# fast_count(tmp, as.numeric(names(tab)), tab)
*/
