#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

//' Empirical Cumulative Distribution Function (ECDF)
//'
//' @description This function calculates the jump points of the empirical
//' cumulative distribution function. It identifies unique sorted values and
//' their corresponding cumulative probabilities.
//'
//' @param x A numeric vector of observations.
//' @return A numeric vector of length 2 * J (where J is the number of unique
//' values in x). The first half of the vector contains the unique sorted
//' values, and the second half contains their cumulative probabilities
//' ranging from (1/n) to 1.0.
//'
//' @details
//' The function sorts the input vector and iterates through it to identify
//' unique values. The probability for each unique value \eqn{x_j} is
//' calculated as \eqn{F(x_j) = i/n}, where \eqn{i} is the count of all
//' elements \eqn{\le x_j}.
//'
//' @examples
//' \dontrun{
//' x <- c(1, 2, 2, 3, 5)
//' res <- c_ecdf(x)
//' n_unique <- length(res) / 2
//' values <- res[1:n_unique]
//' probs <- res[(n_unique + 1):length(res)]
//' }
//' @export
//'
// [[Rcpp::export]]
NumericVector c_ecdf(NumericVector x){

  // 1. Filter out NAs first
  // LogicalVector is_na = is_na(x); // Rcpp way
  // However, for speed, we can do a single pass:
  std::vector<double> clean_x;
  clean_x.reserve(x.size());

  for(int i = 0; i < x.size(); ++i) {
    if (!NumericVector::is_na(x[i])) {
      clean_x.push_back(x[i]);
    }
  }

  int n = clean_x.size();
  if (n == 0) return NumericVector(0);

  // 1. Sort the input
  std::sort(clean_x.begin(),clean_x.end());

  // 2. Identify unique values and counts
  // We'll store unique values and the cumulative index
  NumericVector xu(n);
  NumericVector wu(n);


  int j = 0;
  for (int i = 0; i < n; ++i) {
    // If it's the last occurrence of a specific value
    if (i == n - 1 || clean_x[i] != clean_x[i+1]) {
      xu[j] = clean_x[i];
      // Cumulative probability is (index + 1) / total N
      wu[j] = static_cast<double>(i + 1) / n;
      j++;
    }
  }

  // 3. Generate compact output
  NumericVector out(2 * j);
  for (int i = 0; i < j; ++i) {
    out[i] = xu[i];         // Values
    out[i + j] = wu[i];     // Probabilities
  }

  //std::cout << j << std::endl;

  return out;

}
