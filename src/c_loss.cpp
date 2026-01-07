#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// Forward Declaration
NumericVector c_ecdf_values(NumericVector uni, NumericVector w, NumericVector z, double ofset = 0);

//' Probability of an Unfavorable Pair (Loss Probability)
//'
//' @description
//' Calculates the probability of an unfavorable pair, defined as \eqn{P(X < Y - \lambda)},
//' across a vector of thresholds \eqn{\lambda}. This represents the probability
//' that group Y outperforms group X by at least the margin \eqn{\lambda}.
//'
//' @param x_uni A sorted numeric vector of unique event times for group X.
//' @param x_w Cumulative probabilities (ECDF/Kaplan-Meier) for group X.
//' @param y_uni A sorted numeric vector of unique event times for group Y.
//' @param y_w Cumulative probabilities (ECDF/Kaplan-Meier) for group Y.
//' @param lambda A numeric vector of thresholds (\eqn{\lambda}).
//'
//' @return A numeric vector representing \eqn{P(X < Y - \lambda)} for each
//' threshold provided.
//'
//' @details
//' The function calculates the "Loss Probability" by integrating the probability
//' that Y exceeds \eqn{X + \lambda} against the probability mass of X:
//' \deqn{P(X < Y - \lambda) = \sum_{i} [1 - F_Y(x_i + \lambda)] \cdot f_X(x_i)}
//' where \eqn{F_Y} is the right-continuous ECDF of group Y.
//'
//' @export
// [[Rcpp::export]]
 NumericVector c_loss(NumericVector x_uni, NumericVector x_w,
                     NumericVector y_uni, NumericVector y_w,
                     NumericVector lambda) {

   int n_lam = lambda.size();
   int nx_u = x_uni.size();
   NumericVector out(n_lam);

   if (nx_u == 0) return out;

   // 1. Pre-calculate the 'mass' (jump heights) of X once
   // This avoids re-calculating (x_w[j] - x_w[j-1]) inside the lambda loop
   NumericVector ft(nx_u);
   ft[0] = x_w[0];
   for (int j = 1; j < nx_u; ++j) {
     ft[j] = x_w[j] - x_w[j - 1];
   }

   // 2. Loop through each lambda
   for (int i = 0; i < n_lam; ++i) {
     double current_sum = 0.0;

     // Evaluate Fc for all X points using the current lambda
     // Returns P(Y < x_uni - lambda[i])
     NumericVector Fc = c_ecdf_values(y_uni, y_w, x_uni, -lambda[i]);

     // Calculate the weighted sum for this specific lambda
     for (int j = 0; j < nx_u; ++j) {
       current_sum += (1.0 - Fc[j]) * ft[j];
     }

     out[i] = current_sum;
   }

   return out;
 }
