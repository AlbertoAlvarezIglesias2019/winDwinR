#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;


// Forward Declaration
IntegerVector c_order(NumericVector x);

//' Kaplan-Meier Product-Limit Estimator (C++ Implementation)
//'
//' @description Computes the Kaplan-Meier estimate of the survival function
//' (or cumulative incidence) for right-censored data. This implementation is
//' optimized for speed and correctly handles tied event times.
//'
//' @param ttt A numeric vector of observed times (time-to-event or time-to-censoring).
//' @param sss A numeric vector of status indicators (1 = event/death, 0 = censored).
//'
//' @return A numeric vector of length 2 * K, where K is the number of unique
//' event times.
//' \itemize{
//'   \item The first half contains the unique event times (sorted).
//'   \item The second half contains the cumulative incidence (1 - S(t)) at those times.
//' }
//'
//' @details
//' The function calculates survival probability using the formula:
//' \deqn{\hat{S}(t) = \prod_{i: t_i \le t} \left(1 - \frac{d_i}{n_i}\right)}
//' where \eqn{d_i} is the number of events at time \eqn{t_i} and \eqn{n_i}
//' is the number of subjects at risk just before \eqn{t_i}.
//'
//' Note: If an observation is censored at the same time an event occurs, the
//' censored observation is assumed to remain at risk for that event (standard
//' Kaplan-Meier convention).
//'
//' @examples
//' \dontrun{
//' time <- c(1, 3, 3, 5, 7, 10)
//' status <- c(1, 1, 0, 1, 0, 1)
//' res <- c_ecdf_surv(time, status)
//'
//' # Split the results
//' k <- length(res) / 2
//' unique_times <- res[1:k]
//' cum_incidence <- res[(k+1):(2*k)]
//' }
//'
//' @export
// [[Rcpp::export]]
 NumericVector c_ecdf_surv(NumericVector ttt, IntegerVector sss) {
   //int n = ttt.size();
   //if (n == 0) return NumericVector(0);

   int n_orig = ttt.size();
   if (n_orig == 0) return NumericVector(0);

   // 1. Filter out observations with NA in time (ttt)
   // We keep only observations where ttt is not NA.
   std::vector<double> clean_t;
   std::vector<int> clean_s;
   clean_t.reserve(n_orig);
   clean_s.reserve(n_orig);

   for (int i = 0; i < n_orig; ++i) {
     if (!NumericVector::is_na(ttt[i]) && !IntegerVector::is_na(sss[i])) {
       clean_t.push_back(ttt[i]);
       clean_s.push_back(sss[i]);
     }
   }

   int n = clean_t.size();
   if (n == 0) return NumericVector(0);

   // Convert back to Rcpp types for use with your existing c_order function
   NumericVector t_vec = wrap(clean_t);
   IntegerVector s_vec = wrap(clean_s);


   // 2. Get ordering indices (Based on the cleaned data)
   IntegerVector ord = c_order(t_vec);

   std::vector<double> unique_times;
   std::vector<double> survival_probs;

   double current_surv = 1.0;
   int i = 0;

   while (i < n) {
     double current_time = t_vec[ord[i] - 1];
     int events = 0;
     int n_at_risk = n - i;

     // Handle tied time points
     while (i < n && t_vec[ord[i] - 1] == current_time) {
       if (s_vec[ord[i] - 1] > 0) {
         events++;
       }
       i++;
     }

     // Kaplan-Meier Product-Limit Calculation
     if (events > 0) {
       current_surv *= (1.0 - static_cast<double>(events) / n_at_risk);
       unique_times.push_back(current_time);
       survival_probs.push_back(current_surv);
     }
   }

   // 3. Prepare output (Times in first half, 1-Surv in second half)
   int k = unique_times.size();
   NumericVector out(2 * k);
   for (int j = 0; j < k; ++j) {
     out[j] = unique_times[j];
     out[j + k] = 1.0 - survival_probs[j];
   }

   return out;
 }
