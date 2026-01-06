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
 NumericVector c_ecdf_surv(NumericVector ttt, NumericVector sss) {
   int n = ttt.size();
   if (n == 0) return NumericVector(0);

   // 1. Get ordering indices (Assuming c_order returns 1-based indices)
   IntegerVector ord = c_order(ttt);

   // 2. Storage for unique time points and survival values
   // We use n as max size and shrink it later
   std::vector<double> unique_times;
   std::vector<double> survival_probs;

   double current_surv = 1.0;
   int i = 0;

   while (i < n) {
     double current_time = ttt[ord[i] - 1]; // Note the -1 for C++ indexing
     int events = 0;
     int censored = 0;
     int n_at_risk = n - i;

     // Handle tied time points: count events and censorings at this exact time
     while (i < n && ttt[ord[i] - 1] == current_time) {
       if (sss[ord[i] - 1] > 0) {
         events++;
       } else {
         censored++;
       }
       i++;
     }

     // Product-Limit Calculation
     // Only update survival and record point if there was at least one event
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
     out[j + k] = 1.0 - survival_probs[j]; // Cumulative Incidence
   }

   return out;
 }
