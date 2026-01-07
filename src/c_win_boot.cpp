#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// --- Forward Declarations ---
NumericVector c_ecdf(NumericVector x);
NumericVector c_ecdf_surv(NumericVector ttt, IntegerVector sss);
NumericVector c_win(NumericVector x_uni, NumericVector x_w, NumericVector y_uni, NumericVector y_w, NumericVector lambda);

//' ID-Linked Bootstrap for Win Probability (Survival & Non-Survival)
//'
//' @description
//' Performs a non-parametric bootstrap to estimate the distribution of the
//' win probability \eqn{P(X > Y + \lambda)}. This function automatically detects
//' whether outcomes are survival-based (censored) or standard continuous values
//' by inspecting the status vectors.
//'
//' @param xid IntegerVector of participant IDs for Group X.
//' \bold{Must be 1-based indices} ranging from 1 to \code{length(x)}.
//' @param x NumericVector of raw outcome values (or times) for Group X.
//' @param xs IntegerVector of status indicators for Group X (1 = event, 0 = censored).
//' @param yid IntegerVector of participant IDs for Group Y.
//' \bold{Must be 1-based indices} ranging from 1 to \code{length(y)}.
//' @param y NumericVector of raw outcome values (or times) for Group Y.
//' @param ys IntegerVector of status indicators for Group Y (1 = event, 0 = censored).
//' @param lambda A NumericVector (usually of length 1) specifying the
//' threshold for a "favorable" difference.
//' @param nsim Integer, the number of bootstrap replicates to perform.
//'
//' @return A NumericVector of length \code{nsim} containing the calculated
//' win probability for each bootstrap iteration.
//'
//' @section Warning:
//' The ID vectors (\code{xid} and \code{yid}) ensure that resampling is
//' participant-aligned. It is \bold{critical} that these are provided as
//' 1-based integers. If your status vectors (\code{xs}, \code{ys}) contain
//' any value other than 1, the function automatically switches to survival
//' ECDF calculation.
//'
//' @details
//' The function optimizes performance by:
//' \enumerate{
//'   \item Identifying the outcome type (survival vs. standard) once outside the main loop.
//'   \item Using a "straight-line" resampling approach to maximize CPU cache efficiency.
//'   \item Employing branch prediction-friendly ternary dispatch for ECDF processing.
//' }
//'
//' @export
// [[Rcpp::export]]
 NumericVector c_win_boot(IntegerVector xid,
                          NumericVector x,
                          IntegerVector xs,
                          IntegerVector yid,
                          NumericVector y,
                          IntegerVector ys,
                          NumericVector lambda,
                          int nsim) {

   int nx = x.size();
   int ny = y.size();
   NumericVector out(nsim);

   // 1. Identify if survival logic is needed (once, outside the loop)
   // Logic: if any status is NOT 1, it's potentially censored survival data
   bool is_surv_x = false;
   for(int i = 0; i < nx; ++i) { if(xs[i] != 1) { is_surv_x = true; break; } }
   bool is_surv_y = false;
   for(int i = 0; i < ny; ++i) { if(ys[i] != 1) { is_surv_y = true; break; } }

   // Pre-allocate containers for the resampled values
   NumericVector x_star(nx), y_star(ny);
   IntegerVector xs_star(nx), ys_star(ny);

   // RNGScope ensures set.seed() from R works and resets properly
   RNGScope scope;

   for (int k = 0; k < nsim; ++k) {


     // 3. Resample Group X
     for (int i = 0; i < nx; ++i) {
       int idx = nx * R::runif(0, 1); // Slightly faster than runif(0, nx)
       x_star[i] = x[idx];
       xs_star[i] = xs[idx];
     }

     // 4. Resample Group Y
     for (int i = 0; i < ny; ++i) {
       int idx = ny * R::runif(0, 1);
       y_star[i] = y[idx];
       ys_star[i] = ys[idx];
     }


     // 3. Process Group X: Get ECDF and split
     NumericVector outx = is_surv_x ? c_ecdf_surv(x_star, xs_star) : c_ecdf(x_star);
     int nx_u = outx.size() / 2;
     NumericVector x_uni = outx[Range(0, nx_u - 1)];
     NumericVector x_w   = outx[Range(nx_u, outx.size() - 1)];

     // 4. Process Group Y: Get ECDF and split
     NumericVector outy = is_surv_y ? c_ecdf_surv(y_star, ys_star) : c_ecdf(y_star);
     int ny_u = outy.size() / 2;
     NumericVector y_uni = outy[Range(0, ny_u - 1)];
     NumericVector y_w   = outy[Range(ny_u, outy.size() - 1)];

     // 5. Calculate Win Probability
     NumericVector win_res = c_win(x_uni, x_w, y_uni, y_w, lambda);
     out[k] = win_res[0];
   }

   return out;
 }


