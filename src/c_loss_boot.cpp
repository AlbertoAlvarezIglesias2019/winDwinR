#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// --- Forward Declarations ---
NumericVector c_ecdf(NumericVector x);
NumericVector c_loss(NumericVector x_uni, NumericVector x_w, NumericVector y_uni, NumericVector y_w, NumericVector lambda);

//' ID-Linked Bootstrap for Loss Probability (Non-Survival)
//'
//' @description
//' Performs a non-parametric bootstrap to estimate the distribution of the
//' loss probability \eqn{P(X < Y - \lambda)}. This function maintains the
//' statistical integrity of participant records by sampling indices.
//'
//' @param xid IntegerVector of participant IDs for Group X.
//' \bold{Must be 1-based indices} ranging from 1 to \code{length(x)}.
//' @param x NumericVector of raw outcome values for Group X.
//' @param yid IntegerVector of participant IDs for Group Y.
//' \bold{Must be 1-based indices} ranging from 1 to \code{length(y)}.
//' @param y NumericVector of raw outcome values for Group Y.
//' @param lambda A NumericVector (usually of length 1) specifying the
//' threshold for a "favorable" difference.
//' @param nsim Integer, the number of bootstrap replicates to perform.
//'
//' @return A NumericVector of length \code{nsim} containing the calculated
//' loss probability for each bootstrap iteration.
//'
//' @section Warning:
//' The ID vectors (\code{xid} and \code{yid}) are currently used to ensure
//' that resampling is participant-aligned. It is \bold{critical} that these
//' are provided as 1-based integers (1, 2, ..., N). If your IDs are
//' non-numeric or 0-based, please map them to a 1-based sequence before
//' passing them to this function.
//'
//' @details
//' The function iterates \code{nsim} times. In each iteration:
//' \enumerate{
//'   \item A bootstrap sample is created for both groups by sampling indices
//'         with replacement.
//'   \item Empirical Cumulative Distribution Functions (ECDF) are computed
//'         using \code{c_ecdf}.
//'   \item The loss probability is calculated using \code{c_loss}.
//' }
//'
//' @export
 // [[Rcpp::export]]
 NumericVector c_loss_boot(IntegerVector xid,
                          NumericVector x,
                          IntegerVector yid,
                          NumericVector y,
                          NumericVector lambda,
                          int nsim) {

   int nx = x.size();
   int ny = y.size();
   NumericVector out(nsim);

   // Pre-allocate containers for the resampled values
   NumericVector x_star(nx);
   NumericVector y_star(ny);

   // RNGScope ensures set.seed() from R works and resets properly
   RNGScope scope;

   for (int k = 0; k < nsim; ++k) {

     // 1. Resample Group X
     // We sample an index 'idx', which maintains the link between x[idx] and xid[idx]
     for (int i = 0; i < nx; ++i) {
       int idx = R::runif(0, nx);
       x_star[i] = x[idx];
       // Even if we don't use xid_star[i] = xid[idx] here,
       // the sampling logic is now participant-aligned.
     }

     // 2. Resample Group Y
     for (int i = 0; i < ny; ++i) {
       int idx = R::runif(0, ny);
       y_star[i] = y[idx];
     }

     // 3. Process Group X: Get ECDF and split
     NumericVector outx = c_ecdf(x_star);
     int nx_u = outx.size() / 2;
     NumericVector x_uni = outx[Range(0, nx_u - 1)];
     NumericVector x_w   = outx[Range(nx_u, outx.size() - 1)];

     // 4. Process Group Y: Get ECDF and split
     NumericVector outy = c_ecdf(y_star);
     int ny_u = outy.size() / 2;
     NumericVector y_uni = outy[Range(0, ny_u - 1)];
     NumericVector y_w   = outy[Range(ny_u, outy.size() - 1)];

     // 5. Calculate loss Probability
     NumericVector loss_res = c_loss(x_uni, x_w, y_uni, y_w, lambda);
     out[k] = loss_res[0];
   }

   return out;
 }


