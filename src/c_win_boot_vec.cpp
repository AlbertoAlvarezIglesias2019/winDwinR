#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// --- Forward Declarations ---
NumericVector c_ecdf(NumericVector x);
NumericVector c_ecdf_surv(NumericVector ttt, IntegerVector sss);
NumericVector c_win(NumericVector x_uni, NumericVector x_w, NumericVector y_uni, NumericVector y_w, NumericVector lambda);

//' Multivariate ID-Linked Bootstrap for Win Probability (Survival & Non-Survival)
//'
//' @description
//' Performs a non-parametric bootstrap to estimate the distribution of the
//' win probability \eqn{P(X > Y + \lambda)} across multiple outcome variables
//' simultaneously. The function automatically detects survival (censored)
//' vs. standard outcomes per column and preserves multivariate correlation
//' via participant-aligned resampling.
//'
//' @param xid IntegerVector of participant IDs for Group X.
//' \bold{Must be 1-based indices} ranging from 1 to \code{nrow(x_mat)}.
//' @param x_mat NumericMatrix where each row is a participant and each column
//' is a distinct outcome variable for Group X.
//' @param xs_mat IntegerMatrix of the same dimensions as \code{x_mat}.
//' Values should be 1 for events/continuous data, and 0 for censored observations.
//' @param yid IntegerVector of participant IDs for Group Y.
//' \bold{Must be 1-based indices} ranging from 1 to \code{nrow(y_mat)}.
//' @param y_mat NumericMatrix where each row is a participant and each column
//' is a distinct outcome variable for Group Y.
//' @param ys_mat IntegerMatrix of the same dimensions as \code{y_mat}.
//' Values should be 1 for events/continuous data, and 0 for censored observations.
//' @param lambda A NumericVector of length \code{ncol(x_mat)} specifying the
//' threshold for a "favorable" difference for each outcome variable.
//' @param nsim Integer, the number of bootstrap replicates to perform.
//'
//' @return A NumericMatrix with \code{nsim} rows and \code{ncol(x_mat)} columns.
//' Each element \code{(k, j)} represents the win probability for the \eqn{j}-th
//' outcome in the \eqn{k}-th bootstrap iteration.
//'
//' @details
//' \bold{Preserving Correlation:} To maintain the joint distribution between
//' different outcomes (columns), the function generates a single set of
//' resample indices for all columns within each bootstrap iteration.
//'
//' \bold{Automated Method Selection:} Before the bootstrap loop, the function
//' scans each column of the status matrices (\code{xs_mat}, \code{ys_mat}).
//' If any value in a column is not 1, that variable is treated as a survival
//' outcome using \code{c_ecdf_surv}; otherwise, it uses the standard \code{c_ecdf}.
//'
//' @note
//' It is assumed that \code{x_mat}, \code{xs_mat}, \code{y_mat}, and \code{ys_mat}
//' all have the same number of columns, representing aligned outcome variables.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix c_win_boot_vec(IntegerVector xid,
                              NumericMatrix x_mat,
                              IntegerMatrix xs_mat,
                              IntegerVector yid,
                              NumericMatrix y_mat,
                              IntegerMatrix ys_mat,
                              NumericVector lambda, // Now expected to be length x_mat.ncol()
                              int nsim) {

  int nx = x_mat.nrow();
  int ny = y_mat.nrow();
  int nj = x_mat.ncol();

  // Output: rows = bootstrap iterations, columns = outcome variables
  NumericMatrix out(nsim, nj);

  // Pre-allocate temporary vectors to reuse in the loops
  NumericVector x_star(nx);
  NumericVector y_star(ny);
  IntegerVector xs_star(nx);
  IntegerVector ys_star(ny);

  LogicalVector is_surv_x(nj, false);
  LogicalVector is_surv_y(nj, false);

  for (int j = 0; j < nj; ++j) {
    for(int i = 0; i < nx; ++i) { if (xs_mat(i,j) != 1) { is_surv_x[j] = true; break; } }
    for(int i = 0; i < ny; ++i) { if (ys_mat(i,j) != 1) { is_surv_y[j] = true; break; } }
  }

  IntegerVector x_idx(nx);
  IntegerVector y_idx(ny);

  RNGScope scope;

  for (int k = 0; k < nsim; ++k) {

    // 1. Generate resample indices once per bootstrap iteration
    // This preserves the correlation between columns (outcomes)
    for (int i = 0; i < nx; ++i) {x_idx[i] = nx*R::runif(0, 1); }
    for (int i = 0; i < ny; ++i) {y_idx[i] = ny*R::runif(0,1); }

    // 2. Loop through each outcome column
    for (int j = 0; j < nj; ++j) {

      // Extract resampled values for column j
      for (int i = 0; i < nx; ++i) {
        x_star[i] = x_mat(x_idx[i], j);
        xs_star[i] = xs_mat(x_idx[i], j);
      }
      for (int i = 0; i < ny; ++i) {
        y_star[i] = y_mat(y_idx[i], j);
        ys_star[i] = ys_mat(y_idx[i], j);
      }

      // 3. Process Group X ECDF
      //NumericVector outx = c_ecdf(x_star);
      NumericVector outx = is_surv_x[j] ? c_ecdf_surv(x_star, xs_star) : c_ecdf(x_star);
      int nx_u = outx.size() / 2;

      // 4. Process Group Y ECDF
      //NumericVector outy = c_ecdf(y_star);
      NumericVector outy = is_surv_y[j] ? c_ecdf_surv(y_star, ys_star) : c_ecdf(y_star);
      int ny_u = outy.size() / 2;

      // 5. Calculate win Probability
      // Passing lambda[j] as a temporary NumericVector
      NumericVector win_res = c_win(outx[Range(0, nx_u - 1)],
                                    outx[Range(nx_u, outx.size() - 1)],
                                    outy[Range(0, ny_u - 1)],
                                    outy[Range(ny_u, outy.size() - 1)],
                                    NumericVector::create(lambda[j]));

      out(k, j) = win_res[0];
    }
  }

  return out;
}
