#' Evaluate Empirical Cumulative Distribution Function
#'
#' @description
#' A high-performance wrapper for the C++ function \code{c_ecdf_values}.
#' Evaluates a step-function (ECDF) at specific points, optionally applying an offset.
#'
#' @param x_uni Numeric vector. The unique x-coordinates (breakpoints) of the ECDF,
#' typically returned by \code{r_ecdf} or \code{r_ecdf_surv}.
#' @param x_w Numeric vector. The cumulative probabilities (y-coordinates)
#' corresponding to \code{x_uni}.
#' @param z Numeric vector. The values at which to evaluate the ECDF.
#' @param ofset Single numeric value (default = 0). If non-zero, the function
#' evaluates the ECDF at \code{z - ofset}.
#'
#' @return A numeric vector of the same length as \code{z} containing the
#' probabilities evaluated at \code{z - ofset}.
#'
#' @details
#' This function evaluates the ECDF defined by \code{x_uni} and \code{x_w} at
#' \code{z - ofset}. It uses binary search in C++ (\code{std::upper_bound})
#' to ensure \eqn{O(\log n)} lookup time, making it highly efficient for
#' large vectors.
#'
#' @note
#' While results are identical to the base R \code{ecdf} function, this
#' implementation is optimized for integration into larger C++ workflows and
#' handles numeric offsets without creating temporary R vectors.
#'
#' @examples
#' # Basic Usage
#' x <- c(4, 3, 2, 5, 34, 52, 64, 87, 23)
#' tf <- r_ecdf(x)
#' r_ecdf_values(tf$uni, tf$w, 3)
#' r_ecdf_values(tf$uni, tf$w, 2.99)
#'
#' # Comparison with base R ecdf
#' z <- c(2, 4, 23, 100)
#' base_r <- ecdf(x)(z)
#' custom_cpp <- r_ecdf_values(tf$uni, tf$w, z)
#' all.equal(base_r, custom_cpp)
#'
#' #++++++ TEST SPEED
#' x <- rnorm(2000)
#' z <- rnorm(100)
#' start.time <- Sys.time()
#' out <- ecdf(x)(z)
#' end.time <- Sys.time()
#' time.taken1 <- as.numeric(end.time - start.time)
#'
#' start.time <- Sys.time()
#' tf <- r_ecdf(x)
#' out <- r_ecdf_values(tf$uni,tf$w,z)
#' end.time <- Sys.time()
#' time.taken2 <- as.numeric(end.time - start.time)
#'
#' time.taken1
#' time.taken2
#' time.taken1/time.taken2
#'
#' @export
r_ecdf_values <- function(x_uni, x_w, z, ofset = 0) {
  winDwinR::c_ecdf_values(x_uni, x_w, z, ofset)
}
