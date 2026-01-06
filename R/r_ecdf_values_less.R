#' Evaluate Empirical Cumulative Distribution Function (Strictly Less Than)
#'
#' @description
#' A wrapper for the C++ function \code{c_ecdf_values_less}. Unlike standard ECDF
#' evaluations, this calculates the probability \eqn{P(X < x)} instead of
#' \eqn{P(X \le x)}, resulting in a left-continuous step function.
#'
#' @param x_uni Numeric vector. The unique x-coordinates (breakpoints) of the ECDF.
#' @param x_w Numeric vector. The cumulative probabilities (y-coordinates)
#' corresponding to \code{x_uni}.
#' @param z Numeric vector. The values at which to evaluate the ECDF.
#' @param ofset Single numeric value (default = 0). If non-zero, the function
#' evaluates the ECDF at \code{z - ofset}.
#'
#' @return A numeric vector of the same length as \code{z} containing the
#' probabilities \eqn{P(X < z - \text{ofset})}.
#'
#' @details
#' This function is essential for workflows requiring the value of the distribution
#' "just before" a jump point. For any point \eqn{z} that exactly matches a
#' value in \code{x_uni}, this function returns the value of the previous step,
#' whereas \code{r_ecdf_values} would return the value of the current step.
#'
#'
#'
#' @note
#' The implementation uses \code{std::lower_bound} in C++ to achieve
#' efficient \eqn{O(\log n)} lookups.
#'
#' @examples
#' x <- c(4, 3, 2, 5, 34, 52, 64, 87, 23)
#' tf <- r_ecdf(x)
#'
#' # Create points exactly at, slightly below, and slightly above jumps
#' xx <- sort(c(tf$uni, tf$uni - 1e-5, tf$uni + 1e-5))
#'
#' # Compare inclusive vs strict evaluation
#' v_inclusive <- r_ecdf_values(tf$uni, tf$w, xx)
#' v_strict    <- r_ecdf_values_less(tf$uni, tf$w, xx)
#'
#' data.frame(time = xx, inclusive = v_inclusive, strict = v_strict)
#'
#' x <- c(4,3,2,5,34,52,64,87,23)
#' tf <- r_ecdf(x)
#' xx <- sort(c(tf$uni,tf$uni-0.00001,tf$uni+0.00001))
#' r_ecdf_values_less(tf$uni,tf$w,xx )
#' r_ecdf_values(tf$uni,tf$w,xx )
#'
#' @export
r_ecdf_values_less <- function(x_uni, x_w, z, ofset = 0) {
  winDwinR:::c_ecdf_values_less(x_uni, x_w, z, ofset)
}
