#' Empirical Cumulative Distribution Function Wrapper
#'
#' @description A user-friendly wrapper for the internal C++ ECDF function.
#' It computes the unique sorted values of a numeric vector and their
#' corresponding cumulative probabilities, returning them as a named list.
#'
#' @param x A numeric vector of observations.
#'
#' @return A named list containing two components:
#' \itemize{
#'   \item \code{uni}: A numeric vector of the unique sorted values from \code{x}.
#'   \item \code{w}: A numeric vector of the cumulative probabilities (0 to 1)
#'   corresponding to each unique value.
#' }
#'
#' @details
#' This function acts as a bridge to the \code{c_ecdf} C++ implementation.
#' It handles the memory layout of the C++ output (a concatenated vector)
#' and splits it into a more intuitive list format.
#'
#' @seealso \code{\link{c_ecdf}}
#'
#' @examples
#' data_points <- c(5, 1, 3, 1, 5, 5)
#' result <- r_ecdf(data_points)
#'
#' # Access unique values
#' print(result$uni)
#'
#' # Access cumulative weights
#' print(result$w)
#'
#' # compare to ecdf
#' x <- c(1,2,4,2,2,35,6,3,2,3,1,3)
#' r_ecdf(x)
#' unique(x)
#' ecdf(x)(unique(x))
#'
#' # CHECK SPEED
#' large_data <- rnorm(100000)
#' # Time the C++ Wrapper
#' start_time <- Sys.time()
#'
#'    res_cpp <- r_ecdf(large_data)
#'
#' end_time <- Sys.time()
#' t1 <- as.numeric(end_time - start_time)
#'
#' start_time <- Sys.time()
#'
#'    res_r <- ecdf(large_data)(unique(large_data))
#'
#' end_time <- Sys.time()
#' t2 <- as.numeric(end_time - start_time)
#' cat("\n Speed C++ ",t1,"\n Speed R ",t2,"\n Ratio",t2/t1,"\n")
#'
#' @export
#'

r_ecdf <- function(x) {
  tf <- winDwinR:::c_ecdf(x)
  n <- length(tf)/2
  x_uni <- tf[1:n]
  x_w <- tf[(n+1):(2*n)]
  list(uni=x_uni,w=x_w)
  }
