#' Kaplan-Meier Wrapper for Cumulative Incidence
#'
#' @description A user-friendly R wrapper for the C++ implementation of the
#' Kaplan-Meier (Product-Limit) estimator. It converts the internal
#' concatenated vector output into a structured R list.
#'
#' @param x A numeric vector of observed times (time-to-event or time-to-censoring).
#' @param xs A numeric vector of status indicators (1 = event, 0 = censored).
#'
#' @return A named list containing two components:
#' \itemize{
#'   \item \code{uni}: A numeric vector of unique event times where the
#'   survival curve "steps" occur.
#'   \item \code{w}: A numeric vector of the cumulative incidence (1 - S(t))
#'   at each unique event time.
#' }
#'
#' @details
#' This function calls the underlying C++ code \code{c_ecdf_surv}. Unlike a
#' standard ECDF, this calculation accounts for right-censoring using the
#' Product-Limit method. The \code{w} component represents the probability
#' that an event has occurred by time \eqn{t}.
#'
#' @seealso \code{\link{c_ecdf_surv}}
#'
#' @examples
#' # Example data: times and event status
#' times <- c(10, 20, 20, 35, 40, 50)
#' status <- c(1, 1, 0, 1, 0, 1)
#'
#' km_res <- winDwinR::r_ecdf_surv(times, status)
#'
#' # Plotting the cumulative incidence step function
#' plot(km_res$uni, km_res$w, type = "s",
#'      xlab = "Time", ylab = "Cumulative Incidence",
#'      main = "KM Estimate")
#'
#' x <- c(0.2,1,2,4,2,2,35,6,3,2,3,1,3,1.5,45,23,12)
#' xs <- c(1,1,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0)
#' dat <- data.frame(x=x,xs=xs)
#' library(survival)
#' fit <- survfit(Surv(x,xs)~1)
#' summary(fit)$time
#' 1-summary(fit)$surv
#'
#' fit1 <- r_ecdf_surv(x,xs)
#'
#' data.frame(time1 = summary(fit)$time,
#'            time2 = fit1$uni,
#'            surv1 = 1-summary(fit)$surv,
#'            surv2 = fit1$w)
#'
#' #++++++ TEST SPEED
#' start.time <- Sys.time()
#' fit <- survfit(Surv(x,xs)~1)
#' end.time <- Sys.time()
#' time.taken1 <- as.numeric(end.time - start.time)
#'
#' start.time <- Sys.time()
#' tf <- r_ecdf_surv(x,xs)
#' end.time <- Sys.time()
#' time.taken2 <- as.numeric(end.time - start.time)
#'
#' time.taken1
#' time.taken2
#' time.taken1/time.taken2
#'
#' @export
#'

r_ecdf_surv <- function(x,xs) {
  tf <- winDwinR:::c_ecdf_surv(x,xs)
  n <- length(tf)/2
  x_uni <- tf[1:n]
  x_w <- tf[(n+1):(2*n)]
  list(uni=x_uni,w=x_w)
}
