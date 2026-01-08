#' Calculate Win/Loss Summary Statistics
#'
#' @description
#' Iterates through a list of fitted response objects to calculate wins and losses based
#' on a provided threshold. It supports standard win-proportion logic or the
#' Buyse's power-adjusted method.
#'
#' @param fit A named list of objects. Each object should contain:
#'   \itemize{
#'     \item \code{x_uni}, \code{x_w}, \code{y_uni}, \code{y_w}: Data and weights for groups X and Y.
#'     \item \code{good}: A string ("high" or otherwise) indicating if a high value is a favorable outcome.
#'   }
#' @param lambda A numeric threshold (clinical significance margin) used for comparison.
#' @param buyse Logical. If \code{TRUE}, uses \code{c_win_buyse}; otherwise uses \code{c_win}.
#'   Defaults to \code{FALSE}.
#'
#' @return A list of data frames (one for each response variable) containing:
#' \item{Variable}{The name of the response variable.}
#' \item{Lambda}{The threshold used.}
#' \item{Win}{The calculated win proportion.}
#' \item{Loss}{The calculated loss proportion.}
#' \item{WinDiff}{The Net Benefit (Win - Loss).}
#' \item{WinRat}{The Win Ratio (Win / Loss).}
#' \item{NNT}{The Number Needed to Treat, calculated as \code{ceiling(1 / Net Benefit)}.}
#'
#' @export
get_summ <- function(fit,lambda,buyse=FALSE){

  Responses <- names(fit)
  out <- lapply(1:length(Responses),function(k) {

    ff <- fit[[k]]

    if (buyse) {
      www <- c_win_buyse(ff$x_uni,ff$x_w,ff$y_uni,ff$y_w,lambda)
      lll <- c_loss_buyse(ff$x_uni,ff$x_w,ff$y_uni,ff$y_w,lambda)
    } else {
      www <- c_win(ff$x_uni,ff$x_w,ff$y_uni,ff$y_w,lambda)
      lll <- c_loss(ff$x_uni,ff$x_w,ff$y_uni,ff$y_w,lambda)
    }



    if (!ff$good == "high") {
      tt <- www
      www <- lll
      lll <- tt
    }

    data.frame(Variable = Responses[k],
               Lambda = lambda,
               Win = www,
               Loss = lll,
               WinDiff = www-lll,
               WinRat = www/lll,
               NNT = ceiling(1/(www-lll)) )
  })
  names(out) <- Responses
  out

}
