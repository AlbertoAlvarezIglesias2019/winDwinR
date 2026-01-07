#' Prepare Data for Win/Loss Probability Calculations
#'
#' @description
#' A preprocessing utility that parses formulas (including survival objects),
#' splits data into treatment and control groups, and pre-calculates Empirical
#' Cumulative Distribution Functions (ECDF). It is designed to format data
#' specifically for input into the \code{c_win_boot_vec} and \code{c_loss_boot_vec}
#' C++ functions.
#'
#' @param formu A formula of the type \code{Outcome1 + Surv(time, status) ~ Treatment}.
#'   Multiple outcomes can be provided on the left-hand side separated by \code{+}.
#' @param data A data frame containing the variables specified in \code{formu}.
#' @param treated A string specifying the level of the grouping variable to be
#'   treated as the "Treated" (Group X) arm.
#' @param good A vector of strings (e.g., "high" or "low") indicating whether
#'   higher or lower values of the corresponding outcome are considered favorable.
#'   Defaults to "high".
#'
#' @return A named list where each element corresponds to an outcome variable.
#'   Each element contains:
#' \itemize{
#'   \item \code{xid, yid}: Participant IDs for groups X and Y.
#'   \item \code{x, y}: Outcome values for groups X and Y.
#'   \item \code{xs, ys}: Status indicators (1 = event/continuous, 0 = censored).
#'   \item \code{x_uni, x_w}: Pre-calculated unique values and weights for Group X ECDF.
#'   \item \code{y_uni, y_w}: Pre-calculated unique values and weights for Group Y ECDF.
#'   \item \code{good}: The direction of benefit for that outcome.
#'   \item \code{outtype}: The detected type of outcome (e.g., "Survival", "numeric").
#' }
#'
#' @details
#' \bold{Participant IDs:} The function automatically creates \code{p_i_d} (Participant ID)
#' columns to ensure that when multivariate outcomes are used, the row-wise
#' correlation is maintained during bootstrapping.
#'
#' \bold{Survival Parsing:} If \code{Surv(time, status)} is detected in the formula,
#' the function automatically extracts the time as the outcome and the status
#' as the censoring indicator. For non-survival variables, the status is
#' defaulted to 1.
#'
#' @examples
#' #Example usage:
#' results <- prepwin(Surv(time, death) + weight ~ arm, data = my_data, treated = "Drug A")
#'
#' @export
prepwin <- function(formu,data,treated = "Treated",good="high") {
  data <- data.frame(data)
  #Arm <- as.character(f_rhs(formu))
  Arm <- as.character(formu)[3]
  T1 <- unique(data[,Arm])[1]
  T2 <- unique(data[,Arm])[2]
  wher1 <- data[,Arm]==T1
  wher2 <- data[,Arm]==T2
  data[wher1,"p_i_d"]<-1:sum(wher1)
  data[wher2,"p_i_d"]<-1:sum(wher2)

  #data <- data %>% group_by() %>% mutate(p_i_d = 1:n())

  #Response <- as.character(f_lhs(formu))
  Response <- unlist(str_split(as.character(formu)[2],"[+]"))
  Responses <- str_trim(Response)

  ########################################
  ### Check grouping variable is binary
  ########################################
  temp <- class(data[,Arm])
  if ( !(temp =="character" | temp =="factor")  ) return("ERROR: Grouping variable not character")
  if (length(unique(data[,Arm]))>2 ) return("ERROR: Grouping variable with more than 2 levels")

  temp <- lapply(1:length(Responses),function(k){

    RR <- Responses[k]
    ########################################
    ### First determine the type of outcome
    ########################################
    if (str_detect(RR,"Surv")) outtype <- "Survival"
    if (!str_detect(RR,"Surv") ) outtype <- class(data[,RR])

    if (outtype == "Survival") {
      formu <- as.formula(paste(RR,"~nada",sep=""))
      status <- as.character(formu[[2]][[3]])
      Response <- as.character(formu[[2]][[2]])
    } else {Response <- RR}

    X <- data[data[,Arm]==treated,Response]
    Y <- data[!data[,Arm]==treated,Response]
    Xid <- data[data[,Arm]==treated,"p_i_d"]
    Yid <- data[!data[,Arm]==treated,"p_i_d"]

    if (outtype == "Survival") {
      Xs <- data[data[,Arm]==treated,status]
      Ys <- data[!data[,Arm]==treated,status]
    } else {
      Xs <- rep(1,length(X))
      Ys <- rep(1,length(Y))
    }

    wherx <- is.na(X)
    X <- X[!wherx]
    Xs <- Xs[!wherx]
    Xid <- Xid[!wherx]

    whery <- is.na(Y)
    Y <- Y[!whery]
    Ys <- Ys[!whery]
    Yid <- Yid[!whery]

    if (outtype == "Survival") {
      tempx <- r_ecdf_surv(X,Xs)
      tempy <- r_ecdf_surv(Y,Ys)
    } else {
      tempx <- r_ecdf(X)
      tempy <- r_ecdf(Y)
    }



    list(xid = Xid,x=X,xs=Xs,yid = Yid,y=Y,ys=Ys,
         x_uni=tempx$uni,x_w=tempx$w,y_uni=tempy$uni,y_w=tempy$w,
         good = good[[k]],
         outtype = outtype)
  })

  rnames <- sapply(1:length(Responses),function(k){
    RR <- Responses[k]
    if (str_detect(RR,"Surv")) {
      tem <- str_remove(RR,"Surv[(]")
      tem <- unlist(str_split(tem,","))[1]
      RR <- tem
    }
    RR
  })

  names(temp) <- rnames

  temp

}
