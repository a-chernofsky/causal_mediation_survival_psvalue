#' Calculate the restricted mean survival time for interval censored data
#'
#' @param left numeric vector of left endpoints of censoring interval
#' @param right numeric vector of right endpoints of censoring interval
#' @param tau time horizon
#'
#' @return numeric
#' @export
#'
#' @examples
#' library(interval)
#' data(bcos)
#' rmst_ic(bcos$left, bcos$right, tau = 30)
rmst_ic <- function(left, right, tau, subset = NULL){
  if(length(right) != length(left))stop("the left and right vectors must be of the same length")
  
  #if subset option is provided, subset the data
  if(!is.null(subset)){
    left <- left[subset]
    right <- right[subset]
  }
  
  #fit survival curve
  fit <- interval::icfit(left, right)
  
  #vector of times
  time <- c(0,as.vector(fit$intmap))
  
  #times le tau; if tau is between times add tau as the last time
  if(tau %in% time) timetau <- time[time <= tau]
  else  timetau <- c(time[time <= tau], tau)
  
  #calculate survival function at each time point
  St <- getsurv(timetau, fit)[[1]]$S
  
  #store t + 1 and S(t+1)
  timetau_tp1 <- c(timetau[2:length(timetau)], NA)
  St_tp1 <- c(St[2:length(St)], NA)
  
  #estimate of RMST
  sum(0.5*(St + St_tp1)*(timetau_tp1 - timetau), na.rm = T)
}


