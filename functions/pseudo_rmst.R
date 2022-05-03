#' Calculate pseudo-values for restricted mean survival time
#'
#' @param time
#' @param event
#' @param left
#' @param right
#' @param tau
#' @param censor
#'
#' @return
#' @export
#'
#' @examples
pseudo_rmst <- function(time, event, left = NULL, right = NULL, tau, censor = "right"){
  #check that a valid type of censoring is inputted
  if(!(censor %in% c("right", "interval")))stop("A valid censoring type must be provided: try 'right' or 'interval'")
  
  #calculate pseudo-values for right censored data
  if(censor == "right"){
    n <- length(time)
    theta <- rmst_rc(time, event, tau)
    theta_jk <- numeric(length(time))
    for(i in 1:length(theta_jk)){
      theta_jk[i] <- rmst_rc(time, event, tau, subset = -i)
    }
  }
  
  #calculate pseduo-values for interval censored data
  else if(censor == "interval"){
    n <- length(left)
    theta <- rmst_ic(left, right, tau)
    theta_jk <- numeric(length(left))
    for(i in 1:length(theta_jk)){
      theta_jk[i] <- rmst_ic(left, right, tau, subset = -i)
    }
  }
  n*theta - (n-1) * theta_jk
}
