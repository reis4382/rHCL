
#' Optimize tau_c parameter on sleep midpoint timing
#'
#' This is the function that is passed to optimize() for minimizing the sum
#' of squares between model estimated sleep midpoint timing and the observed value.
#' The difference between timings is calculated as a phase angle.
#'
#' @param tau_c Parameter value for tau_c passed by optimize().
#' @param mid_sleep Observed sleep midpoint timing for calculating residual.
#' @param desolve_args List of arguments passed to odeIter().
#' @param max_iter Max iterations to be passed to odeIter().
#' @param dur_tol Tolerance allowed for sleep duration (hours) to determine convergence.
#' Passed to odeIter().
#' @param mid_tol Tolerance allowed for sleep midpoint (hours) to determine convergence.
#' Passed to odeIter().
#'
#' @returns The squared residual between estimated and observed sleep midpoint timing.
#' @noRd
#'
odeOptim_midpoint <- function(tau_c, mid_sleep, desolve_args, max_iter, dur_tol, mid_tol){

  ## update tau_c in desolve_args ##
  desolve_args[["parms"]][["tau_c"]] <- tau_c

  ## iterate ##
  ode_res <- odeIter(desolve_args=desolve_args, max_iter = max_iter, dur_tol = dur_tol, mid_tol = mid_tol)

  ## check against observed midsleep time ##
  if(ode_res$converge == FALSE){

    # return 12 if the model did not converge for that value - larger than largest possible phase angle
    return(13^2)
  } else{

    # phase angle - take midpoint of final iteration
    diff <- clockAngle(ode_res$sleep_sum$sleep_midpoint[nrow(ode_res$sleep_sum)], mid_sleep)
    diff_squared <- diff^2
    # print(paste("Values for", tau_c, "are: diff = ", sprintf("%.50f", diff), "and diff^2 = ", diff_squared))
    return(diff_squared)
  }

}


#' Optimize mu parameter on sleep duration
#'
#' @param mu Parameter value for mu passed by optimize().
#' @param sleep_dur Observed sleep duration for calculating residual.
#' @param desolve_args List of arguments passed to odeIter().
#' @param max_iter Max iterations to be passed to odeIter().
#' @param dur_tol Tolerance allowed for sleep duration (hours) to determine convergence.
#' Passed to odeIter().
#' @param mid_tol Tolerance allowed for sleep midpoint (hours) to determine convergence.
#' Passed to odeIter().
#'
#' @returns The squared residual between estimated and observed sleep duration.
#' @noRd
#'
odeOptim_duration <- function(mu, sleep_dur, desolve_args, max_iter = 20,
                              dur_tol = 1/60, mid_tol = 1/60){

  ## update mu in desolve_args ##
  desolve_args[["parms"]][["mu"]] <- mu

  ## iterate ##
  ode_res <- odeIter(desolve_args=desolve_args, max_iter = max_iter, dur_tol = dur_tol, mid_tol = mid_tol)

  ## check against observed midsleep time ##
  if(ode_res$converge == FALSE){

    # return 24^2 if the model did not converge for that value - arbitrary but should be greater than any physiological value
    return(24^2)
  } else{

    # phase angle - take midpoint of final iteration
    diff <- ode_res$sleep_sum$sleep_duration[nrow(ode_res$sleep_sum)] - sleep_dur # difference in sleep durations
    diff_squared <- diff^2 # square of sleep duration difference
    # print(paste("Values for", mu, "are: diff = ", sprintf("%.50f", diff), "and diff^2 = ", diff_squared))

    return(diff_squared)
  }
}



#' Check if sum of squared residuals is < .03
#'
#' @param midpoint_res Residual of estimated vs actual sleep midpoint timing
#' @param duration_res Residual of estimated vs actual sleep duration
#'
#' @returns Boolean - True if sum of squared residuals is < .03.
#' @noRd
#'
residualCheck <- function(midpoint_res, duration_res){
  return((midpoint_res^2 + duration_res^2) < .03)
}
