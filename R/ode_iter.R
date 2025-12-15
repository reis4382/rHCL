#' Convergence check for sleep durations and midpoints
#'
#' @param sleep_dur1 First sleep duration value
#' @param sleep_dur2 Second sleep duration value
#' @param sleep_mid1 First sleep midpoint value in 24-hour decimal time
#' @param sleep_mid2 Second sleep midpoint value in 24-hour decimal time
#' @param dur_tol Tolerance allowed for difference between sleep durations. Should be same scale as durations (i.e., hours).
#' @param mid_tol Tolerance allowed for phase angle between midpoints Should be same scale as durations (i.e., hours).
#'
#' @returns
#' @internal
#'
convergeCheck <- function(sleep_dur1, sleep_dur2, sleep_mid1, sleep_mid2, dur_tol, mid_tol){

  ## absolute differences of sleep duration and sleep midpoint##
  # Including NA checks #
  if(is.na(sleep_dur1) | is.na(sleep_dur2)){
    dur_check <- NA
    dur_converge <- FALSE
  } else{
    dur_check <- abs(sleep_dur2 - sleep_dur1)
    ## check against tolerances. Note, to avoid floating point errors, adding a small amount to the defined tolerances ##
    dur_converge <- dur_check <= dur_tol+1e-10
  }

  if(is.na(sleep_mid1) | is.na(sleep_mid2)){
    mid_check <- NA
    mid_converge <- FALSE
  } else{
    mid_check <- abs(clockAngle(sleep_mid1, sleep_mid2))
    mid_converge <- mid_check <= mid_tol+1e-10
  }

  return(list(deviations = data.frame(variable = c("duration", "midpoint"),
                                      devs = c(dur_check, mid_check),
                                      converge = c(dur_converge, mid_converge)),
              overall_converge = dur_converge & mid_converge)
  )

}

#' Checks if initial sleep state is consistent with initial circadian values and homeostatic pressure
#'
#' @param x Initial x state value
#' @param y Initial y state value
#' @param S Initial sleep state value
#' @param h Initial homeostatic sleep pressure value
#' @param hzero Mean sleep pressure parameter
#' @param ca_par Circadian wake propensity amplitude parameter
#' @param delta Separation between thresholds parameter
#'
#' @returns Sleep state value, corrected if necessary
#' @noRD
#'
initialStateCheck <- function(x, y, S, h, hzero, ca_par, delta){
  # new S value #
  new_S <- S # copy

  # circadian propensity #
  circ_prop <- circFunction(x=x, y=y)

  # partial h threshold ignoring sleep v. wake state for now
  h_thresh <- hzero + ca_par * circ_prop

  if(S == 0){
    # threshold if awake
    h_thresh <- h_thresh + 0.5 * delta
    # change new_S if current h is above threshold during wake
    if(h > h_thresh){
      new_S <- 1
    }
  } else if(S == 1){
    # threshold if asleep
    h_thresh <- h_thresh - 0.5 * delta
    # change new_S if current h is below threshold during sleep
    if(h < h_thresh){
      new_S <- 0
    }
  }

  return(new_S) # return starting sleep wake state
}

odeIter <- function(desolve_args, max_iter = 20, dur_tol = 1/60, mid_tol = 1/60){

  ### Check that starting value for sleep pressure is below upper threshold if awake
  desolve_args[["y"]][["S"]] <- initialStateCheck(
    x = desolve_args[["y"]][["x"]],
    y = desolve_args[["y"]][["y"]],
    S = desolve_args[["y"]][["S"]],
    h = desolve_args[["y"]][["h"]],
    hzero = desolve_args[["parms"]][["Hzero"]],
    ca_par = desolve_args[["parms"]][["ca_par"]],
    delta = desolve_args[["parms"]][["delta"]]
  )

  ### Check that initial sleep pressure value isn't greater than mu ###
  desolve_args[["y"]][["h"]] <- min(desolve_args[["y"]][["h"]], desolve_args[["parms"]][["mu"]]) # might not be necessary, but good to be consistent

  ## TODO consider other initial parameter checks, such as negative h
  # - probably not necessary because iterations will remove transients, but could speed up iterations ##

  converge <- FALSE # set converge to FALSE
  iter <- 1 # start iteration counter

  ## prepare results storage ##
  iter_res <- data.frame(
    iteration = c(),
    sleep_midpoint = c(),
    sleep_duration = c())

  ## while loop
  while(!converge & iter <= max_iter){
    ## call deSolve::ode using desolve_args as list of needed arguments
    ode_res <- as.data.frame(do.call(deSolve::ode, desolve_args)) # convert to data.frame

    ## extract sleep summaries ##
    sleep_sum <- sleepSummary(x=ode_res, sleep_var = "S", time_var = "time")

    # add to results data frame #
    iter_res <- rbind(iter_res,
                      data.frame(
                        iteration = iter,
                        sleep_midpoint = timeMean(sleep_sum$sleep_midpoint), # average sleep midpoint
                        sleep_duration = mean(sleep_sum$sleep_duration)
                      ))

    ## Update starting values
    tol <- 1/60/60 # tolerance up to 1 second for matching start time
    tmp.inds <- which((ode_res$time %% 24) > (ode_res$time[1] %% 24) - tol &
                        (ode_res$time %% 24) < (ode_res$time[1] %% 24) + tol) # find matches to time of first observation
    tmp.inds <- tmp.inds[length(tmp.inds)] # take final time index

    # extract new starting state values
    new_y_vals <- c(h = ode_res$h[tmp.inds],
                    n = ode_res$n[tmp.inds],
                    x = ode_res$x[tmp.inds],
                    y = ode_res$y[tmp.inds],
                    S = ode_res$S[tmp.inds])

    ## perform a convergence check ##
    # if first iteration, move on to next, unless only one iteration is requested
    if(iter==1 & max_iter == 1){
      break # break out of while loop without converging

    } else if(iter==1 & max_iter > 1){
      desolve_args[["y"]] <- new_y_vals # update starting y values
      iter <- iter + 1 # increment

    } else {
      ## otherwise, perform convergence check ##
      ode_converge <- convergeCheck(sleep_dur1 = iter_res$sleep_duration[iter-1],
                                    sleep_dur2 = iter_res$sleep_duration[iter],
                                    sleep_mid1 = iter_res$sleep_midpoint[iter-1],
                                    sleep_mid2 = iter_res$sleep_midpoint[iter],
                                    dur_tol = dur_tol,
                                    mid_tol = mid_tol)

      # if iterations converged
      if(ode_converge[["overall_converge"]]){
        converge <- TRUE

      } else if(iter == max_iter){
        break # break out of loop if not converged on max iteration
      } else{
        desolve_args[["y"]] <- new_y_vals # update starting y values
        iter <- iter + 1 # increment
      }
    }
  } # end of while loop

  ## TODO - Return sleep summary and convergence values ##

  ## Convergence messages ##
  if(max_iter == 1){
    conv_message <- "Converge is N/A - Only one iteration requested."
  } else if(converge == FALSE){
    conv_message <- "The model did not converge. You can try increasing the max iterations."
  } else if(converge == TRUE){
    conv_message <- paste("Convergence obtained after", iter, "iterations.")
  }

  ### Final function actions ###

  return(list(ode_res = ode_res, sleepSum = iter_res, converge = converge, conv_message = conv_message, converge_df = ode_converge[["deviations"]], iterations = iter))
}
