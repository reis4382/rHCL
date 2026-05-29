#' Convergence check for sleep durations and midpoints
#'
#' @param sleep_dur1 First sleep duration value
#' @param sleep_dur2 Second sleep duration value
#' @param sleep_mid1 First sleep midpoint value in 24-hour decimal time
#' @param sleep_mid2 Second sleep midpoint value in 24-hour decimal time
#' @param dur_tol Tolerance allowed for difference between sleep durations. Should be same scale as durations (i.e., hours).
#' @param mid_tol Tolerance allowed for phase angle between midpoints Should be same scale as durations (i.e., hours).
#'
#' @returns A list with information on overall convergence and remaining deviations.
#' @noRd
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
#' @noRd
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


#' Iterate through ODEs until results converge.
#'
#' @param desolve_args List of arguments needed by [deSolve::ode()]
#' @param dtime_vec Vector of original POSIXct format datetime values.
#' @param max_iter Maximum number of iterations to run
#' @param dur_tol Tolerance of differences in average sleep duration between
#' iterations to determine convergence (in hours)
#' @param mid_tol Tolerance of differences in average sleep midpoint times
#' between iterations to determine convergence (in hours)
#' @param epoch_length_min Numeric value of the length of each epoch in minutes.
#' @param min_observed_hours Minimum hours of data observed for the day, based on
#' epoch_length_min, required for a day to be considered valid for the calculation
#' of sleep statistics.
#'
#' @returns A list with multiple components, including the final ODE results,
#' the summary of sleep values per iteration, and convergence checks.
#' @noRd
#'
odeIterOld <- function(desolve_args, dtime_vec, max_iter, dur_tol, mid_tol,
                    epoch_length_min, min_observed_hours){

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
    ode_res$dtime <- dtime_vec # add original POSIXct datetimes to results

    ## extract sleep summaries ##
    sleep_sum <- sleepSummary(df=ode_res, sleep_var = "S", time_var = "dtime",
                              epoch_length_min = epoch_length_min,
                              min_observed_hours = min_observed_hours)

    # add to results data frame #
    iter_res <- rbind(iter_res,
                      data.frame(
                        iteration = iter,
                        sleep_midpoint = sleep_sum$summary$sleep_mid,
                        sleep_duration = sleep_sum$summary$sleep_dur_noon_24hr # average sleep per day (noon-to-noon days)
                        ## below is code from an old version of sleepSummary
                        # sleep_midpoint = timeMean(sleep_sum$sleep_midpoint), # average sleep midpoint
                        # sleep_duration = mean(sleep_sum$sleep_duration)
                      ))

    ## Update starting values
    tol <- 1/60/60 # tolerance up to 1 second for matching start time
    tmp_inds <- which((ode_res$time %% 24) > (ode_res$time[1] %% 24) - tol &
                        (ode_res$time %% 24) < (ode_res$time[1] %% 24) + tol) # find matches to time of first observation
    tmp_inds <- tmp_inds[length(tmp_inds)] # take final time index

    # extract new starting state values
    new_y_vals <- c(h = ode_res$h[tmp_inds],
                    n = ode_res$n[tmp_inds],
                    x = ode_res$x[tmp_inds],
                    y = ode_res$y[tmp_inds],
                    S = ode_res$S[tmp_inds])

    ## perform a convergence check ##
    # if first iteration, move on to next, unless only one iteration is requested
    if(iter==1 & max_iter == 1){
      ode_converge = list("deviations" = NULL) # no deviations to carry forward
      break # break out of while loop without converging

    } else if(iter==1 & max_iter > 1){
      desolve_args[["y"]] <- new_y_vals # update starting y values
      iter <- iter + 1 # increment

    } else {
      ## otherwise, perform convergence check ##
      ## TODO - build in way of tracking if deviations increase, suggesting
      # that convergence will not be achieved.
      ode_converge <- convergeCheck(sleep_dur1 = iter_res$sleep_duration[iter-1],
                                    sleep_dur2 = iter_res$sleep_duration[iter],
                                    sleep_mid1 = iter_res$sleep_midpoint[iter-1],
                                    sleep_mid2 = iter_res$sleep_midpoint[iter],
                                    dur_tol = dur_tol,
                                    mid_tol = mid_tol)

      # debugging
      # print(paste("dur resid =", sprintf("%.2f", ode_converge$deviations$devs[1])))
      # print(paste("mid resid =", sprintf("%.2f", ode_converge$deviations$devs[2])))

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

  ## Return sleep summary and convergence values ##

  ## Convergence messages ##
  if(max_iter == 1){
    conv_message <- "Converge is N/A - Only one iteration requested."
  } else if(converge == FALSE){
    conv_message <- "The model did not converge. You can try increasing the max iterations."
  } else if(converge == TRUE){
    conv_message <- paste("Convergence obtained after", iter, "iterations.")
  }

  ### Final function actions ###
  # re-arrange ode_res columns #
  other_col_names <- names(ode_res)[!names(ode_res) %in% c("time", "dtime")] # non-time columns
  ode_res <- ode_res[,c("dtime", "time", other_col_names)]

  return(list(ode_res = ode_res, sleep_sum = iter_res, converge = converge, conv_message = conv_message, converge_df = ode_converge[["deviations"]], iterations = iter))
}

#' Replicate data for use with odeIter()
#'
#' @param df Data.frame with column of cumulative time (24-hour format) and column of light values
#' @param ctime_var Name of column in df with cumulative time values
#' @param light_var Name of column in df with light values
#' @param max_ode_iter Number of iterations to run through ODEs for convergence
#' @param tol Tolerance for matching time of day in ctime_var to identify full days.
#'
#' @returns A list containing: 1) data.frame with replications of the full days in df, with full df
#' appended to the end; 2) length of original data; and 3) full days that were replicated.
#' @noRd
#'
odeIterPrep <- function(df, ctime_var, light_var, max_ode_iter, tol){

  orig_length <- nrow(df) # store original length to pass on to odeIter

  # Identify last epoch that matches start time #
  if(max_ode_iter==1){
    res_df <- df
    full_days <- 0 # should make it so that times are not corrected in odeIter()
    final_ind <- nrow(df) # take all data
  } else{
    tmp_inds <- which((df[[ctime_var]] %% 24) > (df[[ctime_var]][1] %% 24) - tol &
                        (df[[ctime_var]] %% 24) < (df[[ctime_var]][1] %% 24) + tol) # find matches to time of first observation
    final_ind <- tmp_inds[length(tmp_inds)] # final observed time (even days)
    full_days <- (df[[ctime_var]][final_ind] - df[[ctime_var]][1]) / 24 # number of full days getting pulled in

    # piece together new times
    new_times <- df[[ctime_var]][1:(final_ind-1)] # new time vector
    day_adj <- rep(1:(max_ode_iter-1), each = length(new_times))-1 # vector for adjusting times
    day_adj <- day_adj * (full_days) * 24 # hours by which to adjust each value
    new_times <- new_times + day_adj # should create appropriate time steps
    new_times <- c(new_times, df[[ctime_var]] + full_days*(max_ode_iter-1)*24) # add full original data to final iteration

    # piece together new light #
    new_light <- df[[light_var]][1:(final_ind-1)] # new light vector
    new_light <- c(rep(new_light, max_ode_iter - 1), df[[light_var]]) # add original data as final iteration

    # prepare output
    res_df <- data.frame(
      var1 = new_times,
      var2 = new_light
    )

    # rename cols
    names(res_df) <- c(ctime_var, light_var)

    # subtract 1 from final_ind (to make it the end of the replicated data)
    final_ind <- final_ind-1
  }


  return(list(
    df = res_df,
    orig_length = orig_length,
    full_days = full_days,
    final_ind = final_ind
  ))
}

#' Iterate through ODEs until results converge.
#'
#' @param desolve_args List of arguments needed by [deSolve::ode()]
#' @param dtime_vec Vector of original POSIXct format datetime values.
#' @param light_vec Vector of original light values.
#' @param max_ode_iter Maximum number of iterations to run
#' @param orig_length Length of original data prior to replication via [odeIterPrep()]
#' @param full_days Number of full days replicated in the data via [odeIterPrep()]
#' @param final_ind Index of last row in original data used during replication of
#' full days in [odeIterPrep()]
#' @param dur_tol Tolerance of differences in average sleep duration between
#' iterations to determine convergence (in hours)
#' @param mid_tol Tolerance of differences in average sleep midpoint times
#' between iterations to determine convergence (in hours)
#' @param epoch_length_min Numeric value of the length of each epoch in minutes.
#' @param min_observed_hours Minimum hours of data observed for the day, based on
#' epoch_length_min, required for a day to be considered valid for the calculation
#' of sleep statistics.
#'
#' @returns A list with multiple components, including the final ODE results,
#' the summary of sleep values per iteration, and convergence checks.
#' @noRd
#'
odeIter <- function(desolve_args, dtime_vec, light_vec, max_ode_iter, orig_length,
                    full_days, final_ind, dur_tol, mid_tol, epoch_length_min,
                    min_observed_hours, sleep_test = TRUE){

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

  ### call desolve::ode using desolve_args as list of needed arguments ###
  ode_res_all <- as.data.frame(do.call(deSolve::ode, desolve_args))

  ### calculate sleep summaries ###
  # original, full data #
  ode_res <- ode_res_all[(nrow(ode_res_all)-orig_length+1):(nrow(ode_res_all)), ]
  ode_res$time <- ode_res$time - (max_ode_iter-1) * full_days * 24 # correct times
  ode_res$dtime <- dtime_vec # add original POSIXct datetimes to results
  ode_res$light <- light_vec # add original light values to results
  row.names(ode_res) <- 1:nrow(ode_res) # fix row.names

  # sleep summary for full data #
  # full_sleep <- sleepSummary(df=ode_res, sleep_var = "S", time_var = "dtime",
  #                            epoch_length_min = epoch_length_min,
  #                            min_observed_hours = min_observed_hours)
  #
  # # summarize #
  # full_sleep_sum <- data.frame(
  #   sleep_midpoint = full_sleep$summary$sleep_mid,
  #   sleep_duration = full_sleep$summary$sleep_dur_noon_24hr
  # )

  full_sleep <- sleepProcessQuick(df = ode_res, sleep_var = "S", time_var = "dtime",
                                  epoch_length_min = epoch_length_min)

  full_sleep_sum <- data.frame(
    sleep_midpoint = full_sleep$sleep_midpoint,
    sleep_duration = full_sleep$sleep_duration
  )

  ## Check convergence if more than 1 iteration
  if(max_ode_iter > 1){
    ## compare ultimate and penultimate iterations ##
    # ultimate iteration #
    iter_ult <- ode_res[1:(final_ind), ] # drop data not used in iterations

    # penultimate iteration #
    pen_ind_end <- nrow(ode_res_all) - orig_length # end of penultimate iteration
    pen_ind_start <- pen_ind_end - nrow(iter_ult) + 1 # start of penultimate iteration

    iter_pen <- ode_res_all[pen_ind_start:pen_ind_end, ] # subset
    iter_pen$time <- iter_pen$time - (max_ode_iter-2) * full_days * 24 # correct times
    iter_pen$dtime <- iter_ult$dtime # add actual datetimes

    ## Sleep summaries ##
    # sleep_sum_ult <- sleepSummary(df=iter_ult, sleep_var = "S", time_var = "dtime",
    #                               epoch_length_min = epoch_length_min,
    #                               min_observed_hours = min_observed_hours)
    sleep_sum_ult <- sleepProcessQuick(df = iter_ult, sleep_var = "S", time_var = "dtime",
                                    epoch_length_min = epoch_length_min)

    # sleep_sum_pen <- sleepSummary(df=iter_pen, sleep_var = "S", time_var = "dtime",
    #                               epoch_length_min = epoch_length_min,
    #                               min_observed_hours = min_observed_hours)
    sleep_sum_pen <- sleepProcessQuick(df = iter_pen, sleep_var = "S", time_var = "dtime",
                                       epoch_length_min = epoch_length_min)

    ## check convergence ##
    # ode_converge <- convergeCheck(sleep_dur1 = sleep_sum_pen$summary$sleep_dur_noon_24hr,
    #                               sleep_dur2 = sleep_sum_ult$summary$sleep_dur_noon_24hr,
    #                               sleep_mid1 = sleep_sum_pen$summary$sleep_mid,
    #                               sleep_mid2 = sleep_sum_ult$summary$sleep_mid,
    #                               dur_tol = dur_tol,
    #                               mid_tol = mid_tol)
    #
    # using sleepProcessQuick() instead of sleepSummary
    ode_converge <- convergeCheck(sleep_dur1 = sleep_sum_pen$sleep_duration,
                                  sleep_dur2 = sleep_sum_ult$sleep_duration,
                                  sleep_mid1 = sleep_sum_pen$sleep_midpoint,
                                  sleep_mid2 = sleep_sum_ult$sleep_midpoint,
                                  dur_tol = dur_tol,
                                  mid_tol = mid_tol)

    # if iterations converged
    if(ode_converge[["overall_converge"]]){
      converge <- TRUE

    } else{
      converge <- FALSE
    }

    ## build iterations data.frame (only keeping penultimate and ultimate) ##
    iter_res <- data.frame(
      iteration = c(max_ode_iter - 1, max_ode_iter),
      # sleep_midpoint = c(sleep_sum_pen$summary$sleep_mid, sleep_sum_ult$summary$sleep_mid),
      # sleep_duration = c(sleep_sum_pen$summary$sleep_dur_noon_24hr, sleep_sum_ult$summary$sleep_dur_noon_24hr)
      sleep_midpoint = c(sleep_sum_pen$sleep_midpoint, sleep_sum_ult$sleep_midpoint),
      sleep_duration = c(sleep_sum_pen$sleep_duration, sleep_sum_ult$sleep_duration)
    )
  } else{
    ## Sleep summaries ##
    # use full data if not comparing against previous full-day iterations
    # sleep_sum_ult <- sleepSummary(df=ode_res, sleep_var = "S", time_var = "dtime",
    #                               epoch_length_min = epoch_length_min,
    #                               min_observed_hours = min_observed_hours)

    full_sleep <- sleepProcessQuick(df = ode_res, sleep_var = "S", time_var = "dtime",
                                    epoch_length_min = epoch_length_min)

    full_sleep_sum <- data.frame(
      sleep_midpoint = full_sleep$sleep_midpoint,
      sleep_duration = full_sleep$sleep_duration
    )

    iter_res <- data.frame(
      iteration = c(max_ode_iter),
      # sleep_midpoint = sleep_sum_ult$summary$sleep_mid,
      # sleep_duration = sleep_sum_ult$summary$sleep_dur_noon_24hr
      sleep_midpoint = full_sleep_sum$sleep_midpoint,
      sleep_duration = full_sleep_sum$sleep_duration
    )

    ode_converge = list("deviations" = NULL) # no deviations to carry forward
    converge <- FALSE # no convergence possible

  }

  ## Convergence messages ##
  if(max_ode_iter == 1){
    conv_message <- "Converge is N/A - Only one iteration requested."
  } else if(converge == FALSE){
    conv_message <- "The model did not converge. You can try increasing the max iterations."
  } else if(converge == TRUE){
    conv_message <- paste("Convergence obtained with", max_ode_iter, "iterations.")
  }

  ### Final function actions ###
  # re-arrange ode_res columns #
  other_col_names <- names(ode_res)[!names(ode_res) %in% c("time", "dtime", "light")] # non-time columns
  ode_res <- ode_res[,c("dtime", "time", "light", other_col_names)]

  return(list(ode_res = ode_res, sleep_sum = full_sleep_sum, iter_sleep_sum = iter_res, converge = converge, conv_message = conv_message, converge_df = ode_converge[["deviations"]], iterations = max_ode_iter))
}
