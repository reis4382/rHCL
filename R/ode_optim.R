
#' Optimize tau_c parameter on sleep midpoint timing
#'
#' This is the function that is passed to [optimize()] for minimizing the sum
#' of squares between model estimated sleep midpoint timing and the observed value.
#' The difference between timings is calculated as a phase angle.
#'
#' @param tau_c Parameter value for tau_c passed by [optimize()].
#' @param sleep_mid Observed sleep midpoint timing for calculating residual.
#' @param desolve_args List of arguments passed to [odeIter()].
#' @param dtime_vec Vector of original POSIXct format datetime values.
#' @param max_ode_iter Max iterations to be passed to [odeIter()].
#' @param orig_length Length of original data prior to replication via [odeIterPrep()].
#' @param full_days Number of full days replicated in the data via [odeIterPrep()].
#' @param final_ind Index of last row in original data used during replication of
#' full days in [odeIterPrep()].
#' @param dur_tol Tolerance allowed for sleep duration (hours) to determine convergence.
#' Passed to [odeIter()].
#' @param mid_tol Tolerance allowed for sleep midpoint (hours) to determine convergence.
#' Passed to [odeIter()].
#' @param epoch_length_min Numeric value of the length of each epoch in minutes.
#' @param min_observed_hours Minimum hours of data observed for the day, based on
#' epoch_length_min, required for a day to be considered valid for the calculation
#' of sleep statistics.
#'
#' @returns The squared residual between estimated and observed sleep midpoint timing.
#' @noRd
#'
odeOptim_midpoint <- function(tau_c, sleep_mid, desolve_args, dtime_vec, max_ode_iter,
                              orig_length, full_days, final_ind,
                              dur_tol, mid_tol, epoch_length_min, min_observed_hours){

  # print(paste("tau_c", sprintf("%.2f", tau_c))) # debugging

  ## update tau_c in desolve_args ##
  desolve_args[["parms"]][["tau_c"]] <- tau_c

  ## iterate ##
  ode_res <- odeIter(desolve_args=desolve_args, dtime_vec = dtime_vec,
                     max_ode_iter = max_ode_iter, orig_length = orig_length,
                     full_days = full_days, final_ind = final_ind,
                     dur_tol = dur_tol, mid_tol = mid_tol,
                     epoch_length_min = epoch_length_min, min_observed_hours = min_observed_hours)

  ## check against observed midsleep time ##
  if(ode_res$converge == FALSE){

    # return 12 if the model did not converge for that value - larger than largest possible phase angle
    return(13^2)
  } else{

    # phase angle - take midpoint of final iteration
    diff <- clockAngle(ode_res$sleep_sum$sleep_midpoint, sleep_mid)
    diff_squared <- diff^2
    # print(paste("Values for", tau_c, "are: diff = ", sprintf("%.5f", diff), "and diff^2 = ", diff_squared))
    return(diff_squared)
  }

}


#' Optimize mu parameter on sleep duration
#'
#' @param mu Parameter value for mu passed by [optimize()].
#' @param sleep_dur Observed sleep duration for calculating residual.
#' @param desolve_args List of arguments passed to [odeIter()].
#' @param dtime_vec Vector of original POSIXct format datetime values.
#' @param max_ode_iter Max iterations to be passed to [odeIter()].
#' @param orig_length Length of original data prior to replication via [odeIterPrep()].
#' @param full_days Number of full days replicated in the data via [odeIterPrep()].
#' @param final_ind Index of last row in original data used during replication of
#' full days in [odeIterPrep()].
#' @param dur_tol Tolerance allowed for sleep duration (hours) to determine convergence.
#' Passed to [odeIter()].
#' @param mid_tol Tolerance allowed for sleep midpoint (hours) to determine convergence.
#' Passed to [odeIter()].
#' @param epoch_length_min Numeric value of the length of each epoch in minutes.
#' @param min_observed_hours Minimum hours of data observed for the day, based on
#' epoch_length_min, required for a day to be considered valid for the calculation
#' of sleep statistics.
#'
#'
#' @returns The squared residual between estimated and observed sleep duration.
#' @noRd
#'
odeOptim_duration <- function(mu, sleep_dur, desolve_args, dtime_vec, max_ode_iter,
                              orig_length, full_days, final_ind,
                              dur_tol, mid_tol, epoch_length_min, min_observed_hours){

  # print(paste("mu", sprintf("%.2f", mu))) # debugging

  ## update mu in desolve_args ##
  desolve_args[["parms"]][["mu"]] <- mu

  ## iterate ##
  ode_res <- odeIter(desolve_args=desolve_args, dtime_vec = dtime_vec,
                     max_ode_iter = max_ode_iter, orig_length = orig_length,
                     full_days = full_days, final_ind = final_ind,
                     dur_tol = dur_tol, mid_tol = mid_tol,
                     epoch_length_min = epoch_length_min, min_observed_hours = min_observed_hours)

  ## check against observed midsleep time ##
  if(ode_res$converge == FALSE){

    # return 24^2 if the model did not converge for that value - arbitrary but should be greater than any physiological value
    return(24^2)
  } else{

    # phase angle - take midpoint of final iteration
    diff <- ode_res$sleep_sum$sleep_duration - sleep_dur # difference in sleep durations
    diff_squared <- diff^2 # square of sleep duration difference
    # print(paste("Values for", mu, "are: diff = ", sprintf("%.5f", diff), "and diff^2 = ", diff_squared))

    return(diff_squared)
  }
}

# #' Joint optimization of tau and mu
# #'
# #' @param pars Vector of estimates for parameters tau_c and mu
# #' @param sleep_mid Observed sleep midpoint timing for calculating residual.
# #' @param sleep_dur Observed sleep duration for calculating residual.
# #' @param desolve_args List of arguments passed to [odeIter()].
# #' @param max_iter Max iterations to be passed to [odeIter()].
# #' @param dur_tol Tolerance allowed for sleep duration (hours) to determine convergence.
# #' Passed to [odeIter()].
# #' @param mid_tol Tolerance allowed for sleep midpoint (hours) to determine convergence.
# #' Passed to [odeIter()].
# #'
# #' @returns The sum of squared residuals for sleep midpoint and sleep duration estimates.
# #' @noRd
#
# odeOptim_tauAndMu <- function(pars, sleep_mid, sleep_dur, desolve_args,
#                               max_iter, dur_tol, mid_tol){
#
#   ## update tau_c (pars[1]) and mu (pars[2]) in desolve_args ##
#   desolve_args[["parms"]][["tau_c"]] <- pars[1]
#   desolve_args[["parms"]][["mu"]] <- pars[2]
#
#   ## iterate ##
#   ode_res <- odeIter(desolve_args=desolve_args, max_iter = max_iter, dur_tol = dur_tol, mid_tol = mid_tol)
#
#   ## check against observed midsleep time ##
#   if(ode_res$converge == FALSE){
#
#     # return large number if the model did not converge for that value - arbitrary but should be greater than any physiological value
#     return(13^2 + 24^2)
#   } else{
#
#     # return sum of squared residuals #
#     sr_mid <- clockAngle(ode_res$sleep_sum$sleep_midpoint[nrow(ode_res$sleep_sum)], sleep_mid)^2 # squared difference in sleep midpoints
#     sr_dur <- (ode_res$sleep_sum$sleep_duration[nrow(ode_res$sleep_sum)] - sleep_dur)^2 # squared difference in sleep durations
#
#     return(sr_mid+sr_dur)
#   }
# }


#' Check if sum of squared residuals is < .03
#'
#' @param midpoint_res Residual of estimated vs actual sleep midpoint timing
#' @param duration_res Residual of estimated vs actual sleep duration
#'
#' @returns Boolean - True if sum of squared residuals is < .03.
#' @noRd
#'
residualCheck <- function(midpoint_res, duration_res, square){
  if(square){
    midpoint_res <- midpoint_res^2
    duration_res <- duration_res^2
  }
  return((midpoint_res + duration_res) < .03)
}


#' Iterate through nearby parameter values until the Ordinary Differential Equation (ODE)
#' models converge.
#'
#' @param param_val Original value of parameter to be tested.
#' @param param_name Name of parameter being tested (i.e., tau_c or mu)
#' @param lower_bound Lower bound of parameter values to test.
#' @param upper_bound Upper bound of parameter values to test.
#' @param max_steps Maximum number of steps to test between lower and upper bound.
#' @param desolve_args List of arguments passed to [odeIter()].
#' @param dtime_vec Vector of original POSIXct format datetime values.
#' @param max_ode_iter Max iterations to be passed to [odeIter()].
#' @param orig_length Length of original data prior to replication via [odeIterPrep()].
#' @param full_days Number of full days replicated in the data via [odeIterPrep()].
#' @param final_ind Index of last row in original data used during replication of
#' full days in [odeIterPrep()].
#' @param dur_tol Tolerance allowed for sleep duration (hours) to determine convergence.
#' Passed to [odeIter()].
#' @param mid_tol Tolerance allowed for sleep midpoint (hours) to determine convergence.
#' Passed to [odeIter()].
#' @param epoch_length_min Numeric value of the length of each epoch in minutes.
#' @param min_observed_hours Minimum hours of data observed for the day, based on
#' epoch_length_min, required for a day to be considered valid for the calculation
#' of sleep statistics.
#'
#' @returns A list with the results of the converged ODE model and chosen paremeter value. If ODE convergence was not found,
#' returns NAs.
#' @noRd
#'
bisectWhileLoop <- function(param_val, param_name, lower_bound, upper_bound, max_steps,
                            desolve_args, dtime_vec, max_ode_iter,
                            orig_length, full_days, final_ind,
                            dur_tol, mid_tol,
                            epoch_length_min, min_observed_hours){

  # ### check that seq_order is valid ###
  # if(length(seq_order) !=1 | !seq_order %in% c("ascend", "descend", "random")){
  #   stop("seq_order argument needs to be one of 'ascend', 'descend', or 'random'.")
  # }

  ### generate a sequence from the lower to upper bound ###
  jump_seq <- seq(from = lower_bound, to = upper_bound, length.out = max_steps)

  ## If param_val is not in the sequence, add
  if(!param_val %in% jump_seq){
    jump_seq <- c(jump_seq, param_val)
  }

  ## order jumps by distance from param_val, so that smaller changes are tested
  ## before larger ones
  jump_seq <- jump_seq[order(abs(param_val - jump_seq))]

  converge_flag <- FALSE # flag for ODE convergence

  ## iterate through proposed parameter values and test ##
  for(i in jump_seq){
    # print(paste(param_name, sprintf("%.2f", i))) # for debugging

    # update desolve parameters
    desolve_args[["parms"]][[param_name]] <- i
    # run ODEs
    ode_res <- odeIter(desolve_args=desolve_args, dtime_vec = dtime_vec,
                       max_ode_iter = max_ode_iter, orig_length = orig_length,
                       full_days = full_days, final_ind = final_ind,
                       dur_tol = dur_tol, mid_tol = mid_tol,
                       epoch_length_min = epoch_length_min, min_observed_hours = min_observed_hours)

    if(ode_res$converge){
      converge_flag <- TRUE # update flag for later use
      param_val <- i # update parameter value for return
      break # escape for loop if convergence is found
    }
  }

  # ### find a nearby point where the HCL model converges for a given parameter ###
  # jump_count <- 0 # track number of jumps made to find ODE convergence
  # converge_flag <- FALSE # flag for ODE convergence
  # orig_val <- param_val # create variable to track original parameter value
  #
  # ## prepare jump sequence - will sample positive and negative values ##
  # jump_seq <- seq(-1*abs(max_jump), abs(max_jump), length.out = max_steps)
  # # adjust 0s by half the distance between 0 and the next smallest number (with the sign random
  # if(sum(jump_seq==0)>0){
  #   next_smallest_halved <- sort(unique(abs(jump_seq)))[2] / 2 # find half the smallest distance
  #   signs <- rep(sample(c(-1, 1), 1), sum(jump_seq==0)) # randomize signs
  #   jump_seq[jump_seq==0] <- 0 + next_smallest_halved * signs # modify zeros
  # }
  #
  # # arrange sequence order as ascending, descending, or random
  # if(seq_order == "ascend"){
  #   jump_seq <- sort(jump_seq)
  # } else if(seq_order == "descend"){
  #   jump_seq <- sort(jump_seq, decreasing = TRUE)
  # }
  # else if(seq_order == "random"){
  #   jump_seq <- jump_seq[sample(1:length(jump_seq))] # randomize jumps
  # }
  #
  # while(!converge_flag & jump_count <= max_steps){
  #
  #   # update desolve parameters
  #   desolve_args[["parms"]][[param_name]] <- param_val
  #   # run ODEs
  #   ode_res <- odeIter(desolve_args=desolve_args, max_iter = max_ode_iter, dur_tol = dur_tol, mid_tol = mid_tol)
  #
  #   ## check for convergence ##
  #   if(!ode_res$converge){
  #     jump_count <- jump_count + 1 # increment jump_count
  #     param_val <- orig_val + jump_seq[jump_count] # update param_val by modifying original value
  #
  #   } else{
  #     converge_flag <- TRUE
  #   }
  # }

  ### prepare results if non-convergence ###
  if(converge_flag){
    return(list(ode_res = ode_res, param_val = param_val))
  } else {
    return(list(ode_res = NA, param_val = NA))
  }

}


#' Bisection based approach to find roots associated with mu and tau_c
#'
#' This bisection approach has been written to account for non-convergence
#' in the ordinary differential equation iterations during [odeIter()]. Small
#' jumps will be made in the parameter values to attempt to achieve convergence.
#'
#'
#' @param param_lower Lower boundary for estimating parameter.
#' @param param_upper Upper boundary for estimating parameter.
#' @param observed_param Observed value of the parameter for calculating residuals.
#' @param root_stop Value for the squared residual that is considered sufficient
#' for stopping the search. Any parameter value that produces a squared residual
#' less than root_stop will be considered the root.
#' @param max_iter Maximum number of iterations to search for the root.
#' @param abs_tol Absolute tolerance for stopping the root search. The search will
#' stop if half the difference between the new lower and upper bounds is less
#' than abs_tol. No warning is given if search is stopped due to abs_tol in the
#' absence of root_stop being achieved. Absolute tolerance is chosen over relative
#' tolerance because the "midpoint" step c may not actually be the middle of
#' a and b, due to jumps made for [odeIter()] convergence.
#' @param method A value of either "mu" or "tau_c", representing the parameter
#' being estimated.
#' @param num_ode_jumps Maximum number of jumps that will be made in order to
#' address non-convergence of [odeIter()].
#' @param desolve_args List of arguments passed to [odeIter()].
#' @param dtime_vec Vector of original POSIXct format datetime values.
#' @param max_ode_iter Max iterations to be passed to [odeIter()].
#' @param orig_length Length of original data prior to replication via [odeIterPrep()].
#' @param full_days Number of full days replicated in the data via [odeIterPrep()].
#' @param final_ind Index of last row in original data used during replication of
#' full days in [odeIterPrep()].
#' @param dur_tol Tolerance allowed for sleep duration (hours) to determine convergence.
#' Passed to [odeIter()].
#' @param mid_tol Tolerance allowed for sleep midpoint (hours) to determine convergence.
#' Passed to [odeIter()].
#' @param epoch_length_min Numeric value of the length of each epoch in minutes.
#' @param min_observed_hours Minimum hours of data observed for the day, based on
#' epoch_length_min, required for a day to be considered valid for the calculation
#' of sleep statistics.
#' @param update_y0 Experimental argument. Default is TRUE. Will update the starting
#' start after each ODE iteration set to hopefully speed up convergence.
#'
#' @returns A list with two elements designed to copy the output of the [optimize()]
#' function: 1) "minimum" that indicates the parameter
#' value that best matched the root; and 2) "objective" that represents the
#' squared residual of the best fitting parameter.
#' @noRd
#'
odeBisect <- function(param_lower, param_upper, observed_param, root_stop,
                      max_iter, abs_tol, method, num_ode_jumps,
                      desolve_args, dtime_vec, max_ode_iter,
                      orig_length, full_days, final_ind, dur_tol, mid_tol,
                      epoch_length_min, min_observed_hours, update_y0 = TRUE){

  ### check that only mu or tau_c are being checked ###
  if(!method %in% c("mu", "tau_c")){
    stop("method argument for bisection can only be mu or tau_c")
  }

  ### set up function to calculate squared residual based on method ###
  if(method == "mu"){
    # if mu, subtract observed parameter from final sleep_duration summary
    resid_calc <- function(sleep_sum, observed_param){
      resid <- (sleep_sum$sleep_duration - observed_param)
    }
  } else if(method == "tau_c"){
    # if tau_c, use the clockAngle function
    resid_calc <- function(sleep_sum, observed_param){
      resid <- clockAngle(observed_param, sleep_sum$sleep_midpoint)
    }
  }

  ### check that param_higher is greater than param_lower ###
  # appropriate optControl argument for method ##
  if(method == "mu"){
    opt_con <- "duration"
  } else if(method == "tau_c"){
    opt_con <- "midpoint"
  }

  if(param_upper <= param_lower){
    stop(paste("Upper bound (param_upper) for", method, "must be greater than",
               "the lower bound (param_lower) for bisection optimization.",
               "Correct arguments in", paste0(opt_con, "OptControl().")))
  }

  ### initialize f_a and f_b values ###
  f_a <- NA
  f_b <- NA

  ### Check that a value at (or around) the lower bound converges ###
  lower_res <- bisectWhileLoop(param_val = param_lower, param_name = method,
                               lower_bound = param_lower, upper_bound = param_upper,
                               max_steps = num_ode_jumps,
                               desolve_args = desolve_args, dtime_vec = dtime_vec,
                               max_ode_iter = max_ode_iter,
                               orig_length = orig_length, full_days = full_days,
                               final_ind = final_ind,
                               dur_tol = dur_tol, mid_tol = mid_tol,
                               epoch_length_min = epoch_length_min, min_observed_hours = min_observed_hours)

  ## check if lower boundary could be identified ##
  if(is.na(lower_res$param_val)){
    bisect_message <- paste("No value for", method, "at or near the lower boundary could be found that lead to ODE convergence.",
    "Try adjusting the lower boundary (param_lower) in", paste0(opt_con, "OptControl()."),
    "Alternatively, increase the bisect_max_jumps argument.",
    "It is also possible that the ODEs are simply unable to converge at given tau_c value, in",
    "which case look at the arguments dur_tol, and mid_tol, and max_ode_iter.")

    return(list(ode_res = NA, minimum = NA, objective = NA, bisect_message = bisect_message))

  } else{
    val_a <- lower_res$param_val # value for a
    f_a <- resid_calc(lower_res$ode_res$sleep_sum, observed_param) # use function defined at beginning of odeBisect
  }

  if(update_y0){
    ## Update starting values
    tol <- 1/60/60 # tolerance up to 1 second for matching start time
    tmp_inds <- which((lower_res$ode_res$ode_res$time %% 24) > (lower_res$ode_res$ode_res$time[1] %% 24) - tol &
                        (lower_res$ode_res$ode_res$time %% 24) < (lower_res$ode_res$ode_res$time[1] %% 24) + tol) # find matches to time of first observation
    tmp_inds <- tmp_inds[length(tmp_inds)] # take final time index

    # update starting state values
    desolve_args[["y"]] <- c(h = lower_res$ode_res$ode_res$h[tmp_inds],
                    n = lower_res$ode_res$ode_res$n[tmp_inds],
                    x = lower_res$ode_res$ode_res$x[tmp_inds],
                    y = lower_res$ode_res$ode_res$y[tmp_inds],
                    S = lower_res$ode_res$ode_res$S[tmp_inds])
  }

  ### Check that a value at (or around) the upper bound converges ###
  upper_res <- bisectWhileLoop(param_val = param_upper, param_name = method,
                               lower_bound = val_a, upper_bound = param_upper,
                               max_steps = num_ode_jumps,
                               desolve_args = desolve_args, dtime_vec = dtime_vec,
                               max_ode_iter = max_ode_iter,
                               orig_length = orig_length, full_days = full_days,
                               final_ind = final_ind,
                               dur_tol = dur_tol, mid_tol = mid_tol,
                               epoch_length_min = epoch_length_min, min_observed_hours = min_observed_hours)

  ## check if lower boundary could be identified ##
  if(is.na(upper_res$param_val)){
    # I don't believe this error is possible to get, as either the lower bound will
    # cause an error due to lack of identification or the upper bound will be set
    # to the lower bound and therefore pass this check (but fail the next one).
    bisect_message <- paste("No value for", method, "at or near the upper boundary could be found that lead to ODE convergence.",
               "Try adjusting the upper boundary (param_upper) in", paste0(opt_con, "OptControl()."),
               "Alternatively, increase the bisect_max_jumps argument.",
               "It is also possible that the ODEs are simply unable to converge at given tau_c value, in",
               "which case look at the arguments dur_tol, and mid_tol, and max_ode_iter.")

    return(list(ode_res = NA, minimum = NA, objective = NA, bisect_message = bisect_message))

  } else{
    val_b <- upper_res$param_val # value for b
    f_b <- resid_calc(upper_res$ode_res$sleep_sum, observed_param)

  }

  ## check that lower and upper boundaries didn't become identical ##
  if(val_a==val_b){
    bisect_message <- paste("After modification, values for lower and upper bounds became identical.",
               "Try increasing bisect_max_jumps.")

    return(list(ode_res = NA, minimum = NA, objective = NA, bisect_message = bisect_message))
  }

  ### Check if lower and upper boundaries result in opposite signs for residuals ###
  if(sign(f_a) == sign(f_b)){
    if(sign(f_a) == -1){
      bisect_message <- (paste("Boundaries for", method, "parameter estimation result in residuals with the same sign (negative).",
                 "Bisection method will not work. Try increasing the param_upper argument",
                 "in", paste0(opt_con, "OptControl().")))
    } else if(sign(f_a)==1){
      bisect_message <- paste("Boundaries for", method, "parameter estimation result in residuals with the same sign (positive).",
                 "Bisection method will not work. Try decreasing the param_lower argument",
                 "in", paste0(opt_con, "OptControl()."))
    } else if(sign(f_a)==0){
      bisect_message <- paste("Boundaries for", method, "parameter estimation both resulted in residuals of 0.",
                 "This indicates both boundaries match the observed sleep parameter",
                 "That shouldn't happen and means something is not working correctly.")
    }
    return(list(ode_res = NA, minimum = NA, objective = NA, bisect_message = bisect_message))
  }

  ### Check if a root has already been found. If so, return results###
  if(f_a^2 < root_stop){
    return(list(ode_res = lower_res$ode_res, minimum = val_a, objective = f_a^2))
  } else if(f_b^2 < root_stop){
    return(list(ode_res = upper_res$ode_res, minimum = val_b, objective = f_b^2))
  }

  # update starting values if specified #
  if(update_y0){
    ## Update starting values
    tol <- 1/60/60 # tolerance up to 1 second for matching start time
    tmp_inds <- which((upper_res$ode_res$ode_res$time %% 24) > (upper_res$ode_res$ode_res$time[1] %% 24) - tol &
                        (upper_res$ode_res$ode_res$time %% 24) < (upper_res$ode_res$ode_res$time[1] %% 24) + tol) # find matches to time of first observation
    tmp_inds <- tmp_inds[length(tmp_inds)] # take final time index

    # update starting state values
    desolve_args[["y"]] <- c(h = upper_res$ode_res$ode_res$h[tmp_inds],
                             n = upper_res$ode_res$ode_res$n[tmp_inds],
                             x = upper_res$ode_res$ode_res$x[tmp_inds],
                             y = upper_res$ode_res$ode_res$y[tmp_inds],
                             S = upper_res$ode_res$ode_res$S[tmp_inds])
  }

  ### Implement a bisection search for root ###
  root_flag <- FALSE # root flag for while loop
  while_count <- 1 # counter for while loop iterations

  # loop until root is found (within tolerance) or iterations run out
  while(!root_flag & while_count <= max_iter){
    val_c <- (val_a + val_b) / 2 # new value is at midpoint of previous two parameter values

    # if tau_c, set upper and/or lower bounds so that tau_c moves towards 24 in an
    # effort to obtain ode convergence
    if(method == "tau_c"){
      if(val_c < 24){
        new_lower <- val_a + 1e-8 # set lower bound to val_a (slightly above) if val_c is under 24 hours
        new_upper <- min(24, val_b - 1e-8) # set upper bound to either 24 or val_b, whichever is less
      } else if(val_c > 24){
        new_lower <- max(24, val_a + 1e-8) # set lower bound to either 24 or val_a, whichever is greater
        new_upper <- val_b - 1e-8 # set upper bound to val_b if val_c is greater than 24
      } else{
        # if val_c is exactly 24, keep val_a and val_b boundaries
        new_lower <- val_a + 1e-8
        new_upper <- val_b - 1e-8
      }
    } else{
      # if testing mu, keep val_a and val_b as boundaries
      new_lower <- val_a + 1e-8
      new_upper <- val_b - 1e-8
    }

    # test new value for ODE convergence
    c_res <- bisectWhileLoop(param_val = val_c, param_name = method,
                             lower_bound = new_lower, upper_bound = new_upper,
                             max_steps = num_ode_jumps,
                             desolve_args = desolve_args, dtime_vec = dtime_vec,
                             max_ode_iter = max_ode_iter,
                             orig_length = orig_length, full_days = full_days,
                             final_ind = final_ind,
                             dur_tol = dur_tol, mid_tol = mid_tol,
                             epoch_length_min = epoch_length_min, min_observed_hours = min_observed_hours)

    # stop if convergence isn't obtained for c_res
    if(is.na(c_res$param_val)){
      bisect_message <- paste("Convergence could not be obtained for ODEs when testing",
      "c (i.e., values between boundaries) in bisection approach. Look at arguments",
      "that may help ODE convergence, like dur_tol, mid_tol, and max_ode_iter.")
      return(list(ode_res = NA, minimum = NA, objective = NA, bisect_message = bisect_message))
    }

    ## calculate residual ##
    new_val_c <- c_res$param_val # update val_c, based on any ODE convergence jumps
    f_c <- resid_calc(c_res$ode_res$sleep_sum, observed_param)

    ### determine if root has been found, else check for absolute tolerance, else replace either a or b ###
    if(f_c^2 < root_stop){
      root_flag <- TRUE # indicate the root has been found
    } else if(((val_b - val_a)/2) < abs_tol){
      root_flag <- TRUE # indicate that the absolute tolerance has been reached
    } else if(sign(f_a) == sign(f_c)){
      val_a <- val_c # update val_a with val_c if the two have the same sign
      f_a <- f_c # also update residual
    } else{
      val_b <- val_c # update val_b with val_c if the two have the same sign
      f_b <- f_c # also update residual
    }

    # update starting values if specified #
    if(update_y0){
      ## Update starting values
      tol <- 1/60/60 # tolerance up to 1 second for matching start time
      tmp_inds <- which((c_res$ode_res$ode_res$time %% 24) > (c_res$ode_res$ode_res$time[1] %% 24) - tol &
                          (c_res$ode_res$ode_res$time %% 24) < (c_res$ode_res$ode_res$time[1] %% 24) + tol) # find matches to time of first observation
      tmp_inds <- tmp_inds[length(tmp_inds)] # take final time index

      # update starting state values
      desolve_args[["y"]] <- c(h = c_res$ode_res$ode_res$h[tmp_inds],
                               n = c_res$ode_res$ode_res$n[tmp_inds],
                               x = c_res$ode_res$ode_res$x[tmp_inds],
                               y = c_res$ode_res$ode_res$y[tmp_inds],
                               S = c_res$ode_res$ode_res$S[tmp_inds])
    }

  } # end of while loop


  return(list(ode_res = c_res$ode_res, minimum = new_val_c, objective = f_c^2, bisect_message = NULL))
}
