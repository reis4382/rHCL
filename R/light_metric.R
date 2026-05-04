#' Determine the duration of observed time for each day-by-offset
#'
#' @param offset_vec A vector of day-by-offset values corresponding to each epoch.
#' @param epoch_length_min Number of minutes in each epoch.
#'
#' @returns A summary data.frame of each calendar day, its start and end times,
#' and the duration in each day.
#' @noRd
#'
dayDurations <- function(offset_vec, epoch_length_min){

  ### Determine hours present in each calendar day ###
  day_summary <- do.call("rbind", lapply(unique(offset_vec), function(x){

    sub_vec <- offset_vec[offset_vec == x] # subset vector
    sub_duration <- length(sub_vec) * epoch_length_min / 60 # hours of observed data

    return(data.frame(
      offset_day = x,
      day_duration = sub_duration
    ))
  }))

  return(day_summary)
}


#' Novel light metric quantifying the biological effect of light
#'
#' This is an implementation of the novel light metric introduced by
#' Skeldon et al. 2023 (see reference). It estimates the impact of 24 hours of
#' light exposure on the circadian pacemaker of an "average" person. Specifically,
#' for each complete calendar day of light data, it estimates the minutes by which
#' the circadian clock of an average person would be expected to speed up or slow
#' down.
#'
#' @param df A data.frame with a column of datetime stamps (in POSIXct format) and
#' a column of light values.
#' @param time_var The name of the datetime column in df.
#' @param light_var The name of the light values column in df.
#' @param epoch_length_min The length of each epoch in minutes.
#' @param y0 A vector with named elements (n, x, and y) corresponding the the
#' initial state of the ordinary differential equation system. This should be
#' left at the default values. While this argument is included to allow for
#' modification of the initial state values, any changes will no longer correspond
#' to the light metric described by Skeldon et al. 2023.
#' @param ode_parms Parameters used for the ordinary differential equations,
#' constructed using the [hclParms()] function. These should also be left at the
#' default value, as any changes will no longer reflect the Skeldon et al. 2023 metric.
#' Note that [hclParms()] has more parameters than needed for the differential equations
#' used by this function, but all parameters (even unused ones) must be passed
#' along to the underlying C functions used by [deSolve::ode()].
#'
#' @returns A data.frame with one row per complete calendar day. The data.frame
#' will include the following elements:
#'  \item{"offset_date"}{The calendar date for the results.}
#'  \item{"offset_day"}{The day of observed data. Calendar days with incomplete
#'  data will be dropped.}
#'  \item{"phase_diff"}{The estimated phase angle between the start and end of the
#'  calendar day with respect to the circadian clock. The units are in minutes.
#'  Positive values mean that there was a net speed up in the circadian clock
#'  over the day (i.e., a phase advance). Negative values mean there was a net
#'  slowdown (i.e., a phase delay).}
#'
#' @export
#'
#' @details Default values for y0 are based on personal correspondence with
#' Professor Skeldon and correspond to the state of the system at midnight for
#' someone with a mid-sleep of 04:30, a sleep duration of 8.22 hours, and a
#' synthetic light profile (see reference for further details).
#'
#' @references Skeldon AC, Rodriguez Garcia T, Cleator SF, Della Monica C,
#' Ravindran KKG, Revell VL, Dijk DJ. Method to determine whether sleep
#' phenotypes are driven by endogenous circadian rhythms or environmental light
#' by combining longitudinal data and personalised mathematical models. PLoS
#' Comput Biol. 2023 Dec 22;19(12):e1011743. doi: 10.1371/journal.pcbi.1011743.
#' PMID: 38134229; PMCID: PMC10817199.
#'
#' @examples
#' # using rhcl_df example data.frame
#' res <- circLight(df = rhcl_df, time_var = "times", light_var = "light",
#'                  epoch_length_min = 1)
#'
circLight <- function(
    df,
    time_var,
    light_var,
    epoch_length_min,
    y0 = c(n = 0.3181, x = -0.8738, y = -0.5534),
    ode_parms = hclParms()

  ){

  # Starting state values are those provided by Prof. Skeldon in
  # personal correspondence

  ## Alternative values for starting states:
  ## I pulled the starting values by running the HCL model and altering mu and
  ## ca_par until I achieved a sleep_midpoint of 4:30 am and a sleep duration of ~8.22 hours
  ## under the default light profile from Skeldon 2017
  # Final values were mu = 18.875 and ca_par = 2.523. All other parameters used
  # the model defaults (e.g., tau_c = 24.2).
  # n = 0.3194685, x = -0.8931994, y = -0.5350999

  ## Need to feed in the same number of parameters as full HCL model to C code,
  # even though homeostatic parameters won't be used. So easier to re-use the
  # same function.

  ### Prepare data.frame and check input ###
  df <- dfPrep(df, time_var = time_var, light_var = light_var)

  ### Break df into complete calendar days ###
  df$daybyoffset <- dayByOffsetVector(df$dtime, hour_offset = 0)
  df$offset_date <- offsetDates(df$dtime, hour_offset = 0)

  ### TODO - consider providing a check for DST transitions within data ###

  ### Identify observed duration for each calendar day ###
  duration_df <- dayDurations(df$daybyoffset, epoch_length_min = epoch_length_min)

  ### Keep only complete days ###
  complete_days <- unique(duration_df$offset_day[duration_df$day_duration == 24]) # all 24 hours must be present

  ### Run ODEs ###
  light_res <- do.call("rbind", lapply(complete_days, function(x){

    ## subset df ##
    sub_df <- df[df$daybyoffset == x, ]

    offset_date <- sub_df$offset_date[1] # take date of offset for day

    # reset time to start at 0 #
    sub_df$ctime <- (sub_df$ctime %% 24)

    ## Set up deSolve arguments - compiled code ##
    desolve_list <- list(
      y = y0,
      times = c(sub_df$ctime, 24), # append 24 to end so that start at next midnight is derived
      func = "derivsc_forger",
      parms = unlist(ode_parms),
      dllname = "rHCL",
      initforc = "forcc_p",
      forcings = cbind(sub_df$ctime, sub_df[[light_var]]),
      fcontrol = list(method = "linear", rule = 2, f = 0),
      initfunc = "parmsc_p",
      nout = 0
    )

    ## run ODE ##
    ode_res <- as.data.frame(do.call(deSolve::ode, desolve_list)) # convert to data.frame

    ## Extract results of system rotation ##
    start_angle <- atan2(ode_res$y[1], ode_res$x[1])
    end_angle <- atan2(ode_res$y[nrow(ode_res)], ode_res$x[nrow(ode_res)])
    # convert phase difference to minutes - positive means the clock surpassed 1 rotation
    phase_diff <- angleDiffs(start_angle, end_angle, period = 2*pi, lbound = -pi) * 12 / pi * 60

    return(data.frame(
     offset_date = offset_date,
     offset_day = x,
     phase_diff = phase_diff
    ))
  }))

  if(!is.null(light_res)){
    ## Return results ##
    # convert NaN to NA #
    light_res$phase_diff[is.nan(light_res$phase_diff)] <- NA

    # calculate mean of phase diffs #
    phase_diff_mean <- mean(light_res$phase_diff, na.rm = TRUE)
    if(is.nan(phase_diff_mean)){
      phase_diff_mean <- NA # convert to NA if no valid days available (i.e., result is nan)
    }
  } else{
    light_res <- NULL
    phase_diff_mean <- NA
  }

  return(list(
    mean_diff = phase_diff_mean,
    df = light_res))
}
