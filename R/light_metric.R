#' Divide data.frame into calendar days and calculate duration of each one
#'
#' @param df A data.frame with a column of times
#' @param time_var The name of the time column in df. Values must be in 24-hour
#' decimal format.
#'
#' @returns A summary data.frame of each calendar day, its start and end times,
#' and the duration in each day.
#' @noRd
#'
completeDays <- function(df, time_var){

  ### assuming time vector is in 24-hour decimal time and ordered ###

  # TODO - consider a resampling step prior to processing data.frame #

  ### calculate calendar days for data.frame ###
  df$days <- epochDays(df[[time_var]], noon_to_noon = FALSE) # generate calendar days

  ### Determine hours present in each calendar day ###
  # not looking for gaps between successive values #
  day_summary <- do.call("rbind", lapply(unique(df$days), function(x){

    ## extract summary variables ##
    sub_df <- df[df$days==x, ] # subset data.frame
    start_time <- sub_df$time[1] # first observed time in day
    end_time <- sub_df$time[nrow(sub_df)] # second observed time in day
    day_duration <- end_time - start_time # duration between start and end time

    return(data.frame(
      day = x,
      start_time = start_time,
      end_time = end_time,
      duration = day_duration
    ))
  }))

  return(day_summary)
}


circLight <- function(
    df,
    time_var,
    light_var,
    y0 = c(n = 0.3194685, x = -0.8931994, y = -0.5350999),
    ode_parms = hclParms()

  ){

  ## I pulled the starting values by running the HCL model and altering mu and
  ## ca_par until I achieved a sleep_midpoint of 4:30 am and a sleep duration of ~8.22 hours
  ## under the default light profile from Skeldon 2017
  # Final values were mu = 18.875 and ca_par = 2.523. All other parameters used
  # the model defaults (e.g., tau_c = 24.2)


  ## Break df into complete calendar days ##
  df$days <- epochDays(df[[time_var]], noon_to_noon = FALSE) # generate calendar days

  browser()


  ## TODO keep only complete days ##

  ## TODO run ODEs ##

  ## TODO extract results of system rotation ##

  ## TODO return results ##



  ## Need to feed in the same number of parameters as full HCL model to C code,
  # even though homeostatic parameters won't be used. So easier to re-use the
  # same function. Also simplifies things.


}
