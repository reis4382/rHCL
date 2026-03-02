#' Create a summary of sleep runs
#'
#' @param sub_df Data.frame with columns corresponding to sleep/wake states (binary)
#' and datetimes (POSIXct).
#' @param sleep_var Name of sleep column.
#' @param time_var Name of ctime variable.
#' @param epoch_length_min Length of each epoch in minutes.
#'
#' @returns A data.frame with a summary of each identified sleep run.
#' @noRd
#'
sleepRunSummary <- function(df, sleep_var, time_var, epoch_length_min){

  # ### create ctime variable from POSIXct dtime ###
  # df$ctime <- ctimeCalc(df[["time_var"]])

  ## Calculate sleep values ##
  sleep_rle <- rlFunc(df[[sleep_var]]) # rle for sleep within given observation day
  sleep_rle <- sleep_rle[sleep_rle$value == 1, ] # extract sleeping runs only within the 24 hours
  start_inds <- sleep_rle[, "start"] # index of initial sleep epoch
  end_inds <- sleep_rle[, "end"] # index of final sleep epoch

  sleep_ons <- df[start_inds, time_var] # sleep onsets in POSIXct times
  # sleep_ons_ctime <- df[start_inds, "ctime"] # sleep onsets in cumulative time

  sleep_offs <- df[end_inds, time_var] + epoch_length_min * 60 # sleep offsets in POSIXct times; note that sleep offset is after final epoch
  # sleep_offs_ctime <- df[end_inds, "ctime"] + epoch_length_min / 60 # sleep offsets in original scale; note that sleep offset occurs after the index of the final sleep epoch

  ## duration and sleep midpoints ##
  sleep_durs <- as.numeric(difftime(sleep_offs, sleep_ons, unit = "hours")) # differences in hours
  sleep_mids <- sleep_ons + (sleep_durs/2) * 60 * 60 # add half the duration after converting to seconds

  # sleep_durs <- sleep_offs_ctime - sleep_ons_ctime # sleep durations in original time units
  # sleep_mids <- sleep_ons_ctimes + sleep_durs/2 # sleep midpoints in original time scale

  sleep_df <- data.frame(
    sleep_onset = sleep_ons,
    sleep_midpoint = sleep_mids,
    sleep_offset = sleep_offs,
    sleep_duration = sleep_durs,
    start_index = start_inds,
    end_index = end_inds)

  return(sleep_df)
}


#' Summarize 24-hour sleep metrics
#'
#' @param df Dataframe with a vector of times (POSIXct format) and sleep/wake states (0 = wake, 1 = sleep).
#' @param sleep_var String - name of sleep variable in df dataframe
#' @param time_var String - name of time variable in df dataframe
#' @param epoch_length_min Length of each epoch for time_var in minutes.
#' @param hour_offset Numeric value representing the hour of the day upon which
#' to divide days. For example, 12 will create noon-to-noon days, 18 will create
#' 6pm-to-6pm days, etc. Must be a value [0, 24).
#'
#' @returns A data frame with the sleep duration and average sleep midpoint for each day.
#' @noRd
#'
sleep24Summary <- function(df, sleep_var, time_var, epoch_length_min, hour_offset){

  ### construct days-by-offset ###
  df$day <- dayByOffsetVector(df[[time_var]], hour_offset = hour_offset)

  ### Calculate starting date for each day-by-offset ###
  df$offset_date <- offsetDates(dtime = df[[time_var]], hour_offset = hour_offset)

  ### calculate cumulative time ###
  ## TODO - send dtime to sleepRunSummary(), and return POSIXct sleep timings in
  ## addition to numeric values?
  # df$ctime <- ctimeCalc(df[[time_var]])

  ### extract day runs ###
  day_rle <- rlFunc(df$day)

#   ## calculate duration of each day ##
#   # subtract start times from end times
#   day_rle$duration <- df[[time_var]][day_rle$end] -
#     df[[time_var]][day_rle$start] # subtract all start times
#
#   ## identify full days, equal to 24 hours minus the minimum observed difference in time between epochs
#   # round sleep durations to account for earlier rounding of epoch_length
#   # going to add a tiny value to increase tolerance for rounding errors
#   day_rle <- day_rle[round(day_rle$duration, 10) >= 24 - epoch_length - 1e-8,]

  ## Calculate 24-hour sleep summaries for each day ##
  day_stats <- do.call("rbind", lapply(day_rle$value, function(i){
    ## calculate duration of observations for each day based on epoch_length_min ##
    sleep_vec <- df[[sleep_var]][df$day==i] # sub vector of sleep values for the day
    observed_hours <- length(sleep_vec) / (60 / epoch_length_min)  # convert to hours

    ## Calculate sleep values ##
    sleep_sum <- sleepRunSummary(df = df[df$day==i, ], sleep_var = sleep_var,
                                 time_var = time_var, epoch_length_min = epoch_length_min)

    # sleep_rle <- rlFunc(sleep_vec) # rle for sleep within given observation day
    # sleep_rle <- sleep_rle[sleep_rle$value == 1, ] # extract sleeping runs only within the 24 hours
    # start_inds <- sleep_rle[, "start"] # index of initial sleep epoch
    # end_inds <- sleep_rle[, "end"] # index of final sleep epoch
    #
    # sleep_ons <- df[df$day==i,][start_inds, "ctime"] # sleep onsets in original time scale
    # sleep_offs <- df[df$day==i,][end_inds, "ctime"] + epoch_length_min / 60 # sleep offsets in original scale; note that sleep offset occurs after the index of the final sleep epoch
    # sleep_durs <- sleep_offs - sleep_ons # sleep durations in original time units
    #
    # # EXPERIMENTAL - average weighted sleep midpoint for day runs #
    # sleep_mids <- sleep_ons + sleep_durs/2 # sleep midpoints in original time scale
    #
    # sleep_mid_wt <- timeMean(sleep_mids, weights = sleep_durs) # weighted sleep midpoint for day

    ### convert sleep midpoints to time-of-day values ###
    sleep_mids <- timeToTOD(sleep_sum$sleep_midpoint)
    # sleep_mids <- lubridate::hour(sleep_sum$sleep_midpoint) +
    #   lubridate::minute(sleep_sum$sleep_midpoint) / 60 +
    #   lubridate::second(sleep_sum$sleep_midpoint) / 60 / 60

    return(data.frame(
      offset_date = df$offset_date[df$day==i][1],
      offset_day = i,
      hour_offset = hour_offset,
      observed_hours = observed_hours,
      sleep_duration = sum(sleep_sum$sleep_duration),
      sleep_midpoint = timeMean(sleep_mids, weights = sleep_sum$sleep_duration) # EXPERIMENTAL - Average of sleep midpoints weighted by their durations
      ))

  }))

  ## TODO consider checking global sleep runs against day boundaries, to identify sleep
  # runs crossing boundaries.
  return(day_stats)

}


#' Function to processes sleep runs and extract summary values.
#'
#' @param df Dataframe with a vector of times (POSIXct format) and sleep/wake states (0 = wake, 1 = sleep).
#' @param sleep_var String - name of sleep variable in df dataframe
#' @param time_var String - name of time variable in df dataframe
#' @param epoch_length_min Numeric value of the length of each epoch in minutes.
#' @param min_observed_hours Minimum hours of data observed for the day, based on
#' epoch_length_min, required for a day to be considered valid for the calculation
#' of sleep statistics. Default is 18 hours.
#'
#' @returns A list with several summaries. The first is a dataframe with following
#' summary values for all sleep runs (excluding those that hit the beginning or end
#' of the data): Times of sleep onsets, times of sleep midpoints, times of sleep offsets,
#' durations of sleep (in original time scale), start indices of each sleep run,
#' and end indices for each sleep run. Dataframes are also returned for
#' 24-hour summaries for both noon-to-noon and midnight-to-midnight values.
#'
#' @export
#'
#' @examples
#'
sleepSummary <- function(df, sleep_var, time_var, epoch_length_min, min_observed_hours = 18){

  ### Data Checks and df prep ###
  df <- dfPrep(df = df, time_var = time_var, sleep_var = sleep_var)

  # # check that df is a data.frame #
  # if(!is(df, "data.frame")){
  #   stop("df must be a data.frame")
  # }
  #
  # # check that sleep_var and time_var are in df #
  # if(!sleep_var %in% names(df)){
  #   stop("Value for sleep_var is not in data.frame df")
  # }
  # if(!time_var %in% names(df)){
  #   stop("Value for time_var is not in data.frame df")
  # }
  #
  # # check for NAs in sleep or time vectors
  # if(sum(is.na(df[[sleep_var]])) > 0){
  #   stop("sleep_var in data.frame df cannot have NAs")
  # }
  #
  # if(sum(is.na(df[[time_var]])) > 0){
  #   stop("time_var in df cannot have NAs")
  # }
  #
  # # check that sleep_var is in correct format (either 0 or 1)
  # if(sum(!unique(df[[sleep_var]] %in% c(0,1))) > 0){
  #   stop("sleep_var in data.frame df must be in binary format (0 = wake, 1 = sleep)")
  # }
  #
  # ## TODO - Build in checks for time variable. Consider function that prepares
  # # a function-appropriate data.frame (e.g., parsing time variables)
  # if(!is(df[[time_var]], "numeric") || sum(df[[time_var]] < 0) > 0){
  #   stop("time_var in data.frame df must be a numeric vector representing cumulative time in 24-hour decimal format.")
  # }

  # # extract epoch length - return error if epochs are not evenly spaced #
  # epoch_lengths <- unique(round(diff(df[[time_var]]), 10)) # round to avoid floating point error

  # if(sum(epoch_lengths %in% NA) > 0){
  #   stop("Differences in time variable includes NA (sleep24Summary())")
  # } else if(length(epoch_lengths) > 1){
  #   stop("Epoch lengths are not evenly spaced (sleep24Summary())")
  # }

  ## extract sleep runs ##
  sleep_df <- sleepRunSummary(df = df, sleep_var = sleep_var,
                              time_var = "dtime", epoch_length_min = epoch_length_min)

  # sleep_rle <- rlFunc(df[[sleep_var]])
  #
  # ## Drop any run that includes the beginning or end of the data ##
  # sleep_rle <- sleep_rle[sleep_rle$value == 1, ] # take only sleep runs
  # sleep_rle <- sleep_rle[sleep_rle$start != 1 & sleep_rle$end != nrow(df), ] # drop sleep runs that hit the beginning or end of data
  #
  # ## calculate sleep metrics ##
  # start_inds <- sleep_rle[, "start"] # index of initial sleep epoch
  # end_inds <- sleep_rle[, "end"] # index of final sleep epoch
  #
  # sleep_ons <- df[start_inds, time_var] # sleep onsets in original time scale
  # sleep_offs <- df[end_inds, time_var] + epoch_lengths # sleep offsets in original scale; note that sleep offset occurs after the index of the final sleep epoch
  # sleep_durs <- sleep_offs - sleep_ons # sleep durations in original time units
  # sleep_mids <- sleep_ons + sleep_durs/2 # sleep midpoints in original time scale
  #
  # # bind results #
  # sleep_df <- data.frame(
  #   sleep_onset = sleep_ons,
  #   sleep_midpoint = sleep_mids,
  #   sleep_offset = sleep_offs,
  #   sleep_duration = sleep_durs,
  #   start_index = start_inds,
  #   end_index = end_inds)

  ## TODO consider adding ability to consolidate fragmented sleep runs

  ## TODO consider including beginning/end sleep runs but flagging them as potentially incomplete

  ### Extract 24-hour summaries ###
  n2n_df <- sleep24Summary(df=df, sleep_var=sleep_var, time_var="dtime",
                           epoch_length_min=epoch_length_min, hour_offset = 12) # noon-to-noon
  m2m_df <- sleep24Summary(df=df, sleep_var=sleep_var, time_var="dtime",
                           epoch_length_min=epoch_length_min, hour_offset = 0) # midnight-to-midnight

  ### keep only days that have required amount of observed data ###
  n2n_df <- n2n_df[n2n_df$observed_hours >= min_observed_hours, ]
  m2m_df <- m2m_df[m2m_df$observed_hours >= min_observed_hours, ]

  ### Calculate primary summary statistics ###
  ## Confirm that values are present in each summary df so that means can be calculated #
  if(nrow(sleep_df) > 0){
    sleep_mid <- timeToTOD(sleep_df$sleep_midpoint)
    sleep_mid <- timeMean(sleep_mid, weights = sleep_df$sleep_duration)
  }

  if(nrow(n2n_df) > 0){
    rownames(n2n_df) <- 1:nrow(n2n_df)
    sleep_dur_noon_24hr <- mean(n2n_df$sleep_duration)
  }

  if(nrow(m2m_df) > 0){
    rownames(m2m_df) <- 1:nrow(m2m_df)
    sleep_dur_midnight_24hr <- mean(m2m_df$sleep_duration)
  }

  summary_df <- data.frame(
    sleep_mid = sleep_mid,
    sleep_dur_noon_24hr = sleep_dur_noon_24hr,
    sleep_dur_midnight_24hr = sleep_dur_midnight_24hr
  )

  return(list(
    summary = summary_df,
    sleep_runs = sleep_df,
    noon_to_noon = n2n_df,
    midnight_to_midnight = m2m_df)
  )
}





