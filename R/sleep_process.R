
#' Summarize 24-hour sleep metrics
#'
#' @param df Dataframe with a vector of times (24-hour decimal format) and sleep/wake states (0 = wake, 1 = sleep).
#' @param sleep_var String - name of sleep variable in df dataframe
#' @param time_var String - name of time variable in df dataframe
#' @param epoch_length Length of each epoch for time_var in 24-hour decimal format.
#' @param noon_to_noon Boolean. If TRUE, will calculate noon-to-noon summaries. If false, will calculate
#' midnight-to-midnight summaries.
#'
#' @returns A data frame with the sleep duration and average sleep midpoint for each day.
#' @noRd
#'
sleep24Summary <- function(df, sleep_var, time_var, epoch_length, noon_to_noon){

  # construct days of observation
  df$day <- epochDays(df[[time_var]], noon_to_noon = noon_to_noon)

  ## extract day runs ##
  day_rle <- rlFunc(df$day)

  ## calculate duration of each day ##
  # subtract start times from end times
  day_rle$duration <- df[[time_var]][day_rle$end] -
    df[[time_var]][day_rle$start] # subtract all start times

  ## identify full days, equal to 24 hours minus the minimum observed difference in time between epochs
  day_rle <- day_rle[day_rle$duration >= (24 - epoch_length),]



  ## Calculate 24-hour sleep summaries for each day ##
  day_stats <- do.call("rbind", lapply(day_rle$value, function(i){

    sleep_rle <- rlFunc(df[[sleep_var]][df$day==i]) # rle for sleep within given observation day
    sleep_rle <- sleep_rle[sleep_rle$value == 1, ] # extract sleeping runs only within the 24 hours
    start_inds <- sleep_rle[, "start"] # index of initial sleep epoch
    end_inds <- sleep_rle[, "end"] # index of final sleep epoch

    sleep_ons <- df[df$day==i,][start_inds, time_var] # sleep onsets in original time scale
    sleep_offs <- df[df$day==i,][end_inds, time_var] + epoch_length # sleep offsets in original scale; note that sleep offset occurs after the index of the final sleep epoch
    sleep_durs <- sleep_offs - sleep_ons # sleep durations in original time units

    # EXPERIMENTAL - average weighted sleep midpoint for day runs #
    sleep_mids <- sleep_ons + sleep_durs/2 # sleep midpoints in original time scale

    sleep_mid_wt <- timeMean(sleep_mids, weights = sleep_durs) # weighted sleep midpoint for day

    ## label for type of day ##
    if(noon_to_noon){
      type_val = "noon-to-noon"
    } else{
      type_val = "midnight-to-midnight"
    }

    return(data.frame(day = i,
                      type = type_val,
                      sleep_duration = sum(sleep_durs),
                      sleep_midpoint = sleep_mid_wt))


  }))

  ## TODO consider checking global sleep runs against day boundaries, to identify sleep
  # runs crossing boundaries.
  return(day_stats)

}


#' Function to processes sleep runs and extract summary values.
#'
#' @param df Dataframe with a vector of times (24-hour decimal format) and sleep/wake states (0 = wake, 1 = sleep).
#' @param sleep_var String - name of sleep variable in df dataframe
#' @param time_var String - name of time variable in df dataframe
#'
#' @returns A list with several summaries. The first is a dataframe with following
#' summary values for all sleep runs (excluding those that hit the beginning or end
#' of the data): Times of sleep onsets, times of sleep midpoints, times of sleep offsets,
#' durations of sleep (in original time scale), start indices of each sleep run,
#' and end indices for each sleep run. Dataframes are also returned for
#' 24-hour summaries for both noon-to-noon and midnight-to-midnight values.
#' @export
#'
#' @examples
sleepSummary <- function(df, sleep_var, time_var){

  ### Data Checks ###
  # check that df is a data.frame #
  if(!is(df, "data.frame")){
    stop("df must be a data.frame")
  }

  # check that sleep_var and time_var are in df #
  if(!sleep_var %in% names(df)){
    stop("Value for sleep_var is not in data.frame df")
  }
  if(!time_var %in% names(df)){
    stop("Value for time_var is not in data.frame df")
  }

  # check for NAs in sleep or time vectors
  if(sum(is.na(df[[sleep_var]])) > 0){
    stop("sleep_var in data.frame df cannot have NAs")
  }

  if(sum(is.na(df[[time_var]])) > 0){
    stop("time_var in df cannot have NAs")
  }

  # check that sleep_var is in correct format (either 0 or 1)
  if(sum(!unique(df[[sleep_var]] %in% c(0,1))) > 0){
    stop("sleep_var in data.frame df must be in binary format (0 = wake, 1 = sleep)")
  }

  ## TODO - Build in checks for time variable. Consider function that prepares
  # a function-appropriate data.frame (e.g., parsing time variables)
  if(!is(df[[time_var]], "numeric") || sum(df[[time_var]] < 0) > 0){
    stop("time_var in data.frame df must be a numeric vector representing cumulative time in 24-hour decimal format.")
  }

  # extract epoch length - return error if epochs are not evenly spaced #
  epoch_lengths <- unique(round(diff(df[[time_var]]), 10)) # round to avoid floating point error

  if(sum(epoch_lengths %in% NA) > 0){
    stop("Differences in time variable includes NA (sleep24Summary())")
  } else if(length(epoch_lengths) > 1){
    stop("Epoch lengths are not evenly spaced (sleep24Summary())")
  }

  ## extract sleep runs ##
  sleep_rle <- rlFunc(df[[sleep_var]])

  ## Drop any run that includes the beginning or end of the data ##
  sleep_rle <- sleep_rle[sleep_rle$value == 1, ] # take only sleep runs
  sleep_rle <- sleep_rle[sleep_rle$start != 1 & sleep_rle$end != nrow(df), ] # drop sleep runs that hit the beginning or end of data

  ## calculate sleep metrics ##
  start_inds <- sleep_rle[, "start"] # index of initial sleep epoch
  end_inds <- sleep_rle[, "end"] # index of final sleep epoch

  sleep_ons <- df[start_inds, time_var] # sleep onsets in original time scale
  sleep_offs <- df[end_inds, time_var] + epoch_lengths # sleep offsets in original scale; note that sleep offset occurs after the index of the final sleep epoch
  sleep_durs <- sleep_offs - sleep_ons # sleep durations in original time units
  sleep_mids <- sleep_ons + sleep_durs/2 # sleep midpoints in original time scale

  # bind results #
  sleep_df <- data.frame(
    sleep_onset = sleep_ons,
    sleep_midpoint = sleep_mids,
    sleep_offset = sleep_offs,
    sleep_duration = sleep_durs,
    start_index = start_inds,
    end_index = end_inds)

  ## TODO consider adding ability to consolidate fragmented sleep runs

  ## TODO consider including beginning/end sleep runs but flagging them as potentially incomplete

  ### Extract 24-hour summaries ###
  n2n_df <- sleep24Summary(df=df, sleep_var=sleep_var, time_var=time_var,
                           epoch_length=epoch_lengths, noon_to_noon = TRUE) # noon-to-noon
  m2m_df <- sleep24Summary(df=df, sleep_var=sleep_var, time_var=time_var,
                           epoch_length=epoch_lengths, noon_to_noon = FALSE) # midnight-to-midnight

  ### Calculate primary summary statistics ##
  summary_df <- data.frame(
    sleep_mid = timeMean(sleep_df$sleep_midpoint, weights = sleep_df$sleep_duration),
    sleep_dur_noon_24hr = mean(n2n_df$sleep_duration),
    sleep_dur_midnight_24hr = mean(m2m_df$sleep_duration)
  )

  return(list(
    summary = summary_df,
    sleep_runs = sleep_df,
    noon_to_noon = n2n_df,
    midnight_to_midnight = m2m_df)
  )
}





