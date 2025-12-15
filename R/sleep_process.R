#' Function to processes sleep runs and extract summary values.
#'
#' @details Currently, requires non-interrupted sleep runs. Any run that
#' touches the beginning or end of
#'
#' @param x Dataframe of values returned by deSolve.
#' @param sleep_var String - name of sleep variable in x data.frame
#' @param time_var String - name of time variable in x data.frame
#'
#' @returns A dataframe with following summary values: Times of sleep onsets, times of sleep midpoints,
#' times of sleep offsets, durations of sleep (in original time scale), start indices of each sleep run,
#' and end indices for each sleep run.
#' @export
#'
#' @examples
sleepSummary <- function(x, sleep_var, time_var){

  ## extract sleep runs ##
  sleep_rle <- rlFunc(x[[sleep_var]])

  ## Drop any run that includes the beginning or end of the data ##
  sleep_rle <- sleep_rle[sleep_rle$value == 1, ] # take only sleep runs
  sleep_rle <- sleep_rle[sleep_rle$start != 1 & sleep_rle$end != nrow(x), ] # drop sleep runs that hit the beginning or end of data

  ## calculate sleep metrics ##
  start_inds <- sleep_rle[, "start"] # index of initial sleep epoch
  end_inds <- sleep_rle[, "end"] # index of final sleep epoch

  sleep_ons <- x[start_inds, time_var] # sleep onsets in original time scale
  sleep_offs <- x[end_inds+1, time_var] # sleep offsets in original scale; note that sleep offset occurs after the index of the final sleep epoch
  sleep_durs <- sleep_offs - sleep_ons # sleep durations in original time units
  sleep_mids <- sleep_ons + sleep_durs/2 # sleep midpoints in original time scale

  ## TODO add ability to consolidate fragmented sleep runs

  ## TODO consider including beginning/end sleep runs but flagging them as potentially incomplete

  return(data.frame(
    sleep_onset = sleep_ons,
    sleep_midpoint = sleep_mids,
    sleep_offset = sleep_offs,
    sleep_duration = sleep_durs,
    start_index = start_inds,
    end_index = end_inds)
    )
}
