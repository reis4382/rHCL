#' Base function for converting run length encoding into a dataframe
#'
#' @param vec Vector of values
#'
#' @returns A data.frame with each run of values, their lengths, and their
#' start and ending indices in the original vector.
#' @noRd
#'
#' @examples
#' rlFunc_base(c(0, 0, 1, 0, 0, 1, 1, 1))
#'
rlFunc_base <- function(vec){
  new_rle <- as.data.frame(unclass(rle(vec))) # convert rle to data frame
  new_rle$end <- cumsum(new_rle$lengths) # identifying end indices of runs
  new_rle$start <- new_rle$end - new_rle$lengths + 1 # identify start indices of runs

  # re-order and rename cols
  new_rle <- new_rle[,c("values", "lengths", "start", "end")]
  names(new_rle) <- c("value", "length", "start", "end")

  return(new_rle)
}

#' Run length encoding wrapper that handles NAs
#'
#' @param vec Vector of values
#' @param na_collapse If TRUE (default), will collapse runs of NAs, as [rle()]
#' treats each NA as a separate run.
#'
#' @returns A data.frame with each run of values, their lengths, and their
#' start and ending indices in the original vector. NAs are collapsed into
#' runs.
#' @noRd
#'
#' @examples
#' rlFunc_na(c(0, 0, 1, 0, 0, 1, 1, NA, NA, 1), na_collapse = TRUE)
#'
rlFunc_na <- function(vec, na_collapse = TRUE){
  new_rle <- rlFunc_base(vec) # typical rlFunc

  # Collapse NAs into runs, provided there are NAs #
  if(na_collapse & sum(is.na(vec)) > 0){
    rle_no_na <- new_rle[!is.na(new_rle$value), ] # drop NAs
    rle_na <- rlFunc_base(is.na(vec)) # rle on the NA runs
    rle_na <- rle_na[rle_na$value, ] # only keep NAs
    rle_na$value <- NA # switch these back to NA
    # merge run data frames #
    new_rle <- rbind(rle_no_na, rle_na)
    new_rle <- new_rle[order(new_rle$start), ] # order runs
  }

  # reset row_names #
  row.names(new_rle) <- 1:nrow(new_rle)

  return(new_rle)
}

#' Wrapper for run length encoding that can collapse gaps
#'
#' Wrapper around [rle()] that returns a data.frame. Can also collapse runs
#' of NAs ([rle()] treats each NA as a separate run) and can iterate through
#' runs to collapse on certain gaps.
#'
#' @param vec Vector of values
#' @param na_collapse If TRUE (default), will collapse runs of NAs.
#' @param gap_collapse If TRUE, will collapse gaps of collapse_value that
#' have gaps less than or equal to max_gap argument. Default is FALSE.
#' @param collapse_value Value in original vector (vec argument) that will
#' be collapsed if gap_collapse = TRUE.
#' @param max_gap Maximum length of gaps that will be collapsed, if gap_collapse = TRUE.
#'
#' @returns A data.frame with each run of values, their lengths, and their
#' start and ending indices in the original vector. NAs are collapsed into
#' runs.
#' @export
#'
#' @examples
#' rlFunc(c(0, 0, 1, 0, 0, 1, 1, NA, NA, 1), na_collapse = TRUE,
#' gap_collapse = TRUE, collapse_value = 0, max_gap = 1)
#'
rlFunc <- function(vec, na_collapse = TRUE, gap_collapse = FALSE,
                   collapse_value = NULL, max_gap = NULL){

  ### checks ###
  # check that collapse value is present if gap_collapse = TRUE
  # use && and || to stop logic if first value is false
  if(gap_collapse && (is.null(collapse_value) || !all(collapse_value %in% unique(vec)))){
    stop("If gap_collapse = TRUE, collapse_value cannot be NULL and must be a value found in vec.")
  }

  # check that collapse value is a single value if gap_collapse = TRUE
  if(gap_collapse && length(collapse_value) != 1){
    stop("If gap_collapse = TRUE, collapse_value can only be a single value found in vec.")
  }

  # check that max_gap value is present if gap_collapse = TRUE
  if(gap_collapse && (is.null(max_gap) || !is.numeric(max_gap) || max_gap <= 0)){
    stop("If gap_collapse = TRUE, max_gap cannot be NULL and must be a positive number.")
  }

  ### first run length encode vec ###
  rl_df <- rlFunc_na(vec, na_collapse = na_collapse)

  ### Collapse gaps if requested ###
  if(gap_collapse){

    ### if gap_collapse, identify those that exceed limits ###
    merge_vec <- rep(NA, nrow(rl_df)) # to track which indices should be merged. Will leverage fact that rle() keeps NAs separate

    ## First, identify indices of values upon which to collapse ##
    collapse_inds <- which(!is.na(rl_df$value) & rl_df$value==collapse_value)

    ## proceed if more than one run of desired value ##
    if(length(collapse_inds) > 1){

      # loops through collapse indices to determine merge runs
      for(i in 1:(length(collapse_inds)-1)){
        start_ind <- collapse_inds[i] # start index of potential gap merge
        end_ind <- collapse_inds[i+1] # end index of potential gap merge

        gap_length <- sum(rl_df$length[(start_ind+1):(end_ind-1)]) # length of runs between the two indices
        if(gap_length <= max_gap){
          merge_vec[start_ind:end_ind] <- TRUE # if gap meets max_gap, set those indices to merge
        }

      }

      ## run length encode the gap merge runs
      gap_rl <- rlFunc_na(merge_vec, na_collapse = FALSE) # leverage rle() behavior of keeping NAs separate

      ## reconstruct table ##
      new_values <- rl_df$value[gap_rl$start] # value at start of each new run
      new_starts <- rl_df$start[gap_rl$start] # start ind at start of each new run
      new_ends <- rl_df$end[gap_rl$end] # end in at start of each new run
      new_lengths <- new_ends - new_starts + 1 # new lengths

      rl_df <- data.frame(
        value = new_values,
        length = new_lengths,
        start = new_starts,
        end = new_ends
      )

    }
  }

  return(rl_df)
}


#' Function to find local minima/maxima
#'
#' @param y A vector of values
#' @param minima A logical indicating that minima (default) should be found if TRUE
#' or maxima should be found if FALSE.
#'
#' @returns A vector of indices for the minima/maxima in y. The beginning and end
#' of the vector will not be returned as minima/maxima. In the event of duplicates
#' occurring at the minima/maxima, the first index will be returned for the run.
#'
#' @noRd
#'
minMaxFinder <- function(y, minima = TRUE){
  ## TODO - ensure any erroneous arguments are caught
  ## Argument checks ##
  if(!is.logical(minima)){
    stop("minima argument must be TRUE or FALSE")
  }

  ## Find the local minimums/maximums ##
  diff_y <- diff(c(Inf*(c(-1,1)[minima+1]), y))
  # for minima, we look for derivative switch from negative to positive
  # for maxima, positive to negative
  if(minima){
    diff_y <- diff_y < 0
  } else if(!minima){
    diff_y <- diff_y > 0
  }

  diff_y <- cumsum(rle(diff_y)$lengths) # sum the lengths of each run, grabs the end positions of sign switches
  diff_y <- diff_y[seq.int(1, length(diff_y), 2)] # take every other value?
  if(y[[1]] == y[[2]]){
    diff_y <- diff_y[-1] # drop first element if duplicate starts
  }

  ## Drop mins/maxes at either the beginning or end of the sequence
  diff_y <- diff_y[diff_y != 1 & diff_y != length(y)]

  return(diff_y)
}

#' Calculate homeostatic sleep pressure (H) at time t
#'
#' @param mu \eqn{\mu} parameter representing upper asymptote of sleep pressure
#' @param tswitch Time of last switch from sleep to wake or wake to sleep
#' @param h_tswitch Value of sleep pressure at time of last switch from sleep to wake or wake to sleep
#' @param time Current time (single value or vector)
#' @param chi \eqn{\chi} parameter representing the time constant for the rise and decay of sleep pressure
#' @param s Current state of sleep/wake. S=1 represents sleep while S=0 represents wake.
#'
#' @returns A value or values for sleep pressure (H) at current time.
#' @noRd
#'
sleepHomeostasis <- function(mu, tswitch, h_tswitch, time, chi, s){
  if(s==0){
    h_t <- mu + (h_tswitch - mu)*exp(-(time - tswitch)/chi) # homeostatic value during wake
  } else if (s==1){
    h_t <- h_tswitch * exp(-(time - tswitch)/chi) # homeostatic value during sleep
  }

  return(h_t)
}


#' Shortest angle between two times in 24-hour decimal format.
#'
#' @param vec1 Vector of times in 24-hour decimal format
#' @param vec2 Second vector of times in 24-hour decimal format
#' @param lbound Lower bound of values. Defaults to -12, which indicates
#' possible phase angles are [-12, 12) (exclusive of 12).
#'
#' @returns Vector of phase angles in 24-hour time. Positive values means the second time
#' occurs after the first time.
#' @noRd
#'
clockAngle <- function(vec1, vec2, lbound = -12){
  diffs <- vec2 - vec1 - lbound # adjust the differences, with -12 as the lower bound
  diffs <- diffs %% 24 + lbound # modulo and center
  return(diffs)
}

#' Convert POSIXct vector to numeric time-of-day values
#'
#' @param posix_vec A vector of POSIXct objects
#'
#' @returns A vector of time-of-day values in numeric 24-hour format.
#' @noRd
#'
timeToTOD <- function(posix_vec){

  tod <- lubridate::hour(posix_vec) + lubridate::minute(posix_vec) / 60 +
    lubridate::second(posix_vec) / 60 / 60

  return(tod)

}


#' Mean of a vector of time variables
#'
#' Calculates the mean of a vector of time variables (in 24-hour decimal format).
#' Note this is the conventional way of calculating a circular mean, or the
#' arctangent of the averages of the sines and cosines of the angles. This is
#' not the same as the arithmetic mean for values. For example, while the mean
#' of 4, 4, and 7 is 5, the timeMean will be slightly less than 5.
#'
#' Additionally, NA will be returned if the circular mean is undefined, which happens
#' if the times are equally dispersed around the 24 hours (e.g., midnight, 8, and 16).
#'
#'
#' @param vec Vector of times in 24-hour format.
#' @param weights Optional vector of weights. If NULL (default), will give all values
#' the same weight. The intended usage is to weight sleep midpoint timings by the
#' duration of the sleep episode.
#' @param na_rm Logical. If FALSE (default), will return NA if any NAs are in the vector.
#'
#' @returns The mean of the time vector in 24-hour decimal format
#' @noRd
#'
timeMean <- function(vec, weights = NULL, na_rm = FALSE){

  ## return NA if nothing in vector ##
  if(length(vec)==0){
    return(NA)
  }

  ## check for any NAs if na_rm = FALSE
  non_na_indices <- which(!is.na(vec)) # indices of non-NA values (store for use with weights if needed)

  if(!na_rm){
    if(sum(is.na(vec)) > 0){
      return(NA)
    }
  } else if(na_rm){
    ## if na_rm=TRUE, check if sufficient non-NA values (i.e., at least 2)
    if(sum(!is.na(timevec))<2){
      return(NA)
    }
    # remove remaining NAs
    vec <- vec[non_na_indices] # keep non-na values
  }

  ## checks for weight vector ##
  if(!is.null(weights)){
    # check that weights are numeric
    if(!is(weights, "numeric")){
      stop("Weights provided to timeMean() function must be numeric")
    }
    # check that weights match length of input
    if(length(weights) != length(vec)){
      stop("Vector of weights must be the same length as input (timeMean())")
    }
  }

  ### convert times to radians ###
  rvec <- vec * pi / 12

  # ### mean of radians (without weighting) ##
  # part1 <- 1/length(rvec) * sum(sin(rvec)) # first part of atan2
  # part2 <- 1/length(rvec) * sum(cos(rvec)) # second part of atan2

  ### weighted mean of radians ###
  # If no weights are provided, all values given same weight #
  if(is.null(weights)){
    w_vec <- rep(1/length(rvec), length(rvec))
  } else{
    w_vec <- weights[non_na_indices] / sum(weights[non_na_indices])
  }

  part1 <- sum(w_vec * sin(rvec)) # first part of atan2
  part2 <- sum(w_vec * cos(rvec)) # second part of atan2


  ### check if values result in an undefined circular mean (e.g., 0 and 12, which has a center in the center of the circle) ###
  ## I believe the sums of the sines and cosines should not both equal 0 unless the mean is undefined.
  if(all(isTRUE(all.equal(part1, 0)), isTRUE(all.equal(part2, 0)))){
    warning("Mean of times is undefined")
    return(NA)
  }

  ### atan of mean sines and cosines ###
  mang <- atan2(part1, part2) # mean angle (in radians)

  ### for negative signs, find positive equivalent
  if(sign(mang)==-1){
    mang <- 2*pi + mang # 2pi radians in 24 hours
  }

  ### convert back to 24-hour decimal time ###
  mang <- (mang * 12 / pi) %% 24 # modulus to handle edge cases where mang is negative but essentially 0
  return(mang)
}

#' Function to generate average light profile.
#'
#' Equation #11 in supplementary materials of Skeldon et al. 2017. Note that the
#' formula appears to have an error: I believe the difference of the two tanh
#' components should be together in parentheses prior to the multiplication.
#'
#'
#' @param t Vector of times in 24-hour decimal format. Cumulative time is allowed.
#' @param l1 Light value (lux) for bright time of day. Default value of 700 taken
#' from paper and corresponds with summer brightness levels. A value of 300 could
#' be used to correspond to winter levels.
#' @param l2 Light value (lux) for dim parts of day. Default value of 40 taken from
#' paper.
#' @param c Steepness of light profile transitions. Note that paper indicates
#' the default value (Table S1) is 1/6000. This does not work when the time
#' scale is in hours, as it results in practically no change in light values.
#' It could be in seconds
#' @param s1 Time that profile changes from l2 to l1. Default value taken from
#' paper and corresponds to a roughly 12-daylight duration centered on noon.
#' @param s2 Time that profile changes from l1 to l2.Default value taken from
#' paper and corresponds to a roughly 12-daylight duration centered on noon.
#' @param time_scale Scale of the time variable in hours. Three character values are allowed:
#' "hours" (default), "mins", and "secs", which indicate the time variable is
#' scaled to hours, minutes, or seconds of the day respectively.
#'
#'
#' @returns A vector of light values for each time point in t
#'
#' @references Skeldon AC, Phillips AJ, Dijk DJ. The effects of self-selected light-dark
#' cycles and social constraints on human sleep and circadian timing: a modeling
#' approach. Sci Rep. 2017 Mar 27;7:45158. doi: 10.1038/srep45158.
#' PMID: 28345624; PMCID: PMC5366875.
#'
#' @export
#'
#' @example
#' t_vec <- seq(from=0, to=24, by=0.1)
#' light_vals <- lightCycle(t_vec)
#'
lightCycle <- function(t, l1 = 700, l2 = 40, c = 0.6, s1 = 7.5, s2 = 16.5, time_scale = "hours"){

  ## extract time scale adjustment ##
  if(time_scale == "hours"){
    time_scale <- 1
  } else if(time_scale == "mins"){
    time_scale = 60
  } else if(time_scale == "secs"){
    time_scale = 60*60
  } else{
    stop("time_scale argument must be either one of the following: 'hours', 'mins', 'secs'")
  }

  ## adjust parameters for time_scale ##
  c <- c / time_scale
  s1 <- s1 * time_scale
  s2 <- s2 * time_scale
  time_modulo <- 24 * time_scale

  ## l1 and l2 difference component ##
  lvars <- (l1 - l2) / 2

  ## tanh part 1 ##
  tan1 <- tanh(c* ((t%%time_modulo) - s1))

  ## tanh part 2 ##
  tan2 <- tanh(c * ((t%%time_modulo) - s2))

  return(l2 + lvars * (tan1 - tan2))
}

#' Calculate noon-to-noon or midnight-to-midnight days
#'
#' @param time_vec Vector of times in 24-hour decimal format.
#' @param noon_to_noon Logical. If TRUE, days will be constructed for noon-to-noon (i.e., 24 hours resetting at 12).
#' If FALSE, days will be constructed for midnight-to-midnight.
#'
#' @returns A vector of the days of observation. Begins at day 1.
#' @noRd

epochDays <- function(time_vec, noon_to_noon){

  ### specify days for each epoch ###
  if(noon_to_noon){
    offset <- 12 # offset for noon-to-noon days
  } else{
    offset <- 0 # offset for midnight-to-midnight days
  }

  # subtract initial time to start data at day 1
  days <- 1 + floor((time_vec + offset) / 24) # should add 1 at 12, 36, etc for noon-to-noon and 24, 48, etc. for midnights

  ## reset day count at 1
  days <- 1 + (days - days[1])

  return(days)

}


#' Return a daybyoffset vector based on dtime and desire hour-of-day offset
#'
#' @param dtime A vector of POSIXct objects
#' @param hour_offset Desired hour-of-day offset. E.g., 12 would create noon-to-noon
#' day, and 18 would create 6pm-to-6pm days.
#'
#' @returns A vector of length(dtime) with the corresponding day-by-offset.
#' @export
#'
#' @examples
#' \dontrun{
#'   # TBD
#' }
dayByOffsetVector <- function(dtime, hour_offset){

  if(hour_offset < 0 | hour_offset >=24){
    stop("Hour offset for distinguishing days must be >= 0 or less than 24")
  }

  if(!is(dtime, "POSIXct")){
    stop("dtime argument must be a POSIXct object")
  }

  # strip out dates - otherwise, time zone conversions will treat dates as UTC #
  string_dates <- as.character(dtime)
  string_dates <- gsub(" .*$", "", string_dates) # remove all characters from first space until end

  # convert to new dates (which will assume UTC, but that should be fine)
  new_dates <- as.Date(string_dates, format = "%Y-%m-%d")

  start_date <- as.Date(new_dates[1]) # take first date
  date_diffs <- as.numeric(difftime(new_dates, start_date, units = "days")) # differences in calendar dates since initial recording

  clock_times <- lubridate::hour(dtime) + lubridate::minute(dtime)/60 + lubridate::second(dtime) / 60 / 60 # clock times in 24-hour decimal format

  ctime <- clock_times + date_diffs * 24 # add calendar date offsets

  daybyoffset <- floor((ctime - hour_offset) / 24) + 1 # day-by-offset

  # if first time is prior to offset value (e.g., 11 am for a 12 hour offset),
  # it will return a negative value in the previous line (then adjusted to 0), so correct here
  if(floor((ctime[1] - hour_offset)/24) < 0){
    daybyoffset <- daybyoffset + 1 # adding one should fix the issue
  }

  return(daybyoffset)
}

#' Return a starting date of the 24-hour window accounting for desired offset
#'
#' @param dtime A POSIXct object
#' @param hour_offset Numeric. Hours of the day used for offset (e.g., 12 = noon).
#'
#' @returns A date.
#' @noRd
#'
offsetDates <- function(dtime, hour_offset){

  res_date <- as.Date(dtime) # dates
  res_date[lubridate::hour(dtime) < hour_offset] <- res_date[lubridate::hour(dtime) < hour_offset] - lubridate::days(1)

  return(res_date)

}


#' Separate times into 24-hour days (per clock time) based on desire hour offset
#'
#' This function will separate timestamps into "days-by-offset", based on a
#' clock time. For example, an hour_offset = 0 will separate days into
#' midnight-to-midnight days, which corresponds to calendar date. An hour_offset = 12
#' would separate days into noon-to-noon.
#'
#' Note that daylight savings transition dates will be separated in the same way,
#' meaning that these days will have greater than or fewer than 24 hours.
#'
#' @param dtimes A vector of timestamps in POSIXct format.
#' @param hour_offset A numeric value representing the desired hour-of-day upon
#' which to separate days. The value must be between 0 (inclusive) and 24 (exclusive).
#'
#' @returns A data.frame with two columns corresponding to the length of dtimes.
#' The first column ("day_by_offset") will be the day of data. The second column
#' ("offset_date") will the be calendar date for the beginning of the day-by-offset.
#' @noRd
#'
daySplit <- function(dtimes, hour_offset){

  if(hour_offset < 0 | hour_offset >= 24){
    stop("hour_offset must be a value >= 0 and < 24")
  }

  # ## extract time zone ##
  # time_zone <- attr(dtimes, "tzone")

  ### Calculate day-by-offset, based on calendar date/clock time ###
  ## grab calendar day for each value, as number of days since start, + 1 ##
  # calendar_day <- as.numeric(as.Date(dtimes) - as.Date(dtimes[1])) + 1
  ## TODO FIX THE CALENDAR DATE DIFFERENCES. THIS METHOD (using difftime or just - )
  ## WILL ASSUME THE DATES ARE IN UTC, WHICH MEANS THAT CALENDAR DATE WILL BE WRONG
  ## IN OTHER TIME ZONES

  ## Calculate time-of-day based on observed POSIXct time ##
  time_of_day <- lubridate::hour(dtimes) + lubridate::minute(dtimes) / 60 + lubridate::second(dtimes) / 60 / 60

  ## Calculate the "duration" of each time step. Final dtime has a duration of 0 ##
  step_durations <- c(as.numeric(difftime(dtimes[2:length(dtimes)], dtimes[1:(length(dtimes)-1)], units = "hours")), 0)

  # if using midnight as offset, day is simply the calendar day calculation
  if(hour_offset == 0){
    day_by_offset <- calendar_day
  } else{

    ## identify values >= hour_offset
    offset_flag <- time_of_day >= hour_offset

    ## combine calendar day and offset flag values
    day_by_offset <- calendar_day + offset_flag
  }

  # ### identify day_by_offset start and end indices ###
  # dbo_rle <- rlFunc(day_by_offset)
  # # rename column #
  # names(dbo_rle)["value"] <- "day_by_offset"

  ### Calculate calendar date for start of each day-by-offset ##
  offset_start_date <- as.Date(dtimes[1])

  if(time_of_day[1] < hour_offset){
    offset_start_date <- offset_start_date - 1 # subtract a day if less than hour_offset (meaning start began previous day)
  }

  offset_date <- offset_start_date + (day_by_offset - 1) # add the day-by-offset

  # # merge to rle df #
  # dbo_rle$offset_date <- offset_start_date + (dbo_rle$day_by_offset - 1)

  # return results #
  return(data.frame(
    day_by_offset = day_by_offset,
    offset_date = offset_date
  ))
}


#' Convert numbers representing time-of-day to clock time in a string
#'
#' @param time_number A number between 0-24 (exclusive of 24) representing time of day.
#'
#' @returns A vector of length(time_number) representing clock time as a string.
#' @noRd
#'
numberToClockTime <- function(time_number){

  if(sum(time_number < 0 | time_number >= 24) > 0){
    stop("time_number must be a vector of numbers >= 0 and < 24")
  }

  ## extract time components of number ##
  hour_component <- floor(time_number) # hour component of number

  # decimal_component <- (time_number - hour_component) * 60 # minutes and seconds components
  # round to address floating point
  decimal_component <- time_number %% 1 * 60

  minutes_component <- floor(decimal_component) # minutes component

  seconds_component <- decimal_component %% 1 * 60 # seconds component

  ## check for rounding errors resulting in seconds or minutes equivalent to 60 ##
  # seconds #
  round_errors <- which(abs(60-seconds_component) < 1e-8)

  if(length(round_errors) > 0){
    minutes_component[round_errors] <- minutes_component[round_errors] + 1
    seconds_component[round_errors] <- 0
  }

  # minutes #
  round_errors2 <- which(abs(60 - minutes_component) < 1e-8)

  if(length(round_errors2) > 0){
    hour_component[round_errors2] <- hour_component[round_errors2] + 1
    minutes_component[round_errors2] <- 0
  }

  # hours #
  hour_component[hour_component == 24] <- 0

  ## convert to string ##
  time_string <- paste(sprintf("%02d", hour_component),
                       sprintf("%02d", minutes_component),
                       sprintf("%02f", seconds_component),
                       sep = ":")

  return(time_string)
}


dayCompleteness <- function(dtimes, hour_offset){

  ### Calculate day-by-offset values ###
  offset_df <- daySplit(dtimes = dtimes, hour_offset = hour_offset)

  ### Generate run data.frame ###
  ## generate day runs and grab start/end indices ##
  dbo_rle <- rlFunc(offset_df$day_by_offset)
  # rename column #
  names(dbo_rle)[which(names(dbo_rle)=="value")] <- "day_by_offset"
  # add start dates #
  dbo_rle$offset_date <- offset_df$offset_date[dbo_rle$start]

  ### Determine how complete each day is ##
  ## convert hour_offset to string ##
  hour_offset_string <- numberToClockTime(hour_offset)

  ## format start datetime for each day_by_offset ##
  offset_starts <- as.POSIXct(paste(dbo_rle$offset_date, hour_offset_string),
                              format = "%Y-%m-%d %H:%M:%S", tz = attr(dtimes, "tzone"))

  ## Calculate starting gaps ##
  dbo_rle$start_gaps <- as.numeric(difftime(dtimes[dbo_rle$start], offset_starts, units = "hours"))

  ## calculate ending gaps ##
  dbo_rle$end_gaps <- as.numeric(difftime(offset_starts + lubridate::days(1), dtimes[dbo_rle$end], units = "hours"))

  ## calculate durations between each time point ##
  dbo_rle$average_gap <- NA
  # dbo_rle$gap_duration <- NA
  for(i in 1:nrow(dbo_rle)){
    day_slice <- dtimes[dbo_rle$start[i]:dbo_rle$end[i]]
    # gaps between observations
    day_diff_times <- day_slice[2:length(day_slice)] - day_slice[1:(length(day_slice)-1)]
    dbo_rle$average_gap[i] <- mean(day_diff_times)

    # # calculate total sum
    # dbo_rle$gap_duration[i] <- sum(day_diff_times)
  }

  browser()


}


#' Convert POSIXct objects to time-of-day cumulative time (24-hour decimal format)
#'
#' @param dtime_vec A vector of POSIXct datetimes
#'
#' @returns A vector of numeric values indicating cumulative time, starting at first observed time-of-day.
#' @noRd
#'
ctimeCalc <- function(dtime_vec){

  ## check POSIXct format ##
  if(!is(dtime_vec, "POSIXct")){
    stop("dtime_vec must be a vector of datetime stamps in POSIXct format")
  }

  # start time #
  start_dtime <- dtime_vec[1]
  start_hour <- lubridate::hour(start_dtime)
  start_min <- lubridate::minute(start_dtime)
  start_sec <- lubridate::second(start_dtime)

  # convert to UNIX seconds #
  new_dtimes <- as.numeric(dtime_vec)

  # add to starting time (with starting time = time-of-day)
  new_dtimes2 <- (start_hour + start_min / 60 + start_sec / 60 / 60) + (new_dtimes - new_dtimes[1]) / 60 / 60

  return(new_dtimes2)
}


#' Prepare user-provided data.frame for use with rHCL functions
#'
#' @param df A data.frame.
#' @param time_var The name of the column in df that has date/time-stamps in
#' POSIXct format.
#' @param light_var The name of the column in df that has numeric values for
#' light. Default is NULL, meaning this column won't be checked/processed.
#' @param sleep_var The name of the column in df that has binary values
#' for sleep/wake states (0 = wake, 1 = sleep). Default is NULL, meaning this column won't be checked/processed.
#'
#' @returns A data.frame with the original time_var and light_var columns, as well
#' as a new processed "cumulative time" (ctime) column that represents time as
#' cumulative hours since the first observation, with the first observation represented
#' as time-of-day. If sleep_var was provided, it will also be added to the processed
#' data.frame.
#' @noRd
#'
dfPrep <- function(df, time_var, light_var=NULL, sleep_var = NULL){

  ### TODO consider adding imputation options for light_var and sleep_var

  ### Check df format ###
  # check that df is a data.frame #
  if(!is(df, "data.frame")){
    stop("df must be a data.frame")
  }

  ### Check timestamp format and convert to a continuous measure of time ###
  ## check that time_var is in df ##
  if(!time_var %in% names(df)){
    stop("time_var must be the name of a column in df")
  }

  ## check POSIXct format ##
  if(!is(df[[time_var]], "POSIXct")){
    stop("The time_var column in df must be a column of datetime stamps in POSIXct format")
  }

  if(sum(is.na(df[[time_var]])) > 0){
    stop("The time_var column in df must not have missing values")
  }

  ## Pre-process timestamps ##
  # convert to UNIX seconds #
  new_dtimes <- as.numeric(df[[time_var]])

  # check that there are no duplicate timestamps #
  if(sum(duplicated(new_dtimes)) > 0){
    stop("The time_var column in df must not have duplicate timestamps")
  }

  # check that time is always increasing #
  time_steps <- new_dtimes[2:length(new_dtimes)] - new_dtimes[1:(length(new_dtimes) - 1)]
  if(sum(time_steps <= 0 ) > 0){
    stop("The time_var column in df must have timestamps that are only increasing in time")
  }

  # ## Check resolution of data and provide warning if needed ##
  # if(mean(time_steps) / 60 > 10){
  #   warning(paste0("Average distance between timesteps is greater than 10 minutes. ",
  #                  "Accuracy of results may be decreased. Consider providing ",
  #                  "higher resolution data (e.g., interpolate 1-minute epochs)."))
  # }

  # create new time variable, calculated as hours since first observation, with
  # the first observation represented as time-of-day
  new_dtimes2 <- ctimeCalc(df[[time_var]])

  ### prep data.frame ###
  res_df <- data.frame(
    dtime = df[[time_var]],
    ctime = new_dtimes2
  )

  ### If light column is specified, check it and return ###
  if(!is.null(light_var)){

    ## Check light_var column ##
    ## Checks ##
    if(!light_var %in% names(df)){
      stop("light_var must be the name of a column in df")
    }

    if(!is(df[[light_var]], "numeric")){
      stop("The light_var column in df must be in numeric format")
    }

    if(sum(is.na(df[[light_var]])) > 0){
      stop("The light_var column in df must not have missing values")
    }

    if(sum(df[[light_var]] < 0) > 0){
      stop("The light_var column in df must not have negative values")
    }

    ## add light column to processed data.frame ##
    res_df[[light_var]] <- df[[light_var]]
  }


  ### If sleep column is specified, check it and return ###
  if(!is.null(sleep_var)){

    ## Checks ##
    if(!sleep_var %in% names(df)){
      stop("sleep_var must be the name of a column in df")
    }

    # check that sleep_var isn't missing data #
    if(sum(is.na(df[[sleep_var]])) > 0){
      stop("The sleep_var column in df must not have missing values")
    }

    # check that sleep_var is in correct format (either 0 or 1)
    if(sum(!unique(df[[sleep_var]] %in% c(0,1))) > 0){
      stop("sleep_var in df must be in binary format (0 = wake, 1 = sleep)")
    }

    ## add sleep column to processed data.frame ##
    res_df[[sleep_var]] <- df[[sleep_var]]

  }

  ### Return results ###
  return(res_df)

}
