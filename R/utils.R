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


#' Run length encoding wrapper that identifies start and end indices
#'
#' @param vec Vector of values for which runs are to be extracted.
#'
#' @returns Data frame with run values, lengths, and start and end indices in original vector
#' @internal
#'
rlFunc <- function(vec){
  new_rle <- as.data.frame(unclass(rle(vec))) # convert rle to data frame
  new_rle$end <- cumsum(new_rle$lengths) # identifying end indices of runs
  new_rle$start <- new_rle$end - new_rle$lengths + 1 # identify start indices of runs

  # re-order and rename cols
  new_rle <- new_rle[,c("values", "lengths", "start", "end")]
  names(new_rle) <- c("value", "length", "start", "end")

  return(new_rle)
}

#' Shortest angle between two times in 24-hour decimal format.
#'
#' @param vec1 Vector of times in 24-hour decimal format
#' @param vec2 Second vector of times in 24-hour decimal format
#' @param lbound Lower bound of values. Defaults to -12, which indicates
#' possible phase angles are [-12, 12) (exclusive of 12).
#'
#' @returns Vector of phase angles in time. Positive values means the second time
#' occurs after the first time.
#' @internal
#'
clockAngle <- function(vec1, vec2, lbound = -12){
  diffs <- vec2 - vec1 - lbound # adjust the differences, with -12 as the lower bound
  diffs <- diffs %% 24 + lbound # modulo and center
  return(diffs)
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
#' if the times are equally dispersed around the 24 hours (e.g., midnight, 8, and 12).
#'
#'
#' @param vec Vector of times in 24-hour format.
#' @param na_rm Logical. If FALSE (default), will return NA if any NAs are in the vector.
#'
#' @returns The mean of the time vector in 24-hour decimal format
#' @internal
#'
timeMean <- function(vec, na_rm = FALSE){

  ## return NA if nothing in vector ##
  if(length(vec)==0){
    return(NA)
  }

  ## check for any NAs if na_rm = FALSE
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
    vec <- vec[!is.na(vec)]
  }

  ### convert times to radians ###
  rvec <- vec * pi / 12

  ### mean of radians ##
  part1 <- 1/length(rvec) * sum(sin(rvec)) # first part of atan2
  part2 <- 1/length(rvec) * sum(cos(rvec)) # second part of atan2

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
