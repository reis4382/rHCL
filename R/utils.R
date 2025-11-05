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
