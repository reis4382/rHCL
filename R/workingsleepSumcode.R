sleepSum <- function(x, sleep_var, time_var) {
  x[order(x[[time_var]]), ]

  sleep_first <- match("1", x[[sleep_var]])
  sleep_on <- x[sleep_first, "time"]

  sleep_last <- max(which(x[[sleep_var]] == 1)) + 1
  sleep_off <- x[sleep_last, "time"]

  sleep_dur <- sleep_off - sleep_on

  sleep_mid <- sleep_on + (sleep_dur/2)

  return(list(sleep_on = sleep_on, sleep_mid = sleep_mid, sleep_off = sleep_off, sleep_dur = sleep_dur))

}
