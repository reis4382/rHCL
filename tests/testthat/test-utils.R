
# rlFunc tests ------------------------------------------------------------

test_that("rlFunc_base() works correctly", {

  expect_equal(rlFunc_base(c(0, 0, 1, 1, 1)),
               data.frame("value" = c(0, 1),
                          "length" = c(2, 3),
                          "start" = c(1, 3),
                          "end" = c(2, 5)))

})

test_that("rlFunc_na() works correctly", {

  expect_equal(rlFunc_na(c(0, 0, 1, 0, 0, 1, 1, NA, NA, 1), na_collapse = TRUE),
               data.frame("value" = c(0, 1, 0, 1, NA, 1),
                          "length" = c(2, 1, 2, 2, 2, 1),
                          "start" = c(1, 3, 4, 6, 8, 10),
                          "end" = c(2, 3, 5, 7, 9, 10)))

})

test_that("rlFunc() works correctly", {

  expect_equal(rlFunc(c(0, 0, 1, 0, 0, 1, 1, NA, NA, 1), na_collapse = TRUE, gap_collapse = FALSE),
               data.frame("value" = c(0, 1, 0, 1, NA, 1),
                          "length" = c(2, 1, 2, 2, 2, 1),
                          "start" = c(1, 3, 4, 6, 8, 10),
                          "end" = c(2, 3, 5, 7, 9, 10)))

  expect_equal(rlFunc(c(0, 0, 1, 0, 0, 1, 1, NA, NA, 1), na_collapse = TRUE,
                      gap_collapse = TRUE, collapse_value = 0, max_gap = 1),
               data.frame("value" = c(0, 1, NA, 1),
                          "length" = c(5, 2, 2, 1),
                          "start" = c(1, 6, 8, 10),
                          "end" = c(5, 7, 9, 10)))

})

test_that("rlFunc() argument errors are caught", {

  expect_error(rlFunc(c(0, 1, 1), gap_collapse = TRUE), regexp = "If gap_collapse = TRUE, collapse_value cannot")
  expect_error(rlFunc(c(0, 1, 1), gap_collapse = TRUE, collapse_value = 3), regexp = "If gap_collapse = TRUE, collapse_value cannot")
  expect_error(rlFunc(c(0, 1, 1), gap_collapse = TRUE, collapse_value = c(0, 1)), regexp = "If gap_collapse = TRUE, collapse_value can only be a single")
  expect_error(rlFunc(c(0, 1, 1), gap_collapse = TRUE, collapse_value = 0), regexp = "If gap_collapse = TRUE, max_gap cannot be NULL")
  expect_error(rlFunc(c(0, 1, 1), gap_collapse = TRUE, collapse_value = 0, max_gap = -1), regexp = "If gap_collapse = TRUE, max_gap cannot be NULL")
  expect_error(rlFunc(c(0, 1, 1), gap_collapse = TRUE, collapse_value = 0, max_gap = "test"), regexp = "If gap_collapse = TRUE, max_gap cannot be NULL")

})


# Tests for minMaxFinder() ------------------------------------------------
test_that("minMaxFinder() finds minima and maxima", {

  vec <- c(1.1, 2.2, 4.4, 4.4, 2.2, -1.87, -4, -2, -2, -8.2, -8.2, -8.2, 0)

  expect_equal(minMaxFinder(y=vec, minima=TRUE), c(7, 10))
  expect_equal(minMaxFinder(y=vec, minima=FALSE), c(3, 8))
  expect_error(minMaxFinder(y=vec, minima="TEST"), regexp = "minima argument must be")
})


# Tests for sleepHomeostasis() --------------------------------------------
test_that("sleepHomeostasis() returns the expected values", {
  res1 <- sleepHomeostasis(mu=17.87, tswitch=0, h_tswitch=11, time=6, chi=45, s=0)
  expect_equal(res1, 17.87 + (11 - 17.87)*exp(-(6-0)/45))

  res2 <- sleepHomeostasis(mu=17.87, tswitch=0, h_tswitch=18, time=8, chi=45, s=1)
  expect_equal(res2, 18*exp(-(8-0)/45))
})

test_that("sleepHomeostasis() works with a vector of times", {
  times <- c(0, 6, 12, 18, 24)
  mu <- 17.87
  h_tswitch <- 10
  chi = 45
  s = 0

  res <- sleepHomeostasis(mu=mu, tswitch=0, h_tswitch=10, time=times, chi=chi, s=s)
  expect_equal(length(res), length(times))
})

test_that("sleepHomeostasis() returns the same values for different time scales and handles sleep/wake state switching", {
  hour_vec <- c(0, 4, 8, 12, 16, 20, 24) # time in hours
  min_vec <- hour_vec * 60 # time in minutes
  sleep_vec <- c(1, 1, 0, 0, 0, 1, 1) # sleep states

  mu <- 17.87
  start_h <- 16 # starting sleep pressure
  chi_hours <- 45 # chi for hours scale
  chi_minutes <- 45*60 # chi for minutes scale

  # separate sleep and wake runs, assuming person just fell asleep at beginning of data
  sleep_runs <- as.data.frame(unclass(rle(sleep_vec))) # create run data.frame
  sleep_runs$ends <- cumsum(sleep_runs$lengths) # ends of runs
  sleep_runs$starts <- sleep_runs$ends - sleep_runs$lengths + 1 # starts of runs

  # create data frame for storing data #
  res_df <- data.frame(sleep = sleep_runs$values,
                       hour_switch = hour_vec[sleep_runs$starts],
                       hour_end = hour_vec[sleep_runs$ends],
                       min_switch = min_vec[sleep_runs$starts],
                       min_end = min_vec[sleep_runs$ends]
                       )
  res_df$h_hour_start <- rep(NA, nrow(res_df)) # start of interval in hours
  res_df$h_min_start <- rep(NA, nrow(res_df)) # start of interval in minutes
  res_df[1,c("h_hour_start", "h_min_start")] <- start_h # starting sleep pressure

  res_df$h_hour_end <- rep(NA, nrow(res_df)) # end of interval in hours
  res_df$h_min_end <- rep(NA, nrow(res_df)) # end of interval in minutes

  # run results for hours and minutes time scales #
  for(i in 1:nrow(res_df)){
    # starting sleep pressure (in hours)
    if(is.na(res_df$h_hour_start[i])){
      res_df$h_hour_start[i] <- res_df$h_hour_end[i-1]
      h_tswitch_hour <- res_df$h_hour_end[i-1]
    } else{
      h_tswitch_hour <- res_df$h_hour_start[i]
    }

    # starting sleep pressure (in minutes)
    if(is.na(res_df$h_min_start[i])){
      res_df$h_min_start[i] <- res_df$h_min_end[i-1]
      h_tswitch_min <- res_df$h_min_end[i-1]
    } else{
      h_tswitch_min <- res_df$h_min_start[i]
    }

    # calculate ending h values
    h_hour <- sleepHomeostasis(mu=mu, tswitch=res_df$hour_switch[i],
                               h_tswitch = h_tswitch_hour, time = res_df$hour_end[i],
                               chi = chi_hours, s = res_df$sleep[i])
    res_df$h_hour_end[i] <- h_hour # replace h for end of interval (hours scale)

    h_min <- sleepHomeostasis(mu=mu, tswitch=res_df$min_switch[i],
                               h_tswitch = h_tswitch_min, time = res_df$min_end[i],
                               chi = chi_minutes, s = res_df$sleep[i])
    res_df$h_min_end[i] <- h_min # replace h for end of interval (minutes scale)
  }

  expect_equal(res_df$h_hour_end, res_df$h_min_end)
})

# angleDiff tests ---------------------------------------------------------

test_that("angleDiff() works", {

  vec1 <- c(90, 0, -10, 270) * pi / 180# degrees
  vec2 <- c(100, 192, 0, 90) * pi / 180 # degrees
  expected_rads <- c(10, -(360-192), 10, -180) * pi / 180

  expect_equal(angleDiffs(vec1, vec2, period = 2*pi, lbound = -pi), expected_rads)

})



# Test clockAngle calculations --------------------------------------------

test_that("clockAngle() returns the correct phase angles in 24-hour time", {

  expect_equal(clockAngle(vec1=c(4.5, 2.3, 6.2, 1.5, 0, 0, 0),
                          vec2=c(1.25, 4.6, 5.2, 23.5, 12, 11.5, 12.5)),
               c(-3.25, 2.3, -1, -2, -12, 11.5, -11.5))

})


# Test timeToTOD function -------------------------------------------------

test_that("timeToTOD() works", {

  vals <- lubridate::ymd_hms(c("2025-01-01 00:00:00", "2025-04-01 05:11:00",
                               "2025-07-07 18:03:02", "2025-03-25 23:59:59"), tz = "America/Denver")
  expect_equal(timeToTOD(vals), c(0, 5+11/60, 18 + 3/60 + 2/60/60, 23 + 59/60 + 59/60/60))

})

# Test timeMean calculations ----------------------------------------------

test_that("timeMean() correctly calculates the mean of times", {
  expect_equal(timeMean(c(23, 22, 1, 2)), 0) # ensure midnight returns as 0
  expect_equal(timeMean(c(1, 2, 4, 5)), 3)
  expect_equal(timeMean(c(23, 3, 7)), 3)
})

test_that("timeMean() returns warnings for undefined times",{
  expect_warning(timeMean(c(0, 12)), regexp = "Mean of times is undefined")
  expect_warning(timeMean(c(0, 8, 16)), regexp = "Mean of times is undefined")
  expect_warning(timeMean(c(1, 9, 17)), regexp = "Mean of times is undefined")
  expect_warning(timeMean(c(3, 9, 15, 21)), regexp = "Mean of times is undefined")
})

test_that("timeMean() returns NA if empty vector provided as input", {
  expect_equal(timeMean(c()), NA)
})

test_that("timeMean() correctly weights angles", {
  expect_equal(round(timeMean(c(4, 5, 13), weights = c(8, 8, 1)), 2), 4.70)
  expect_equal(timeMean(c(0, 5), weights = c(4, 4)), 2.5)
  expect_equal(round(timeMean(c(0, 4), weights = c(6, 2)), 2), 0.93)
})

test_that("timeMean() returns error if weights are provided and not numeric",{
  expect_error(timeMean(c(23, 22, 1), weights = c("a", "b", "c")), regexp = "Weights provided to timeMean\\(\\) function must be numeric")
})

test_that("timeMean() returns error if weight vector is not the same length as data vector", {
  expect_error(timeMean(c(23, 22, 1), weights = c(1, 2)), regexp = "Vector of weights must be the same length as input \\(timeMean\\(\\)\\)")
})

test_that("timeMean() handles NA removal", {

  expect_equal(timeMean(c(23, 22, 1, 2, NA), na_rm = TRUE), 0)
  expect_equal(timeMean(c(23, NA, NA, NA, NA), na_rm = TRUE), NA)
})

# Tests for lightCycle function -------------------------------------------

test_that("lightCycle() correctly generates a light profile", {

  reps <- 300 # number of days of data
  l1 <- 700 # max light
  l2 <- 40 # min light
  times <- seq(0, 24*reps, by = .2)
  light <- lightCycle(t = times, time_scale = "hours")

  # ensure all max light values are centered at noon
  expect_equal(times[light==max(light)], 12 + 24*(0:(reps-1)))
  # ensure (rounded) max and min values are close to l1 and l2
  expect_equal(abs(round(max(light))-l1) < 10, TRUE)
  expect_equal(abs(round(min(light))-l2) < 1, TRUE)

  ## check seconds time scale ##
  times2 <- times * 60*60 # seconds
  light2 <- lightCycle(t = times2, time_scale = "secs")

  # ensure light values are identical #
  expect_equal(light, light2)

})


# Tests for calculating noon-to-noon or midnight-to-midnight days ---------

test_that("epochDays() correctly calculates noon-to-noon or midnight-to-midnight days", {

  expect_equal(epochDays(c(0, 3, 15, 26, 37, 48), noon_to_noon = TRUE), c(1, 1, 2, 2, 3, 3))
  expect_equal(epochDays(c(0, 3, 15, 26, 37, 48), noon_to_noon = FALSE), c(1, 1, 1, 2, 2, 3))
  expect_equal(epochDays(c(36, 37, 60, 72, 73), noon_to_noon = TRUE), c(1, 1, 2, 2, 2))

})


# Tests for daybyoffset functions -----------------------------------------

test_that("dayByOffsetVector() works", {

  dtimes <- seq(from=as.POSIXct("2025-05-20 17:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                to=as.POSIXct("2025-05-21 13:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                by = "hours")

  res <- dayByOffsetVector(dtimes, hour_offset = 12) # noon-to-noon days

  expect_equal(res, c(rep(1, 19), rep(2,2)))

  res2 <- dayByOffsetVector(dtimes, hour_offset = 18) # noon-to-noon days

  expect_equal(res2, c(1, rep(2, 20)))

})

test_that("dayByOffsetVector() works with dtimes from a different TZ", {

  dtimes <- seq(from=as.POSIXct("2025-05-20 17:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/New_York"),
                to=as.POSIXct("2025-05-21 13:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/New_York"),
                by = "hours")

  res <- dayByOffsetVector(dtimes, hour_offset = 12) # noon-to-noon days

  expect_equal(res, c(rep(1, 19), rep(2,2)))

  res2 <- dayByOffsetVector(dtimes, hour_offset = 18) # noon-to-noon days

  expect_equal(res2, c(1, rep(2, 20)))

})

test_that("offsetDates() works", {

  start_time <- as.POSIXct("2025-01-01 11:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  end_time <- as.POSIXct("2025-01-01 18:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

  dtime_vec <- seq(start_time, end_time, by = "1 hour")

  res <- offsetDates(dtime_vec, hour_offset = 12)
  res2 <- offsetDates(dtime_vec, hour_offset = 18)

  expect_equal(res, c(as.Date("2024-12-31"), rep(as.Date("2025-01-01"), 7)))
  expect_equal(res2, c(rep(as.Date("2024-12-31"), 7), as.Date("2025-01-01")))

})

test_that("offsetDates() works with different time zones", {
  start_time <- as.POSIXct("2025-01-01 11:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/New_York")
  end_time <- as.POSIXct("2025-01-01 23:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/New_York")

  dtime_vec <- seq(start_time, end_time, by = "1 hour")

  res <- offsetDates(dtime_vec, hour_offset = 12)
  res2 <- offsetDates(dtime_vec, hour_offset = 18)

  expect_equal(res, c(as.Date("2024-12-31"), rep(as.Date("2025-01-01"), 12)))
  expect_equal(res2, c(rep(as.Date("2024-12-31"), 7), rep(as.Date("2025-01-01"), 6)))


})


# Tests for splitting days using POSIXct variables ------------------------

test_that("numberToClockTime() works", {

  times <- c(1.5, 5.75, 13 + 42/60, 22 + (12.5/60))
  expected <- c("01:30:0.000000", "05:45:0.000000", "13:42:0.000000", "22:12:30.000000")

  expect_equal(numberToClockTime(times), expected)

  ## check for error ##
  expect_error(numberToClockTime(25), regexp = "time_number must be a vector")

})


# Tests for ctimeCalc() ---------------------------------------------------

test_that("ctimeCalc() works", {

  dtime_vec <- c("2025-01-01 00:00:00", "2025-01-01 10:30:00", "2025-01-01 23:15:00",
                 "2025-01-02 12:00:00", "2025-01-04 18:00:00")

  dtime_vec <- as.POSIXct(dtime_vec, format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")

  res <- ctimeCalc(dtime_vec)

  expect_equal(res, c(0, 10.5, 23.25, 36, 90))

  expect_error(ctimeCalc("hi"), regexp = "dtime_vec must be a vector of datetime stamps in POSIXct format")

})


# Tests for preparing data.frame ------------------------------------------
test_that("dfPrep() works", {

  start_dtime <- as.POSIXct("2025-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_dtime <- as.POSIXct("2025-01-04 10:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  dtimes <- seq(start_dtime, end_dtime, by = "1 hour")

  df <- data.frame(
    dtime = dtimes,
    lux = abs(rnorm(length(dtimes))),
    sleep = rbinom(length(dtimes), size = 1, prob = .3)
  )

  res1 <- dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = NULL)

  expect_equal(res1[,c("dtime", "lux")], df[,c("dtime", "lux")])

  res2 <- dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = "sleep")
  expect_equal(res2[,c("dtime", "lux", "sleep")], df[,c("dtime", "lux", "sleep")])

  expect_error(dfPrep(df = 5, time_var = "dtime", light_var = "lux"))

})

test_that("dfPrep() catches errors in timestamp formatting", {

  start_dtime <- as.POSIXct("2025-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_dtime <- as.POSIXct("2025-01-04 10:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  dtimes <- seq(start_dtime, end_dtime, by = "1 hour")

  df <- data.frame(
    dtime = dtimes,
    lux = abs(rnorm(length(dtimes))),
    sleep = rbinom(length(dtimes), size = 1, prob = .3)
  )

  expect_error(dfPrep(df = df, time_var = "dtime_col", light_var = "lux", sleep_var = NULL),
               regexp = "time_var must be the name of a column in df")

  df$dtime <- "HI"
  expect_error(dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = NULL),
               regexp = "The time_var column in df must be a column of datetime stamps")

  df$dtime <- dtimes
  df$dtime[1] <- NA
  expect_error(dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = NULL),
               regexp = "The time_var column in df must not have missing values")

  df$dtime[1] <- df$dtime[2]
  expect_error(dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = NULL),
               regexp = "The time_var column in df must not have duplicate timestamps")

  df$dtime[1] <- dtimes[1]
  df$dtime[5] <- df$dtime[4] - lubridate::minutes(1)

  expect_error(dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = NULL),
               regexp = "The time_var column in df must have timestamps that are only increasing in time")

})

test_that("dfPrep() catches non-timestamp errors", {

  start_dtime <- as.POSIXct("2025-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_dtime <- as.POSIXct("2025-01-04 10:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  dtimes <- seq(start_dtime, end_dtime, by = "1 hour")
  lux_vals <- abs(rnorm(length(dtimes)))
  sleep_vals <- rbinom(length(dtimes), size = 1, prob = .3)

  df <- data.frame(
    dtime = dtimes,
    lux = lux_vals,
    sleep = sleep_vals
  )

  expect_error(dfPrep(df = df, time_var = "dtime", light_var = "lux_col", sleep_var = "sleep"),
               regexp = "light_var must be the name of a column in df")

  df$lux <- "HI"
  expect_error(dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = NULL),
               regexp = "The light_var column in df must be in numeric format")

  df$lux <- lux_vals

  df$lux[1] <- NA
  expect_error(dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = NULL),
               regexp = "The light_var column in df must not have missing values")

  df$lux[1] <- -50
  expect_error(dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = NULL),
               regexp = "The light_var column in df must not have negative values")

  df$lux <- lux_vals
  expect_error(dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = "sleep_col"),
               regexp = "sleep_var must be the name of a column in df")

  df$sleep <- "HI"
  expect_error(dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = "sleep"),
               regexp = "sleep_var in df must be in binary format")

  df$sleep <- sleep_vals

  df$sleep[2] <- NA
  expect_error(dfPrep(df = df, time_var = "dtime", light_var = "lux", sleep_var = "sleep"),
               regexp = "The sleep_var column in df must not have missing values")
})

