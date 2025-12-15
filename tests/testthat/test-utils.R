# Tests for minMaxFinder() ------------------------------------------------
test_that("minMaxFinder() finds minima and maxima", {

  vec <- c(1.1, 2.2, 4.4, 4.4, 2.2, -1.87, -4, -2, -2, -8.2, -8.2, -8.2, 0)

  expect_equal(minMaxFinder(y=vec, minima=TRUE), c(7, 10))
  expect_equal(minMaxFinder(y=vec, minima=FALSE), c(3, 8))
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


# Test rlFunc -------------------------------------------------------------

test_that("rlFunc() correctly identifies runs", {
  expect_equal(rlFunc(c(0,0,1,1,1,0,0,1)),
               data.frame(value = c(0, 1, 0, 1),
                          length = c(2, 3, 2, 1),
                          start = c(1, 3, 6, 8),
                          end = c(2, 5, 7, 8))
               )

})


# Test clockAngle calculations --------------------------------------------

test_that("clockAngle() returns the correct phase angles in 24-hour time", {

  expect_equal(clockAngle(vec1=c(4.5, 2.3, 6.2, 1.5, 0, 0, 0),
                          vec2=c(1.25, 4.6, 5.2, 23.5, 12, 11.5, 12.5)),
               c(-3.25, 2.3, -1, -2, -12, 11.5, -11.5))

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
