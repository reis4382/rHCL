test_that("sleepRunSummary() works", {
  start_time <- as.POSIXct("2025-01-01 12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_time <- as.POSIXct("2025-01-03 11:59:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")

  s_df <- data.frame("time"=seq(from = start_time, to = end_time, by = "1 min"),
                     "S" = 0)

  ctimes <- ctimeCalc(s_df$time) # cumulative time of day to assist test data prep

  ## set sleep runs ##
  s_df[ctimes >= 22 & ctimes < 30, "S"] <- 1 # first run
  s_df[ctimes >= 46 & ctimes < 53, "S"] <- 1 # second run

  res <- sleepRunSummary(df = s_df, sleep_var = "S", time_var = "time",
                         epoch_length_min = 1)

  sleep_onset <- lubridate::ymd_hms(c("2025-01-01 22:00:00", "2025-01-02 22:00:00"), tz = "America/Denver")
  sleep_midpoint <- lubridate::ymd_hms(c("2025-01-02 02:00:00", "2025-01-03 01:30:00"), tz = "America/Denver")
  sleep_offset <- lubridate::ymd_hms(c("2025-01-02 06:00:00", "2025-01-03 05:00:00"), tz = "America/Denver")
  sleep_duration <- c(8, 7)

  expected_df <- data.frame(sleep_onset = sleep_onset,
                            sleep_midpoint = sleep_midpoint,
                            sleep_offset = sleep_offset,
                            sleep_duration = sleep_duration)

  expect_equal(res[,c("sleep_onset", "sleep_midpoint", "sleep_offset", "sleep_duration")], expected_df)

})


test_that("sleep24Summary() correctly calculates 24-hour metrics", {
  start_time <- as.POSIXct("2025-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_time <- as.POSIXct("2025-01-04 11:54:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")

  s_df <- data.frame("time"=seq(from = start_time, to = end_time, by = "6 min"),
                     "S" = 0)

  ctimes <- ctimeCalc(s_df$time) # cumulative time of day to assist test data prep

  ## set sleep runs ##
  s_df[ctimes >= 0 & ctimes < 7, "S"] <- 1 # first run
  s_df[ctimes >= 22 & ctimes < 29.5, "S"] <- 1
  s_df[ctimes >= 47 & ctimes < 55, "S"] <- 1
  s_df[ctimes >= 72 & ctimes < 80, "S"] <- 1

  # noon-to-noon parsing
  expect_equal(sleep24Summary(df=s_df, sleep_var = "S", time_var = "time", epoch_length_min = 6, hour_offset = 12),
  data.frame(
    offset_date = as.Date(c("2024-12-31", "2025-01-01", "2025-01-02", "2025-01-03")),
    offset_day = c(1, 2, 3, 4),
             hour_offset = rep(12, 4),
             observed_hours = c(12, 24, 24, 24),
             sleep_duration = c(7, 7.5, 8, 8),
             sleep_midpoint = c(7/2, 25.75%%24, 51%%24, 76%%24)))

  # midnight-to-midnight parsing
  m2m_df <- sleep24Summary(df=s_df, sleep_var = "S", time_var = "time", epoch_length_min = 6, hour_offset = 0)
  m2m_df$sleep_midpoint <- round(m2m_df$sleep_midpoint, 2) # round midpoint values

  expect_equal(m2m_df,
               data.frame(
                 offset_date = as.Date(c("2025-01-01", "2025-01-02", "2025-01-03", "2025-01-04")),
                 offset_day = c(1, 2, 3, 4),
                 hour_offset = rep(0, 4),
                 observed_hours = c(24, 24, 24, 12),
                 sleep_duration = c(9, 6.5, 7, 8),
                 sleep_midpoint = c(2.61, 2.29, 3.5, 4)))

})

test_that("sleep24Summary() handles runs without sleep", {
  start_time <- as.POSIXct("2025-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_time <- as.POSIXct("2025-01-01 23:59:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")

  s_df <- data.frame("time"=seq(from = start_time, to = end_time, by = "6 min"),
                     "S" = 0)

  res <- sleep24Summary(df = s_df, sleep_var = "S", time_var = "time", epoch_length = 6, hour_offset = 0)

  expect_equal(nrow(res), 1)
  expect_equal(res$sleep_midpoint, NA)
  expect_equal(res$sleep_duration, 0)
})


test_that("sleepSummary() correctly summarizes sleep runs", {

  ## simulated data ##
  start_time <- as.POSIXct("2025-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_time <- as.POSIXct("2025-01-07 23:59:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")

  s_df <- data.frame("time"=seq(from = start_time, to = end_time, by = "6 min"),
                     "S" = 0)

  ctimes <- ctimeCalc(s_df$time) # cumulative time of day to assist test data prep

  ## set sleep runs ##
  s_df[(ctimes %% 24) >= 22 | (ctimes %% 24) < 6, "S"] <- 1

  n2n_df1 <- sleep24Summary(s_df, sleep_var = "S", time_var = "time",
                            epoch_length_min = 6, hour_offset = 12)
  n2n_df1 <- n2n_df1[n2n_df1$observed_hours >= 18,]
  row.names(n2n_df1) <- 1:nrow(n2n_df1)

  m2m_df1 <- sleep24Summary(s_df, sleep_var = "S", time_var = "time",
                            epoch_length_min = 6, hour_offset = 0)
  m2m_df1 <- m2m_df1[m2m_df1$observed_hours >= 18,]
  row.names(m2m_df1) <- 1:nrow(m2m_df1)


  ## prep expected results ##
  sleep_on_inds <- c(1, which(ctimes%%24 == 22))
  sleep_off_inds <- c(which(ctimes%%24 == 6)-1, nrow(s_df))

  sleep_onsets <- s_df$time[sleep_on_inds]
  sleep_offsets <- s_df$time[sleep_off_inds] + 6 * 60 # add one epoch, as sleep offset ends after final row of sleep run

  sleep_durations <- as.numeric(difftime(sleep_offsets, sleep_onsets, units = "hours"))
  sleep_midpoints <- sleep_onsets + sleep_durations / 2 * 60 * 60


  ## dropping any sleep run that hits the start or end of the data
  expect_equal(sleepSummary(df=s_df, sleep_var="S", time_var="time",
                            epoch_length_min = 6, min_observed_hours = 18),
               list(
                 summary = data.frame(
                   sleep_mid = timeMean(c(3, 2, 2, 2, 2, 2, 2, 23), c(6, 8, 8, 8, 8, 8, 8, 2)),
                   sleep_dur_noon_24hr = mean(n2n_df1$sleep_duration),
                   sleep_dur_midnight_24hr = mean(m2m_df1$sleep_duration)
                 ),
                 sleep_runs = data.frame(
                   sleep_onset = sleep_onsets,
                   sleep_midpoint = sleep_midpoints,
                   sleep_offset = sleep_offsets,
                   sleep_duration = sleep_durations,
                   start_index = sleep_on_inds,
                   end_index = sleep_off_inds
                 ),
                 noon_to_noon = n2n_df1,
                 midnight_to_midnight = m2m_df1
               )

  )

  ### Varied sleep windows ###
  ## simulated data ##
  start_time2 <- as.POSIXct("2025-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_time2 <- as.POSIXct("2025-01-03 12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")

  s_df2 <- data.frame("time"=seq(from = start_time2, to = end_time2, by = "6 min"),
                     "S" = 0)

  ctimes2 <- ctimeCalc(s_df2$time) # cumulative time of day to assist test data prep


  ## set sleep runs ##
  s_df2[ctimes2 >= 0 & ctimes2 < 6, "S"] <- 1 # sleep run # 1
  s_df2[ctimes2 >= 22.5 & ctimes2 < 30, "S"] <- 1 # sleep run # 2
  s_df2[ctimes2 >= 47 & ctimes2 < 58.3, "S"] <- 1 # sleep run # 3

  n2n_df2 <- sleep24Summary(s_df2, sleep_var = "S", time_var = "time",
                            epoch_length_min = 6, hour_offset = 12)
  n2n_df2 <- n2n_df2[n2n_df2$observed_hours >= 18,]
  row.names(n2n_df2) <- 1:nrow(n2n_df2)

  m2m_df2 <- sleep24Summary(s_df2, sleep_var = "S", time_var = "time",
                            epoch_length_min = 6, hour_offset = 0)
  m2m_df2 <- m2m_df2[m2m_df2$observed_hours >= 18,]
  row.names(m2m_df2) <- 1:nrow(m2m_df2)

  ## prep expected results ##
  sleep_on_inds2 <- which(ctimes2 %in% c(0, 22.5, 47))
  sleep_off_inds2 <- which(ctimes2 %in% c(6, 30, 58.3)) - 1 # subtract an index

  sleep_onsets2 <- s_df2$time[sleep_on_inds2]
  sleep_offsets2 <- s_df2$time[sleep_off_inds2] + 6 * 60 # add one epoch, as sleep offset ends after final row of sleep run

  sleep_durations2 <- as.numeric(difftime(sleep_offsets2, sleep_onsets2, units = "hours"))
  sleep_midpoints2 <- sleep_onsets2 + sleep_durations2 / 2 * 60 * 60

  res <- sleepSummary(df=s_df2, sleep_var="S", time_var="time",
                      epoch_length_min = 6, min_observed_hours = 18)

  expect_equal(
    res,
    # expected data
    list(
      summary = data.frame(
        sleep_mid = timeMean(c(3, 26.25, 52.65), c(6, 7.5, 11.3)),
        sleep_dur_noon_24hr = mean(n2n_df2$sleep_duration),
        sleep_dur_midnight_24hr = mean(m2m_df2$sleep_duration)
      ),
      sleep_runs = data.frame(
        sleep_onset = sleep_onsets2,
        sleep_midpoint = sleep_midpoints2,
        sleep_offset = sleep_offsets2,
        sleep_duration = sleep_durations2,
        start_index = sleep_on_inds2,
        end_index = sleep_off_inds2
      ),
      noon_to_noon = n2n_df2,
      midnight_to_midnight = m2m_df2
    )
  )

})

test_that("sleepSummary() correctly handles data with no sleep", {

  ## simulated data ##
  start_time2 <- as.POSIXct("2025-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_time2 <- as.POSIXct("2025-01-03 12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")

  s_df2 <- data.frame("time"=seq(from = start_time2, to = end_time2, by = "6 min"),
                      "S" = 0)

  res <- sleepSummary(df=s_df2, sleep_var="S", time_var="time",
                      epoch_length_min = 6, min_observed_hours = 18)

  expect_equal(res$summary$sleep_mid, NA)
  expect_equal(res$summary$sleep_dur_noon_24hr, 0)
  expect_equal(res$summary$sleep_dur_midnight_24hr, 0)

})

test_that("sleepSummary() handles short datasets less than a day without sleep", {

  ## simulated data ##
  start_time <- as.POSIXct("2025-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_time <- as.POSIXct("2025-01-01 00:03:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")

  s_df <- data.frame("time"=seq(from = start_time, to = end_time, by = "1 min"),
                      "S" = 0)

  res <- sleepSummary(df=s_df, sleep_var="S", time_var="time",
                      epoch_length_min = 1, min_observed_hours = 18)

  expect_equal(res$summary$sleep_mid, NA)
  expect_equal(res$summary$sleep_dur_noon_24hr, NA)
  expect_equal(res$summary$sleep_dur_midnight_24hr, NA)

})

test_that("sleepProcessQuick() works", {


  ## simulated data ##
  start_time <- as.POSIXct("2025-01-01 12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_time <- as.POSIXct("2025-01-02 11:59:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")

  df1 <- data.frame(
    dtime = seq(start_time, end_time, by = "1 min"),
    sleep = 0
  )

  df1$sleep[lubridate::hour(df1$dtime) < 6 | lubridate::hour(df1$dtime) >= 22] <- 1

  # drop a few indices #
  df1 <- df1[-c(100, 586, 1000), ]

  res1 <- sleepProcessQuick(df1, sleep_var = "sleep", time_var = "dtime", epoch_length_min = 1,
                            min_observed_hours = 18)

  expect_equal(abs(res1[["sleep_duration"]] - 8) < .02, TRUE)
  expect_equal(abs(res1[["sleep_midpoint"]] - 2) < .01, TRUE)

})

test_that("sleepProcessQuick() handles epochs of different lengths", {

  df1 <- data.frame(
    dtime = as.POSIXct(c("2025-01-01 12:00:00", "2025-01-01 12:00:30", "2025-01-01 12:01:00"), format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver"),
    sleep = c(0, 1, 0)
  )

  res1 <- sleepProcessQuick(df1, sleep_var = "sleep", time_var = "dtime", epoch_length_min = 0.5,
                            min_observed_hours = 18)

  # w/o 24 hours of data, sleep duration will be NA given the current code
  # note - midpoints are in the middle of a given epoch
  expect_equal(res1$sleep_midpoint, 12 + ((30/60) + (30 / 60 / 2))/60)


  df2 <- data.frame(
    dtime = as.POSIXct(c("2025-01-01 12:00:00", "2025-01-01 12:01:31", "2025-01-01 12:03:02"), format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver"),
    sleep = c(0, 1, 0)
  )

  res2 <- suppressWarnings(sleepProcessQuick(df2, sleep_var = "sleep", time_var = "dtime", epoch_length_min = 91/60))

  # w/o 24 hours of data, sleep duration will be NA given the current code
  expect_equal(res2$sleep_midpoint, 12 + ((91/60) + (91 / 60 / 2))/60)

  # expect_warning(sleepProcessQuick(df2, sleep_var = "sleep", time_var = "dtime", epoch_length_min = 91/60,
  #                                  min_observed_hours = 18),
  #                regex = "Epoch_length_min is not an factor of 60 minutes")


})

test_that("sleepDurationQuick2() works", {

  ## Easily classified "average sleep per 24 hours"
  start_dtime1 <- as.POSIXct("2025-01-01 12:00:00")

  ### one day with sleep between 10-6 am ###
  ctimes1 <- (0:(60*24*2-1) * 30) # seconds of cumulative time
  times1 <- start_dtime1 + ctimes1


  sleep1 <- rep(0, length(times1))
  sleep1[(ctimes1 / 60 / 60 + 12) %% 24 >= 22 | (ctimes1 / 60 / 60 + 12) %% 24 < 6] <- 1

  res1 <- sleepDurationQuick2(df = data.frame(dtime = times1, sleep = sleep1),
                             time_var = "dtime", sleep_var = "sleep", epoch_length_min = 0.5,
                             min_observed_hours = 18)

  expect_equal(8, res1)

  ### two days with sleep between 10-6 am ###
  ctimes2 <- (0:(60*24*2*2-1) * 30) # seconds of cumulative time
  times2 <- start_dtime1 + ctimes2


  sleep2 <- rep(0, length(times2))
  sleep2[(ctimes2 / 60 / 60 + 12) %% 24 >= 22 | (ctimes2 / 60 / 60 + 12) %% 24 < 6] <- 1

  res2 <- sleepDurationQuick2(df = data.frame(dtime = times2, sleep = sleep2),
                             time_var = "dtime", sleep_var = "sleep", epoch_length_min = 0.5,
                             min_observed_hours = 18)

  expect_equal(8, res2)


  ### two days with sleeps of 7 and 8 hours ###
  ctimes3 <- (0:(60*24*2*2-1) * 30) # seconds of cumulative time
  times3 <- start_dtime1 + ctimes3

  sleep3 <- rep(0, length(times3))
  sleep3[(ctimes3 / 60 / 60 + 12) >= 23 & (ctimes3 / 60 / 60 + 12) < 30] <- 1
  sleep3[(ctimes3 / 60 / 60 + 12) >= 47 & (ctimes3 / 60 / 60 + 12) < 55] <- 1

  res3 <- sleepDurationQuick2(df = data.frame(dtime = times3, sleep = sleep3),
                             time_var = "dtime", sleep_var = "sleep", epoch_length_min = 0.5,
                             min_observed_hours = 18)

  expect_equal(7.5, res3)

  ### two days with sleeps of 8 and 7 hours - minute epochs ###
  ctimes4 <- (0:(60*24*2-1) * 60) # seconds of cumulative time
  times4 <- start_dtime1 + ctimes4

  sleep4 <- rep(0, length(times4))
  sleep4[(ctimes4 / 60 / 60 + 12) >= 22 & (ctimes4 / 60 / 60 + 12) < 30] <- 1
  sleep4[(ctimes4 / 60 / 60 + 12) >= 47 & (ctimes4 / 60 / 60 + 12) < 54] <- 1

  res4 <- sleepDurationQuick2(df = data.frame(dtime = times4, sleep = sleep4),
                             time_var = "dtime", sleep_var = "sleep", epoch_length_min = 1,
                             min_observed_hours = 18)

  expect_equal(7.5, res4)

  ### 2 days with random sleep  - minute epochs ###
  ctimes5 <- (0:(60*24*2-1) * 60) # seconds of cumulative time
  times5 <- start_dtime1 + ctimes5

  sleep5 <- sample(0:1, size = length(times5), replace = T)

  res5 <- sleepDurationQuick2(df = data.frame(dtime = times5, sleep = sleep5),
                             time_var = "dtime", sleep_var = "sleep", epoch_length_min = 1,
                             min_observed_hours = 18)

  expect_equal(mean(sleep5)*24, res5)

  ### 5 days with sleeps of 9 (x3) and 7 (x2) hours - minute epochs ###
  ctimes6 <- (0:(60*24*5-1) * 60) # seconds of cumulative time
  times6 <- start_dtime1 + ctimes6

  sleep6 <- rep(0, length(times6))
  sleep6[(ctimes6 / 60 / 60 + 12) >= (21+24*0) & (ctimes6 / 60 / 60 + 12) < (30+24*0)] <- 1
  sleep6[(ctimes6 / 60 / 60 + 12) >= (21+24*1) & (ctimes6 / 60 / 60 + 12) < (30+24*1)] <- 1
  sleep6[(ctimes6 / 60 / 60 + 12) >= (21+24*2) & (ctimes6 / 60 / 60 + 12) < (30+24*2)] <- 1
  sleep6[(ctimes6 / 60 / 60 + 12) >= (23+24*3) & (ctimes6 / 60 / 60 + 12) < (30+24*3)] <- 1
  sleep6[(ctimes6 / 60 / 60 + 12) >= (23+24*4) & (ctimes6 / 60 / 60 + 12) < (30+24*4)] <- 1


  res6 <- sleepDurationQuick2(df = data.frame(dtime = times6, sleep = sleep6),
                             time_var = "dtime", sleep_var = "sleep", epoch_length_min = 1,
                             min_observed_hours = 18)

  expect_equal(((9*3)+(7*2)) / 5, res6)


  ### More than an integer number of days ###
  ## Fifth day goes past noon, but still same amount of sleep
  ctimes7 <- (0:(60*24*5+200) * 60) # seconds of cumulative time
  times7 <- start_dtime1 + ctimes7

  sleep7 <- rep(0, length(times7))
  sleep7[(ctimes7 / 60 / 60 + 12) >= (21+24*0) & (ctimes7 / 60 / 60 + 12) < (30+24*0)] <- 1
  sleep7[(ctimes7 / 60 / 60 + 12) >= (21+24*1) & (ctimes7 / 60 / 60 + 12) < (30+24*1)] <- 1
  sleep7[(ctimes7 / 60 / 60 + 12) >= (21+24*2) & (ctimes7 / 60 / 60 + 12) < (30+24*2)] <- 1
  sleep7[(ctimes7 / 60 / 60 + 12) >= (23+24*3) & (ctimes7 / 60 / 60 + 12) < (30+24*3)] <- 1
  sleep7[(ctimes7 / 60 / 60 + 12) >= (23+24*4) & (ctimes7 / 60 / 60 + 12) < (30+24*4)] <- 1


  res7 <- sleepDurationQuick2(df = data.frame(dtime = times7, sleep = sleep7),
                             time_var = "dtime", sleep_var = "sleep", epoch_length_min = 1,
                             min_observed_hours = 18)

  expect_equal(((9*3)+(7*2)) / 5, res7)

  ### Sleep crosses day "boundaries" with variable sleep times ###
  ctimes8 <- (0:(60*24*5+200) * 60) # seconds of cumulative time
  times8 <- start_dtime1 + ctimes8 # start at noon

  sleep8 <- rep(0, length(times8))
  # Day 1: 1 hour nap + 6.5 hour main
  sleep8[(ctimes8 / 60 / 60 + 12) >= (13+24*0) & (ctimes8 / 60 / 60 + 12) < (14+24*0)] <- 1 # 1 hour nap
  sleep8[(ctimes8 / 60 / 60 + 12) >= (23.5+24*0) & (ctimes8 / 60 / 60 + 12) < (30+24*0)] <- 1 # 6.5 hour sleep, ending in morning

  # Day 2: 0.5 hour nap + 7 hour sleep ending at 1 pm next day
  sleep8[(ctimes8 / 60 / 60 + 12) >= (13.5+24*1) & (ctimes8 / 60 / 60 + 12) < (14+24*1)] <- 1 # 0.5 hour nap
  sleep8[(ctimes8/ 60 / 60 + 12) >= (6+24*2) & (ctimes8 / 60 / 60 + 12) < (13+24*2)] <- 1 # 7 hour late sleep

  # Day 3: 6 hour main sleep
  sleep8[(ctimes8 / 60 / 60 + 12) >= (23.5+24*2) & (ctimes8 / 60 / 60 + 12) < (5.5+24*3)] <- 1

  # Day 4: 3 3-hour naps/sleeps
  sleep8[(ctimes8 / 60 / 60 + 12) >= (14+24*3) & (ctimes8 / 60 / 60 + 12) < (17+24*3)] <- 1
  sleep8[(ctimes8 / 60 / 60 + 12) >= (21+24*3) & (ctimes8 / 60 / 60 + 12) < (24+24*3)] <- 1
  sleep8[(ctimes8 / 60 / 60 + 12) >= (2+24*4) & (ctimes8 / 60 / 60 + 12) < (5+24*4)] <- 1

  # Day 5: 10 hour sleep
  sleep8[(ctimes8 / 60 / 60 + 12) >= (20+24*4) & (ctimes8 / 60 / 60 + 12) < (6+24*5)] <- 1

  res8 <- sleepDurationQuick2(df = data.frame(dtime = times8, sleep = sleep8),
                             time_var = "dtime", sleep_var = "sleep", epoch_length_min = 1,
                             min_observed_hours = 18)

  expect_equal(abs(mean(c(7.5, 7.5, 6, 9, 10)) - res8)*60 < 10, TRUE)


})
