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
})



