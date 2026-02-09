test_that("sleep24Summary() correctly calculates 24-hour metrics", {
  s_df <- data.frame("time"=seq(from = 0, to = 83.9, by = .1),
                     "S" = 0)

  ## set sleep runs ##
  s_df[s_df$time >= 0 & s_df$time < 7, "S"] <- 1 # first run
  s_df[s_df$time >= 22 & s_df$time < 29.5, "S"] <- 1
  s_df[s_df$time >= 47 & s_df$time < 55, "S"] <- 1
  s_df[s_df$time >= 72 & s_df$time < 80, "S"] <- 1

  # noon-to-noon parsing
  expect_equal(sleep24Summary(df=s_df, sleep_var = "S", time_var = "time", epoch_length = .1, noon_to_noon = TRUE),
  data.frame(day = c(2, 3, 4),
             type = rep("noon-to-noon", 3),
             sleep_duration = c(7.5, 8, 8),
             sleep_midpoint = c(25.75%%24, 51%%24, 76%%24)))

  # midnight-to-midnight parsing
  m2m_df <- sleep24Summary(df=s_df, sleep_var = "S", time_var = "time", epoch_length = .1, noon_to_noon = FALSE)
  m2m_df$sleep_midpoint <- round(m2m_df$sleep_midpoint, 2) # round midpoint values
  expect_equal(m2m_df,
               data.frame(day = c(1, 2, 3),
                          type = rep("midnight-to-midnight", 3),
                          sleep_duration = c(9, 6.5, 7),
                          sleep_midpoint = c(2.61, 2.29, 3.5)))

})

test_that("sleep24Summary() handles floating point error w/ respect to epoch lengths and durations", {

  s_df <- data.frame("time"=seq(from = 0, to = 60, by = .5/60), "S" = 0) # 30-second epochs

  ## set sleep runs ##
  s_df[s_df$time >= 0 & s_df$time < 7, "S"] <- 1 # first run
  s_df[s_df$time >= 22 & s_df$time < 29.5, "S"] <- 1
  s_df[s_df$time >= 47 & s_df$time < 55, "S"] <- 1

  epoch_lengths <- unique(round(diff(s_df$time), 10)) # round to avoid floating point error

  # noon-to-noon parsing
  expect_equal(sleep24Summary(df=s_df, sleep_var = "S", time_var = "time", epoch_length = epoch_lengths, noon_to_noon = TRUE),
               data.frame(day = c(2, 3),
                          type = rep("noon-to-noon", 2),
                          sleep_duration = c(7.5, 8),
                          sleep_midpoint = c(25.75%%24, 51%%24)))

})



test_that("sleepSummary() correctly summarizes sleep runs", {

  ## simulated data ##
  s_df <- data.frame("time"=seq(from = 0, to = 168, by = .1),
                     "S" = 0)

  ## set sleep runs ##
  s_df[(s_df$time %% 24) >= 22 | (s_df$time %% 24) < 6, "S"] <- 1

  n2n_df1 <- sleep24Summary(s_df, "S", "time", .1, TRUE)
  m2m_df1 <- sleep24Summary(s_df, "S", "time", .1, FALSE)

  ## dropping any sleep run that hits the start or end of the data
  expect_equal(sleepSummary(df=s_df, sleep_var="S", time_var="time"),
               list(
                 summary = data.frame(
                   sleep_mid = timeMean(c(26, 50, 74, 98, 122, 146), c(8, 8, 8, 8, 8, 8)),
                   sleep_dur_noon_24hr = mean(n2n_df1$sleep_duration),
                   sleep_dur_midnight_24hr = mean(m2m_df1$sleep_duration)
                 ),
                 sleep_runs = data.frame(
                   sleep_onset = c(22, 46, 70, 94, 118, 142),
                   sleep_midpoint = c(26, 50, 74, 98, 122, 146),
                   sleep_offset = c(30, 54, 78, 102, 126, 150),
                   sleep_duration = c(8, 8, 8, 8, 8, 8),
                   start_index = c(221, 461, 701, 941, 1181, 1421),
                   end_index = c(300, 540, 780, 1020, 1260, 1500)
                 ),
                 noon_to_noon = n2n_df1,
                 midnight_to_midnight = m2m_df1
               )

  )

  ### Varied sleep windows ###
  ## simulated data ##
  s_df2 <- data.frame("time"=seq(from = 0, to = 60, by = .1),
                     "S" = 0)

  ## set sleep runs ##
  s_df2[s_df2$time >= 0 & s_df2$time < 6, "S"] <- 1 # sleep run # 1
  s_df2[s_df2$time >= 22.5 & s_df2$time < 30, "S"] <- 1 # sleep run # 2
  s_df2[s_df2$time >= 47 & s_df2$time < 58.3, "S"] <- 1 # sleep run # 3

  n2n_df2 <- sleep24Summary(s_df2, "S", "time", .1, TRUE)
  m2m_df2 <- sleep24Summary(s_df2, "S", "time", .1, FALSE)

  expect_equal(sleepSummary(df=s_df2, sleep_var="S", time_var="time"),
               list(
                 summary = data.frame(
                   sleep_mid = timeMean(c(26.25, 52.65), c(7.5, 11.3)),
                   sleep_dur_noon_24hr = mean(n2n_df2$sleep_duration),
                   sleep_dur_midnight_24hr = mean(m2m_df2$sleep_duration)
                 ),
                 sleep_runs = data.frame(
                   sleep_onset = c(22.5, 47),
                   sleep_midpoint = c(26.25, 52.65),
                   sleep_offset = c(30, 58.3),
                   sleep_duration = c(7.5, 11.3),
                   start_index = c(226, 471),
                   end_index = c(300, 583)
                 ),
                 noon_to_noon = n2n_df2,
                 midnight_to_midnight = m2m_df2
               )
  )

})

## TODO - Need to build in tests for error checking ##



