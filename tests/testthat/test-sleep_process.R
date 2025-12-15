test_that("sleepSummary() correctly summarizes sleep runs", {

  ## simulated data ##
  s_df <- data.frame("time"=seq(from = 0, to = 168, by = .1),
                     "S" = 0)

  ## set sleep runs ##
  s_df[(s_df$time %% 24) >= 22 | (s_df$time %% 24) < 6, "S"] <- 1

  ## Same test, but dropping any sleep run that hits the start or end of the data
  expect_equal(sleepSummary(x=s_df, sleep_var="S", time_var="time"),
               data.frame(
                 sleep_onset = c(22, 46, 70, 94, 118, 142),
                 sleep_midpoint = c(26, 50, 74, 98, 122, 146),
                 sleep_offset = c(30, 54, 78, 102, 126, 150),
                 sleep_duration = c(8, 8, 8, 8, 8, 8),
                 start_index = c(221, 461, 701, 941, 1181, 1421),
                 end_index = c(300, 540, 780, 1020, 1260, 1500)
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

  expect_equal(sleepSummary(x=s_df2, sleep_var="S", time_var="time"),
               data.frame(
                 sleep_onset = c(22.5, 47),
                 sleep_midpoint = c(26.25, 52.65),
                 sleep_offset = c(30, 58.3),
                 sleep_duration = c(7.5, 11.3),
                 start_index = c(226, 471),
                 end_index = c(300, 583)
               )
  )

})
