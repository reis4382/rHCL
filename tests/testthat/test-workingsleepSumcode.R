test_that("sleepSum()works", {
  s_df1 <- data.frame("time"=seq(from = 12, to = 36, by = .1),
                      "S" = 0)
  s_df1$S[s_df1$time >= 22 & s_df1$time < 30.5] <- 1
  s_df2 <- data.frame("time"=seq(from = 12, to = 36, by = .2),
                      "S" = 0)
  s_df2$S[s_df2$time >= 24 & s_df2$time < 31.25] <- 1
  expect_equal(sleepSum(x=s_df1, sleep_var="S", time_var="time"), list(sleep_on = 22, sleep_mid = 26.25, sleep_off = 30.5, sleep_dur = 8.5))
  expect_equal(sleepSum(x=s_df2, sleep_var="S", time_var="time"), list(sleep_on = 24, sleep_mid = 27.7, sleep_off = 31.4, sleep_dur = 7.4))
})


