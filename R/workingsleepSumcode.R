#TEST DATAFRAME

## first simulated data ##
s_df1 <- data.frame("time"=seq(from = 12, to = 36, by = .1),
                    "S" = 0)

# convert to sleep starting at 10pm and switching to wake at 6:30am
s_df1$S[s_df1$time >= 22 & s_df1$time < 30.5] <- 1

## second simulated data ##
s_df2 <- data.frame("time"=seq(from = 12, to = 36, by = .2),
                    "S" = 0)

# convert to sleep times
s_df2$S[s_df2$time >= 24 & s_df2$time < 31.25] <- 1



##FUNCTION

sleepSum <- function(x) {
  x[order(x$time), ] #sort ascending

  sleep_first <- match("1", x$S) #find first instance of sleep
  sleep_on <- x[sleep_first, "time"] #time of sleep onset

  sleep_last <- max(which(x$S == 1)) + 1 #find last instance of sleep
  sleep_off <- x[sleep_last, "time"]  #awake onset

  sleep_dur <- sleep_off - sleep_on #duration of sleep in hours

  sleep_mid <- sleep_on + (sleep_dur/2) #midpoint of sleep

  return(list(sleep_on = sleep_on, sleep_mid = sleep_mid, sleep_off = sleep_off, sleep_dur = sleep_dur))

}


#TESTTHAT CODE, DATAFRAME 1
library(testthat)

##TEST S_DF1

test_that("sleepSum", {
  expect_equal(sleepSum(s_df1), list(sleep_on = 22, sleep_mid = 26.25, sleep_off = 30.5, sleep_dur = 8.5))
})


##NEW TEST S_DF2
test_that("sleepSum", {
  expect_equal(sleepSum(s_df2), list(sleep_on = 24, sleep_mid = 27.7, sleep_off = 31.4, sleep_dur = 7.4))
})
