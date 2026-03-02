# Tests of novel light metric ---------------------------------------------

test_that("completeDays() works", {
  ## generate fake data
  offset_vec <- c(rep(1, 12), rep(2, 24), rep(3, 24))

  res <- dayDurations(offset_vec = offset_vec, epoch_length_min = 60)

  expected_df <- data.frame(
    offset_day = c(1, 2, 3),
    day_duration = c(12, 24, 24)
  )

  expect_equal(res, expected_df)

})


test_that("circLight() works", {

  tryCatch(
    {
      ## check if .dll is loaded, load if needed (will unload after test)
      ## TODO - is there a better way of loading c code functions for testing?
      if(!"rHCL" %in% names(getLoadedDLLs())){
        # using here package to find root of rstudio project directory b/c
        # when running test suite the working directory switches to test folder
        dyn.load(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }

      ### generate test data ###
      time <- seq(12, 96, by = .1) # start at 12 pm and have a 3.5 day duration (so first day gets dropped by function)
      time <- time[1:(length(time)-1)] # drop final index, which would technically be a new day
      light <- rep(0, length(time))
      light[time>=30 & time <=44] <- 1000 # should speed up cycle on first complete day
      # second complete day will stay at 0, which should slow down the cycle
      light[time>=72 & time <=75] <-  1000 # should slow down cycle on day 3
      light[time>=90] <- 1000 # should also slow down cycle on day 3

      start_dtime <- as.POSIXct("2025-06-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
      dtime <- start_dtime + time * 60 * 60

      test_df <- data.frame(dtime = dtime, light = light)

      res <- circLight(df = test_df, time_var = "dtime", light_var = "light", epoch_length_min = 6)

    },
    finally = {

      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )

  ## expect speed up on day 1 and slowdown on days 2-3 ##
  # Not sure the actual minutes of change expected #
  expect_equal(sign(res[["df"]]$phase_diff), c(1, -1, -1))

  # check that the expected dates are being returned ##
  expect_equal(res[["df"]]$offset_day, c(2, 3, 4))
  expect_equal(res[["df"]]$offset_date, as.Date(c("2025-06-02", "2025-06-03", "2025-06-04"), tz = "America/Denver"))

  # check that there was more phase delay than advance, on average #
  expect_equal(res[["mean_diff"]] < 0, TRUE)

})

