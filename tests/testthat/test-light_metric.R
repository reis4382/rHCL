# Tests of novel light metric ---------------------------------------------

test_that("completeDays() works", {
  ## generate fake data
  # start at 12 pm and have a 3.5 day duration (so first day gets dropped by function)
  # also, last row is technically a new day
  test_df <- data.frame(time = seq(12, 96, by = .1))

  res_df <- data.frame(
    day = c(1, 2, 3, 4, 5),
    start_time = c(12, 24, 48, 72, 96),
    end_time = c(23.9, 47.9, 71.9,  95.9, 96),
    duration = c(11.9, 23.9, 23.9, 23.9, 0)
  )

  expect_equal(res_df, completeDays(test_df, "time"))

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

      # ## generate test profile - initial state conditions are those at midnight of someone
      # ## with midsleep of 4:30 and sleep duration of 8.22 under typical Skeldon 2017 light profile
      #
      # NOTE: per personal correspondence with Prof. Skeldon, the actual parameters
      # used were s1=8 and s2=17 to account for DST
      #
      # # I'm going to try only altering mu and ca_par (circadian amplitude), as they are
      # # not in the Forger model and that would allow sleep duration and timing shifts without affecting tau
      #
      # times <- seq(0, 24*30, by = .5/60) # .5-minute intervals
      # light <- lightCycle(times, s1=8, s2=17) # generate light profile in skeldon 2017 paper (see function documentation for ref)
      #
      # desolve_list <- list(
      #   y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
      #   times = times,
      #   func = "derivsc_p",
      #   parms = unlist(hclParms(mu = 18.875, ca_par = 2.523)),
      #   dllname = "rHCL",
      #   initforc = "forcc_p",
      #   forcings = cbind(times, light),
      #   fcontrol = list(method = "linear", rule=2, f=0),
      #   initfunc = "parmsc_p",
      #   nout = 0,
      #   events = list(func = "eventc_p", root = TRUE),
      #   rootfun = "rootc_p",
      #   nroot = 1
      # )
      #
      # sol <- odeIter(desolve_args = desolve_list, max_iter = 40, dur_tol=1/60, mid_tol = 1/60)


      # ### generate test data ###
      # time <- seq(12, 96, by = .1) # start at 12 pm and have a 3.5 day duration (so first day gets dropped by function)
      # time <- time[1:(length(time)-1)] # drop final index, which would technically be a new day
      # light <- rep(0, length(time))
      # light[time>=30 & time <=44] <- 1000 # should speed up cycle on first complete day
      # # second complete day will stay at 0, which should slow down the cycle
      # light[time>=72 & time <=75] <-  1000 # should slow down cycle on day 3
      # light[time>=90] <- 1000 # should also slow down cycle on day 3
      #
      # test_df <- data.frame(time = time, light = light)
      #
      # res <- circLight(df = test_df, time_var = "time", light_var = "light")

    },
    finally = {

      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )


  browser()
  sol$sleep_sum


})

