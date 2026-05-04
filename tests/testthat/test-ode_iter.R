test_that("convergeCheck() works", {

  # convergence check passes
  expect_equal(convergeCheck(sleep_dur1=8.3, sleep_dur2=8.2,
                             sleep_mid1=3.1, sleep_mid2=3.15,
                             dur_tol = .1, mid_tol = .1),
               list(deviations = data.frame(variable = c("duration", "midpoint"),
                                            devs = c(.1, .05),
                                            converge = c(TRUE, TRUE)),
                    overall_converge = TRUE))

  # convergence check fails duration check
  expect_equal(convergeCheck(sleep_dur1=7.5, sleep_dur2=8.2,
                             sleep_mid1=3.1, sleep_mid2=3.15,
                             dur_tol = .1, mid_tol = .1),
               list(deviations = data.frame(variable = c("duration", "midpoint"),
                                            devs = c(0.7, .05),
                                            converge = c(FALSE, TRUE)),
                    overall_converge = FALSE))

  # convergence check fails midpoint check
  expect_equal(convergeCheck(sleep_dur1=8.15, sleep_dur2=8.2,
                             sleep_mid1=23.5, sleep_mid2=3.15,
                             dur_tol = .1, mid_tol = .1),
               list(deviations = data.frame(variable = c("duration", "midpoint"),
                                            devs = c(0.05, 3.65),
                                            converge = c(TRUE, FALSE)),
                    overall_converge = FALSE))

  # convergence check fails both checks
  expect_equal(convergeCheck(sleep_dur1=6.15, sleep_dur2=10.5,
                             sleep_mid1=1.2, sleep_mid2=4,
                             dur_tol = .1, mid_tol = .1),
               list(deviations = data.frame(variable = c("duration", "midpoint"),
                                            devs = c(4.35, 2.8),
                                            converge = c(FALSE, FALSE)),
                    overall_converge = FALSE))


})

## TODO - test that convergeCheck catches NAs


# Test initial state checking for sleep and h -----------------------------

test_that("initialStateCheck() correctly adjusts invalid combinations of initial state variables", {

  expect_equal(initialStateCheck(x=.5, y=-.5, S=0, h=18, hzero=13, ca_par=1.72, delta=1), 1)
  expect_equal(initialStateCheck(x=.5, y=-.5, S=0, h=13, hzero=13, ca_par=1.72, delta=1), 0)
  expect_equal(initialStateCheck(x=.5, y=-.5, S=1, h=11, hzero=13, ca_par=1.72, delta=1), 0)
  expect_equal(initialStateCheck(x=.5, y=-.5, S=1, h=14, hzero=13, ca_par=1.72, delta=1), 1)

})


# Test ode iteration functions --------------------------------------------

test_that("ode_iter() correctly iterates over data until convergence", {

  ## set up times and light entrainment profile
  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")

  times <- seq(0, 24*30, by = .2) # 12-minute intervals, ctime format
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps

  light <- rep(0, length(times)) # light vector
  light[(times %% 24) > 8 & (times %%24) < 22] <- 1000 # 1000 lux exposure from 8 am - 10 pm

  # create a list for deSolve::ode arguments #
  desolve_list <- list(
    y = c(h = 13, n = 0, x = 0, y = 1, S = 0),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(times, light),
    fcontrol = list(method = "linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func = "eventc_p", root = TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  sol <- odeIter(desolve_args = desolve_list, dtime_vec = dtimes, max_iter = 20,
                 dur_tol = 1/60, mid_tol = 1/60, epoch_length_min = 12, min_observed_hours = 18)

  ## No light exposure ##
  desolve_list2 <- desolve_list
  desolve_list2[["forcings"]] <- cbind(times, rep(0, length(times)))

  sol2 <- odeIter(desolve_args = desolve_list2,  dtime_vec = dtimes, max_iter = 21,
                  dur_tol = 1/60, mid_tol = 1/60, epoch_length_min = 12, min_observed_hours = 18)

  ## Check - model converges when sufficient light ##
  expect_equal(sol$converge, TRUE)
  expect_match(sol$conv_message, regexp = "Convergence obtained after 3 iterations.")

  ## Check - model does not converge when no light ##
  expect_equal(sol2$converge, FALSE)
  expect_match(sol2$conv_message, regexp = "The model did not converge.")
  expect_equal(sol2$iterations, 21)
})

test_that("odeIter() returns correct sleep midpoint", {
  ## Skeldon 2017 paper default light values ##
  ## presumably, default light profile w/ default model parameters should generate
  ## a midsleep time ~3:16 am (per skeldon 2023 paper text). CORRECTION:
  ## Per personal correspondence w/ Prof. Skeldon, midsleep timing should be 3:25 am
  ## when using s1 = 8, s2 = 17 or 2:55 am when using s1 = 7.5, s2 = 16.5
  times <- seq(0, 24*30, by = 1/60) # 1-minute intervals

  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps

  # generate light profile in skeldon 2017 paper (see function documentation for ref). Note,
  # per personal correspondence w/ Prof. Skeldon, s1 should equal 8 and s2 should equal 17
  # (to account for "half" a year of DST).
  light <- lightCycle(times*60*60, l1=700, l2=40, c=.6, s1 = 8, s2 = 17, time_scale = "secs")
  # light <- lightCycle(times, l1=700, l2=40, s1 = 7.5, s2 = 16.5, time_scale = "secs") # I may have c parameter wrong here

  # # Light was noted to be gated as in Figure 1, which appears to be between midnight and ~8 am
  # light2 <- light
  # # light2[(times%%24) < 8] <- 0 # gated between midnight and 8 am
  # light2[(times%%24) < 7.75] <- 0 # gated between midnight and 7:45 am

  # create a list for deSolve::ode arguments - C code#
  desolve_list <- list(
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = times,
    func = "derivsc_p",
    # parms = unlist(hclParms()),
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(times, light),
    fcontrol = list(method = "linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func = "eventc_p", root = TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  sol <- odeIter(desolve_args = desolve_list, dtime_vec = dtimes,
                 max_iter = 40, mid_tol = 1/60, dur_tol = 1/60,
                 epoch_length_min = 1, min_observed_hours = 18)

  # # attempt with some manual light gating #
  # desolve_list2 <- desolve_list
  # desolve_list2[["forcings"]] <- cbind(times, light2)
  # sol2 <- odeIter(desolve_args = desolve_list2, max_iter = 40)


  expect_equal(round(sol$sleep_sum$sleep_midpoint[nrow(sol$sleep_sum)], 2), 3.43) # Midpoint ~ 3:26 am (after rounding)
  # expect_equal(round(sol2$sleep_sum$sleep_midpoint[nrow(sol2$sleep_sum)], 2), 3.55)

})

test_that("odeIter() returns the same final results for different starting values", {

  ## Skeldon 2017 paper default light values ##
  ## presumably, default light profile w/ default model parameters should generate
  ## a midsleep time ~3:16 am (per skeldon 2023 paper text).
  times <- seq(0, 24*30, by = .2) # 12-minute intervals
  light <- lightCycle(times, l1=700, l2=40) # generate light profile in skeldon 2017 paper (see function documentation for ref)

  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps

  # create a list for deSolve::ode arguments - C code#
  desolve_list <- list(
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(times, light),
    fcontrol = list(method = "linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func = "eventc_p", root = TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  sol <- odeIter(desolve_args = desolve_list, dtime_vec = dtimes,
                 max_iter = 20, mid_tol = 1/60, dur_tol = 1/60,
                 epoch_length_min = 12, min_observed_hours = 18)

  ## alternative starting values ##
  desolve_list2 <- desolve_list
  desolve_list2[["y"]] <- c(h = 15, n = .3, x = -1, y = -0, S = 0)
  sol2 <- odeIter(desolve_args = desolve_list2, dtime_vec = dtimes,
                  max_iter = 20, mid_tol = 1/60, dur_tol = 1/60,
                  epoch_length_min = 12, min_observed_hours = 18)

  # compare final sleep summary results, as specific ODE values may have slight differences
  # also, may have different number of iterations, so reset row name
  res1 <- sol$sleep_sum[nrow(sol$sleep_sum), c("sleep_midpoint", "sleep_duration")] # extract relevant columns
  row.names(res1) <- 1:nrow(res1)

  res2 <- sol2$sleep_sum[nrow(sol2$sleep_sum), c("sleep_midpoint", "sleep_duration")] # extract relevant columns
  row.names(res2) <- 1:nrow(res2)

  expect_equal(res1, res2)

})
