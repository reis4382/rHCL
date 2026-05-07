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
test_that("odeIterPrep() correctly replicates data", {

  df <- data.frame(
    ctime = 0:27,
    light = 0:27*100
  )

  res1 <- odeIterPrep(df, ctime_var = "ctime", light_var = "light", max_ode_iter = 3, tol = 1/60/60)

  expected_df <- list(
    df = data.frame(
    ctime = c(0:(27+48)),
    light = c(rep(0:23*100, 2), 0:27*100)
  ),
  orig_length = nrow(df),
  full_days = 1,
  final_ind = 24
  )
  expect_equal(res1, expected_df)

  res2 <- odeIterPrep(df, ctime_var = "ctime", light_var = "light", max_ode_iter = 1, tol = 1/60/60)
  expect_equal(res2, list(df = df, orig_length = nrow(df), full_days = 0, final_ind = 28))

})

test_that("ode_iter() correctly iterates over data until convergence", {

  ## set up times and light entrainment profile
  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")

  times <- seq(0, 24*30, by = .05) # 3-minute intervals, ctime format
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps

  light <- rep(0, length(times)) # light vector
  light[(times %% 24) > 8 & (times %%24) < 22] <- 1000 # 1000 lux exposure from 8 am - 10 pm

  ode_df <- data.frame(
    times = times,
    light = light
  )

  ode_df <- odeIterPrep(ode_df, ctime_var = "times", light_var = "light", max_ode_iter = 5, tol = 1/60/60)

  # create a list for deSolve::ode arguments #
  desolve_list <- list(
    y = c(h = 13, n = 0, x = 0, y = 1, S = 0),
    times = ode_df[["df"]]$times,
    func = "derivsc_p",
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(ode_df[["df"]]$times, ode_df[["df"]]$light),
    fcontrol = list(method = "linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func = "eventc_p", root = TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  sol <- odeIter(desolve_args = desolve_list, dtime_vec = dtimes, max_ode_iter = 5,
                 orig_length = ode_df[["orig_length"]], full_days = ode_df[["full_days"]],
                 final_ind = ode_df[["final_ind"]],
                 dur_tol = 3/60, mid_tol = 3/60, epoch_length_min = 12, min_observed_hours = 18)

  ## No light exposure ##
  desolve_list2 <- desolve_list
  desolve_list2[["forcings"]] <- cbind(ode_df[["df"]]$times, rep(0, nrow(ode_df[["df"]])))

  sol2 <- odeIter(desolve_args = desolve_list2,  dtime_vec = dtimes, max_ode_iter = 5,
                  orig_length = ode_df[["orig_length"]], full_days = ode_df[["full_days"]],
                  final_ind = ode_df[["final_ind"]],
                  dur_tol = 3/60, mid_tol = 3/60, epoch_length_min = 12, min_observed_hours = 18)

  ## Check - model converges when sufficient light ##
  expect_equal(sol$converge, TRUE)
  expect_match(sol$conv_message, regexp = "Convergence obtained with 5 iterations.")
  # sleep summaries should be identical for these fake data
  expect_equal(sol$sleep_sum$sleep_midpoint[1], sol$sleep_sum$sleep_midpoint[1])
  expect_equal(sol$sleep_sum$sleep_duration[1], sol$sleep_sum$sleep_duration[1])

  ## Check - model does not converge when no light ##
  expect_equal(sol2$converge, FALSE)
  expect_match(sol2$conv_message, regexp = "The model did not converge.")
  expect_equal(sol2$iterations, 5)
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

  ode_df <- odeIterPrep(df = data.frame(times = times, light = light), ctime_var = "times",
                        light_var = "light", max_ode_iter = 20, tol = 1/60/60)

  # create a list for deSolve::ode arguments - C code#
  desolve_list <- list(
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = ode_df[["df"]]$times,
    func = "derivsc_p",
    # parms = unlist(hclParms()),
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(ode_df[["df"]]$times, ode_df[["df"]]$light),
    fcontrol = list(method = "linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func = "eventc_p", root = TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  sol <- odeIter(desolve_args = desolve_list, dtime_vec = dtimes, orig_length = ode_df[["orig_length"]],
                 full_days = ode_df[["full_days"]], final_ind = ode_df[["final_ind"]],
                 max_ode_iter = 20, mid_tol = 1/60, dur_tol = 1/60,
                 epoch_length_min = 1, min_observed_hours = 18)

  # # attempt with some manual light gating #
  # desolve_list2 <- desolve_list
  # desolve_list2[["forcings"]] <- cbind(times, light2)
  # sol2 <- odeIter(desolve_args = desolve_list2, max_ode_iter = 40)


  expect_equal(round(sol$sleep_sum$sleep_midpoint[nrow(sol$sleep_sum)], 2), 3.43) # Midpoint ~ 3:26 am (after rounding)
  # expect_equal(round(sol2$sleep_sum$sleep_midpoint[nrow(sol2$sleep_sum)], 2), 3.55)

})

test_that("odeIter() returns the same final results for different starting values", {

  ## Skeldon 2017 paper default light values ##
  times <- seq(0, 24*30, by = .2) # 12-minute intervals
  light <- lightCycle(times, l1=700, l2=40) # generate light profile in skeldon 2017 paper (see function documentation for ref)

  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps

  ode_df <- odeIterPrep(df = data.frame(times = times, light = light), ctime_var = "times",
                        light_var = "light", max_ode_iter = 20, tol = 1/60/60)

  # create a list for deSolve::ode arguments - C code#
  desolve_list <- list(
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = ode_df[["df"]]$times,
    func = "derivsc_p",
    # parms = unlist(hclParms()),
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(ode_df[["df"]]$times, ode_df[["df"]]$light),
    fcontrol = list(method = "linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func = "eventc_p", root = TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  sol <- odeIter(desolve_args = desolve_list, dtime_vec = dtimes, orig_length = ode_df[["orig_length"]],
                 full_days = ode_df[["full_days"]], final_ind = ode_df[["final_ind"]],
                 max_ode_iter = 20, mid_tol = 1/60, dur_tol = 1/60,
                 epoch_length_min = 1, min_observed_hours = 18)

  ## alternative starting values ##
  desolve_list2 <- desolve_list
  desolve_list2[["y"]] <- c(h = 15, n = .3, x = -1, y = -0, S = 0)
  sol2 <- odeIter(desolve_args = desolve_list2, dtime_vec = dtimes, orig_length = ode_df[["orig_length"]],
                 full_days = ode_df[["full_days"]], final_ind = ode_df[["final_ind"]],
                 max_ode_iter = 20, mid_tol = 1/60, dur_tol = 1/60,
                 epoch_length_min = 1, min_observed_hours = 18)

  # compare final sleep summary results, as specific ODE values may have slight differences
  # also, may have different number of iterations, so reset row name
  res1 <- sol$sleep_sum[nrow(sol$sleep_sum), c("sleep_midpoint", "sleep_duration")] # extract relevant columns
  row.names(res1) <- 1:nrow(res1)

  res2 <- sol2$sleep_sum[nrow(sol2$sleep_sum), c("sleep_midpoint", "sleep_duration")] # extract relevant columns
  row.names(res2) <- 1:nrow(res2)

  expect_equal(res1, res2)

})

