
# Midpoint optimization tests ---------------------------------------------

test_that("odeOptim_midpoint() returns 13^2 on non-convergence", {

  ## set up times and light entrainment profile
  times <- seq(0, 24*30, by = .2) # 12-minute intervals
  light <- lightCycle(times) # default light profile
  light2 <- rep(0, length(times))

  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps

  # standard set up
  desolve_list <- list(
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(times, light),
    fcontrol = list(method="linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func="eventc_p", root=TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  # no light
  desolve_list2 <- desolve_list
  desolve_list2$forcings <- cbind(times, light2)

  # run through optim_midpoint function #
  res1 <- odeOptim_midpoint(tau_c = 27, sleep_mid = 5.25, desolve_args = desolve_list,
                            dtime_vec = dtimes, max_iter = 2, dur_tol = 1/60, mid_tol = 1/60,
                            epoch_length_min = 12, min_observed_hours = 18) # if testing excessively large tau that won't converge

  res2 <- odeOptim_midpoint(tau_c = 24.2, sleep_mid = 5.25, desolve_args = desolve_list2,
                            dtime_vec = dtimes, max_iter = 2, dur_tol = 1/60, mid_tol = 1/60,
                            epoch_length_min = 12, min_observed_hours = 18) # if model won't converge because of insufficient light
  expect_equal(res1, 13^2)
  expect_equal(res2, 13^2)

})



test_that("odeOptim_midpoint() works", {

  ## set up times and light entrainment profile
  times <- seq(0, 24*7, by = .2) # 12-minute intervals
  light <- lightCycle(times) # generate standard light profile

  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps

  ## set up a synthetic sleep wake cycle for tauc = 23.99##

  # create a list for deSolve::ode arguments #
  desolve_list <- list(
    # initial values, should correspond to a 4:24 am ymin time and ~7.4 hours of nightly sleep
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms(tau_c = 23.99)), # set desired tau value
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

  # extract sleep midpoint summary of synthetic data #
  syn_sol <- odeIter(desolve_args = desolve_list,
                     dtime_vec = dtimes, max_iter = 20, dur_tol = 1/60, mid_tol = 1/60,
                     epoch_length_min = 12, min_observed_hours = 18)
  syn_mid <- syn_sol$sleep_sum$sleep_midpoint[nrow(syn_sol$sleep_sum)]

  # optimize - optimize if for 1D optimization.
  # lowering tolerance to speed up test
  opt_midpoint <- optimize(f = odeOptim_midpoint, interval = c(23.8, 24.1), sleep_mid = syn_mid,
                           desolve_args = desolve_list, dtime_vec = dtimes, max_iter = 20,
                           dur_tol = 1/60, mid_tol = 1/60,
                           epoch_length_min = 12, min_observed_hours = 18, tol = .01)


  # ## alternative tau_c value - tests take too long, so limiting to one value ##
  # desolve_list2 <- desolve_list
  # desolve_list2[["forcings"]] <- cbind(times, lightCycle(times, l2=10)) # decrease night light exposure to aid more extreme entrainment
  # desolve_list2$parms[["tau_c"]] <- 24.4 # alternative tau_c
  #
  #
  # # generate synthetic results and extract sleep midpoint
  # syn_sol2 <- odeIter(desolve_args = desolve_list2,
  #                     dtime_vec = dtimes, max_iter = 20, dur_tol = 1/60, mid_tol = 1/60,
  #                     epoch_length_min = 12, min_observed_hours = 18)
  # # check for convergence on synthetic data
  # if(syn_sol2$converge==FALSE){
  #   stop("Synthetic syn_sol2 data did not converge")
  # }
  # syn_mid2 <- syn_sol2$sleep_sum$sleep_midpoint[nrow(syn_sol2$sleep_sum)]
  #
  # # optimize
  # opt_midpoint2 <- optimize(f = odeOptim_midpoint, interval = c(24.3, 24.5), sleep_mid = syn_mid2,
  #                          desolve_args = desolve_list2, dtime_vec = dtimes,
  #                          max_iter = 20, dur_tol = 1/60, mid_tol = 1/60,
  #                          epoch_length_min = 12, min_observed_hours = 18)


  ## check if optimizer got close to values when rounded ##
  expect_equal(round(opt_midpoint$minimum, 2), 23.99)
  # expect_equal(round(opt_midpoint2$minimum, 2), 24.40)

  # alternative checks allowing for close values that may round differently
  # expect_equal(abs(opt_midpoint$minimum - 23.99) <= .01, TRUE) # test that estimated minimum for first data is close to 23.99
  # expect_equal(abs(opt_midpoint2$minimum - 24.4) <= .01, TRUE) # test that estimated minimum for second data is close to 24.4

  ## check that deviances (already squared) for minimum values are below .03 ##
  expect_equal(opt_midpoint$objective < .03, TRUE)
  # expect_equal(opt_midpoint2$objective < .03, TRUE)

})


# Duration optimization tests ---------------------------------------------

test_that("odeOptim_duration() returns 24^2 upon non-convergence", {

  ## set up times and light entrainment profile
  times <- seq(0, 24*30, by = .2) # 12-minute intervals
  light <- lightCycle(times) # default light profile
  light2 <- rep(0, length(times))

  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps


  # standard set up
  desolve_list <- list(
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(times, light),
    fcontrol = list(method="linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func="eventc_p", root=TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  # no light
  desolve_list2 <- desolve_list
  desolve_list2$forcings <- cbind(times, light2)

  # run through optim_midpoint function #
  res1 <- odeOptim_duration(mu = 1000, sleep_dur = 7.4, desolve_args = desolve_list,
                            dtime_vec = dtimes, max_iter = 2, dur_tol = 1/60, mid_tol = 1/60,
                            epoch_length_min = 12, min_observed_hours = 18) # if testing excessively large tau that won't converge
  res2 <- odeOptim_duration(mu = 17.87, sleep_dur = 7.4, desolve_args = desolve_list2,
                            dtime_vec = dtimes, max_iter = 2, dur_tol = 1/60, mid_tol = 1/60,
                            epoch_length_min = 12, min_observed_hours = 18) # if model won't converge because of insufficient light


  expect_equal(res1, 24^2)
  expect_equal(res2, 24^2)
})



test_that("odeOptim_duration() works", {

  ## set up times and light entrainment profile
  times <- seq(0, 24*7, by = .1) # 12-minute intervals
  light <- lightCycle(times) # generate standard light profile

  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps

  ## set up a synthetic sleep wake cycle ##

  # create a list for deSolve::ode arguments #
  desolve_list <- list(
    # initial values, should correspond to a 4:24 am ymin time and ~7.4 hours of nightly sleep
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms(mu = 16.5)), # set desired mu value
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

  # extract sleep duration summary of synthetic data #
  syn_sol <- odeIter(desolve_args = desolve_list,
                     dtime_vec = dtimes, max_iter = 20, dur_tol = 1/60, mid_tol = 1/60,
                     epoch_length_min = 6, min_observed_hours = 18)
  syn_dur <- syn_sol$sleep_sum$sleep_duration[nrow(syn_sol$sleep_sum)]

  # optimize - optimize if for 1D optimization.
  min_mu <- desolve_list[["parms"]][["Hzero"]] + desolve_list[["parms"]][["ca_par"]] + desolve_list[["parms"]][["delta"]]*.05 # constrain mu to be greater than this

  # lowering tolerance to speed up convergence
  opt_duration <- optimize(f = odeOptim_duration, interval = c(16.3, 16.7), sleep_dur = syn_dur,
                           desolve_args = desolve_list, dtime_vec = dtimes,
                           max_iter = 20, dur_tol = 1/60, mid_tol = 1/60,
                           epoch_length_min = 6, min_observed_hours = 18, tol = .1)


  # ## alternative mu value - tests take too long, so limiting to one value##
  # desolve_list2 <- desolve_list
  # desolve_list2[["forcings"]] <- cbind(times, lightCycle(times, l2=10)) # decrease night light exposure to aid more extreme entrainment
  # desolve_list2$parms[["mu"]] <- 21.04 # alternative mu
  #
  #
  # # generate synthetic results and extract sleep midpoint
  # syn_sol2 <- odeIter(desolve_args = desolve_list2,
  #                     dtime_vec = dtimes, max_iter = 20, dur_tol = 1/60, mid_tol = 1/60,
  #                     epoch_length_min = 6, min_observed_hours = 18)
  # # check for convergence on synthetic data
  # if(syn_sol2$converge==FALSE){
  #   stop("Synthetic syn_sol2 data did not converge")
  # }
  # syn_dur2 <- syn_sol2$sleep_sum$sleep_duration[nrow(syn_sol2$sleep_sum)]
  #
  # # optimize
  # min_mu2 <- desolve_list2[["parms"]][["Hzero"]] + desolve_list2[["parms"]][["ca_par"]] + desolve_list2[["parms"]][["delta"]]*.05 # constrain mu to be greater than this
  #
  # opt_duration2 <- optimize(f = odeOptim_duration, interval = c(min_mu2 , 30), sleep_dur = syn_dur2,
  #                           desolve_args = desolve_list2, dtime_vec = dtimes,
  #                           max_iter = 20, dur_tol = 1/60, mid_tol = 1/60,
  #                           epoch_length_min = 6, min_observed_hours = 18)


  # ## check if optimizer fairly close to actual values ##
  expect_equal(abs(opt_duration$minimum - 16.5) < .15, TRUE)
  # expect_equal(abs(opt_duration2$minimum - 21.04) < .15, TRUE)

  ## check that deviances (already squared) for minimum values are below .03 ##
  expect_equal(opt_duration$objective < .03, TRUE)
  # expect_equal(opt_duration2$objective < .03, TRUE)

})


# Joint parameter optimization --------------------------------------------

# test_that("odeOptim_tauAndMu() returns 13^2 + 24^2 upon non-convergence", {
#       ## set up times and light entrainment profile
#       times <- seq(0, 24*30, by = .2) # 12-minute intervals
#       light <- lightCycle(times) # default light profile
#       light2 <- rep(0, length(times))
#
#       # standard set up
#       desolve_list <- list(
#         y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
#         times = times,
#         func = "derivsc_p",
#         parms = unlist(hclParms()),
#         dllname = "rHCL",
#         initforc = "forcc_p",
#         forcings = cbind(times, light),
#         fcontrol = list(method="linear", rule=2, f=0),
#         initfunc = "parmsc_p",
#         nout = 0,
#         events = list(func="eventc_p", root=TRUE),
#         rootfun = "rootc_p",
#         nroot = 1
#       )
#
#       # no light
#       desolve_list2 <- desolve_list
#       desolve_list2$forcings <- cbind(times, light2)
#
#       # run through optim_tauAndMu function #
#       res1 <- odeOptim_tauAndMu(pars=c(24.2, 1000), sleep_mid = 3.5, sleep_dur = 7.4, desolve_args = desolve_list,
#                                 max_iter = 20, dur_tol = 1/60, mid_tol = 1/60) # if testing excessively large mu that won't converge
#       res2 <- odeOptim_tauAndMu(pars=c(35, 17), sleep_mid = 3.5, sleep_dur = 7.4, desolve_args = desolve_list,
#                                 max_iter = 20, dur_tol = 1/60, mid_tol = 1/60) # if testing excessively large tau that won't converge
#       res3 <- odeOptim_tauAndMu(pars=c(24.2, 17.87), sleep_dur = 7.4, desolve_args = desolve_list2,
#                                 max_iter = 20, dur_tol = 1/60, mid_tol = 1/60) # if model won't converge because of insufficient light
#
#
#   expect_equal(res1, 13^2 + 24^2)
#   expect_equal(res2, 13^2 + 24^2)
#   expect_equal(res3, 13^2 + 24^2)
# })
#
#
#
# test_that("odeOptim_tauAndMu() works", {
#
#       ## set up times and light entrainment profile
#       times <- seq(0, 24*7, by = .2) # 12-minute intervals
#       light <- lightCycle(times, l2=10) # generate standard light profile
#
#       ## set up a synthetic sleep wake cycle ##
#
#       # create a list for deSolve::ode arguments #
#       desolve_list <- list(
#         # initial values, should correspond to a 4:24 am ymin time and ~7.4 hours of nightly sleep
#         y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
#         times = times,
#         func = "derivsc_p",
#         parms = unlist(hclParms(tau_c = 24.4, mu = 16.5)), # set desired mu value
#         dllname = "rHCL",
#         initforc = "forcc_p",
#         forcings = cbind(times, light),
#         fcontrol = list(method = "linear", rule=2, f=0),
#         initfunc = "parmsc_p",
#         nout = 0,
#         events = list(func = "eventc_p", root = TRUE),
#         rootfun = "rootc_p",
#         nroot = 1
#       )
#
#       # extract sleep outcomes of synthetic data #
#       syn_sol <- odeIter(desolve_args = desolve_list,
#                          max_iter = 20, dur_tol = 1/60, mid_tol = 1/60)
#       if(syn_sol$converge==FALSE){
#         stop("Synthetic syn_sol data did not converge")
#       }
#
#       syn_mid <- syn_sol$sleep_sum$sleep_midpoint[nrow(syn_sol$sleep_sum)]
#       syn_dur <- syn_sol$sleep_sum$sleep_duration[nrow(syn_sol$sleep_sum)]
#
#       # optimize for multiple parameters
#       opt_optim1 <- optim(par = c(24.2, 17), fn = odeOptim_tauAndMu,
#                           sleep_mid = syn_mid, sleep_dur = syn_dur,
#                           desolve_args = desolve_list, max_iter = 20,
#                           dur_tol = 1/60, mid_tol = 1/60,
#                           control = list(maxit = 1000))
#
#       # # try bobyqa for gradient free optimization
#       # min_mu <- desolve_list[["parms"]][["Hzero"]] + desolve_list[["parms"]][["ca_par"]] + desolve_list[["parms"]][["delta"]]*.05 # constrain mu to be greater than this
#       #
#       # opt_bobyqa1 <- nloptr::bobyqa(x0 = c(24.2, 17.87), fn = odeOptim_tauAndMu,
#       #                               lower = c(22, min_mu),
#       #                               upper = c(26, 30),
#       #                               sleep_mid = syn_mid, sleep_dur = syn_dur,
#       #                               desolve_args = desolve_list, max_iter = 20,
#       #                               dur_tol = 1/60, mid_tol = 1/60,
#       #                               control = list(xtol_rel = 1e-8))
#
#
#       ## alternative mu value ##
#       desolve_list2 <- desolve_list
#       desolve_list2[["forcings"]] <- cbind(times, lightCycle(times, l2=10)) # decrease night light exposure to aid more extreme entrainment
#       desolve_list2$parms[["tau_c"]] <- 23.8
#       desolve_list2$parms[["mu"]] <- 19.9 # alternative mu
#
#
#       # generate synthetic results and extract sleep midpoint
#       syn_sol2 <- odeIter(desolve_args = desolve_list2,
#                           max_iter = 20, dur_tol = 1/60, mid_tol = 1/60)
#       # check for convergence on synthetic data
#       if(syn_sol2$converge==FALSE){
#         stop("Synthetic syn_sol2 data did not converge")
#       }
#       syn_mid2 <- syn_sol2$sleep_sum$sleep_midpoint[nrow(syn_sol2$sleep_sum)]
#       syn_dur2 <- syn_sol2$sleep_sum$sleep_duration[nrow(syn_sol2$sleep_sum)]
#
#       # optimize
#       opt_optim2 <- optim(par = c(24.2, 18.5), fn = odeOptim_tauAndMu,
#                           sleep_mid = syn_mid2, sleep_dur = syn_dur2,
#                           desolve_args = desolve_list2, max_iter = 20,
#                           dur_tol = 1/60, mid_tol = 1/60)
# #
# #       # try bobyqa for gradient free optimization
# #       min_mu2 <- desolve_list2[["parms"]][["Hzero"]] + desolve_list2[["parms"]][["ca_par"]] + desolve_list2[["parms"]][["delta"]]*.05 # constrain mu to be greater than this
# #
# #       opt_bobyqa2 <- nloptr::bobyqa(x0 = c(24.2, 18.5), fn = odeOptim_tauAndMu,
# #                                     lower = c(22, min_mu2),
# #                                     upper = c(26, 30),
# #                                     sleep_mid = syn_mid2, sleep_dur = syn_dur2,
# #                                     desolve_args = desolve_list2, max_iter = 20,
# #                                     dur_tol = 1/60, mid_tol = 1/60,
# #                                     control = list(xtol_rel = 1e-8))
#
#   # ## check if optimizer close to actual values - Not the most accurate ##
#   expect_equal(abs(opt_optim1$par[1] - 24.4) < .15, TRUE)
#   expect_equal(abs(opt_optim1$par[2] - 16.5) < .2, TRUE)
#
#   expect_equal(abs(opt_optim2$par[1] - 23.8) < .15, TRUE)
#   expect_equal(abs(opt_optim2$par[2] - 19.9) < .15, TRUE)
#
#   ## check that deviances (squared) for minimum values are below .03 ##
#   expect_equal(opt_optim1$value < .03, TRUE)
#   expect_equal(opt_optim2$value < .03, TRUE)
# })


# residualCheck function --------------------------------------------------

test_that("residualCheck() works", {
  expect_equal(residualCheck(1, .05, square = TRUE), FALSE)
  expect_equal(residualCheck(-1, .05, square = TRUE), FALSE)
  expect_equal(residualCheck(.0001, .05, square = TRUE), TRUE)
  expect_equal(residualCheck(-.03, -.03, square = TRUE), TRUE)
  expect_equal(residualCheck(.015, .014, square = FALSE), TRUE)
  expect_equal(residualCheck(.1, .014, square = FALSE), FALSE)
})


# bisectWhileLoop tests ---------------------------------------------------

test_that("bisectWhileLoop() correctly adjusts non-convergence of ODEs", {
  ## set up times and light entrainment profile
  times = seq(0, 24*7, by = .2) # 12-minute intervals
  light1 <- lightCycle(times) # generate standard light profile
  light2 <- rep(0, length(times)) # generate no light

  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps

  # create a list for deSolve::ode arguments #
  desolve_list <- list(
    # initial values, arbitrary
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(times, light1),
    fcontrol = list(method = "linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func = "eventc_p", root = TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  desolve_list2 <- list(
    # initial values, arbitrary
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(times, light2),
    fcontrol = list(method = "linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func = "eventc_p", root = TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  # lowering tolerance to speed up convergence
  res1 <- bisectWhileLoop(24.2, "tau_c", lower_bound = 24, upper_bound = 26,
                          max_steps = 6,
                          desolve_args = desolve_list, dtime_vec = dtimes,
                          max_ode_iter = 20,
                          dur_tol = 1/60, mid_tol = 1/60,
                          epoch_length_min = 12, min_observed_hours = 18)

  res2 <- bisectWhileLoop(24.2, "tau_c", lower_bound = 24.1, upper_bound = 26,
                          max_steps = 6,
                          desolve_args = desolve_list2, dtime_vec = dtimes,
                          max_ode_iter = 2,
                          dur_tol = 1/60, mid_tol = 1/60,
                          epoch_length_min = 12, min_observed_hours = 18)

  expect_equal(res1$param_val, 24.2)
  expect_equal(res2$param_val, NA)

})


# odeBisect() tests -------------------------------------------------------
test_that("bisectBisect() works", {

  ## set up times and light entrainment profile
  times = seq(0, 24*7, by = .05) # 3-minute intervals
  light <- lightCycle(times, l1=1000, l2=5) # generate standard light profile

  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps

  # create a list for deSolve::ode arguments #
  desolve_list <- list(
    # initial values, arbitrary
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms(mu = 16.5, tau = 24.5)), # set desired mu and tau values
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

  # extract sleep duration summary of synthetic data #
  syn_sol <- odeIter(desolve_args = desolve_list, dtime_vec = dtimes,
                     max_iter = 20, dur_tol = 1/60, mid_tol = 1/60,
                     epoch_length_min = 3, min_observed_hours = 18)

  # check for convergence on synthetic data
  if(syn_sol$converge==FALSE){
    stop("Synthetic syn_sol data did not converge")
  }
  # extract results
  syn_mid <- syn_sol$sleep_sum$sleep_midpoint[nrow(syn_sol$sleep_sum)] # sleep midpoint
  syn_dur <- syn_sol$sleep_sum$sleep_duration[nrow(syn_sol$sleep_sum)] # sleep duration

  ## estimate mu ##
  min_mu <- desolve_list[["parms"]][["Hzero"]] + desolve_list[["parms"]][["ca_par"]] + desolve_list[["parms"]][["delta"]]*.05 # constrain mu to be greater than this
  desolve_list[["parms"]][["tau_c"]] <- 24.2 # for estimating mu, tau_c can be fixed to 24 to aid convergence

  # lowering tolerance (root_stop) to speed up convergence
  res_duration <- odeBisect(
    param_lower = 16.3,
    param_upper = 16.7,
    observed_param = syn_dur,
    root_stop = .01,
    max_iter = 100,
    abs_tol = 1e-8,
    method = "mu",
    num_ode_jumps = 10,
    desolve_args = desolve_list,
    dtime_vec = dtimes,
    max_ode_iter = 20,
    dur_tol = 1/60,
    mid_tol = 1/60,
    epoch_length_min = 3,
    min_observed_hours = 18
  )


  ## estimate tau ##
  desolve_list[["parms"]][["mu"]] <- res_duration$minimum # estimated value of mu
  # desolve_list[["parms"]][["mu"]] <- 17.87 # default value of mu

  res_midpoint <- odeBisect(
    param_lower = 24.3,
    param_upper = 24.7,
    observed_param = syn_mid,
    root_stop = .01,
    max_iter = 100,
    abs_tol = 1e-8,
    method = "tau_c",
    num_ode_jumps = 20,
    desolve_args = desolve_list,
    dtime_vec = dtimes,
    max_ode_iter = 20,
    dur_tol = 1/60,
    mid_tol = 1/60,
    epoch_length_min = 3,
    min_observed_hours = 18
  )


  expect_equal(abs(res_duration$minimum - 16.5) < .05, TRUE)
  expect_equal(res_duration$objective < .03, TRUE)
  expect_equal(abs(res_midpoint$minimum - 24.5) < .05, TRUE)
  expect_equal(res_midpoint$objective < .03, TRUE)

})

test_that("odeBisect() correctly returns errors", {


  ## set up times and light entrainment profile
  times = seq(0, 24*7, by = .05) # 3-minute intervals
  light <- lightCycle(times, l1=1000, l2=5) # generate standard light profile

  start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
  dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps

  # create a list for deSolve::ode arguments #
  desolve_list <- list(
    # initial values, arbitrary
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms()), # set desired mu and tau values
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

  ## Both boundaries overestimate sleep ##
  expect_error(odeBisect(
    param_lower = 20,
    param_upper = 30,
    observed_param = 7.5,
    root_stop = 1e-4,
    max_iter = 100,
    abs_tol = 1e-8,
    method = "mu",
    num_ode_jumps = 10,
    desolve_args = desolve_list,
    dtime_vec = dtimes,
    max_ode_iter = 20,
    dur_tol = 1/60,
    mid_tol = 1/60,
    epoch_length_min = 3,
    min_observed_hours = 18

  ), regexp = "Boundaries for mu.*same sign \\(positive\\).*")

  ## Both boundaries underestimate sleep ##
  expect_error(odeBisect(
    param_lower = 15,
    param_upper = 16,
    observed_param = 7.5,
    root_stop = 1e-4,
    max_iter = 100,
    abs_tol = 1e-8,
    method = "mu",
    num_ode_jumps = 10,
    desolve_args = desolve_list,
    dtime_vec = dtimes,
    max_ode_iter = 20,
    dur_tol = 1/60,
    mid_tol = 1/60,
    epoch_length_min = 3,
    min_observed_hours = 18
  ), regexp = "Boundaries for mu.*same sign \\(negative\\).*")

  ## Upper bound is not greater than lower bound
  expect_error(odeBisect(
    param_lower = 17.87,
    param_upper = 16.5,
    observed_param = 7.5,
    root_stop = 1e-4,
    max_iter = 100,
    abs_tol = 1e-8,
    method = "mu",
    num_ode_jumps = 10,
    desolve_args = desolve_list,
    dtime_vec = dtimes,
    max_ode_iter = 20,
    dur_tol = 1/60,
    mid_tol = 1/60,
    epoch_length_min = 3,
    min_observed_hours = 18
  ), regexp = "Upper bound \\(param_upper\\) for mu must be greater than.*")

  ## Lower bound won't converge
  expect_error(odeBisect(
    param_lower = 20,
    param_upper = 22,
    observed_param = 2.5,
    root_stop = 1e-4,
    max_iter = 100,
    abs_tol = 1e-8,
    method = "tau_c",
    num_ode_jumps = 1,
    desolve_args = desolve_list,
    dtime_vec = dtimes,
    max_ode_iter = 20,
    dur_tol = 1/60,
    mid_tol = 1/60,
    epoch_length_min = 3,
    min_observed_hours = 18
  ), regexp = "No value for tau_c at or near the lower boundary.*")

  # getting the lower bound to converge, but not the upper bound, doesn't
  # seem possible, because the sequencing approach, as is, will also try
  # setting the upper bound to the lower bound at the end. Another error
  # was created to catch situations where the two bounds become identical.

  ## Upper bound becomes lower bound
  expect_error(odeBisect(
    param_lower = 24,
    param_upper = 35,
    observed_param = 2.5,
    root_stop = 1e-4,
    max_iter = 100,
    abs_tol = 1e-8,
    method = "tau_c",
    num_ode_jumps = 1,
    desolve_args = desolve_list,
    dtime_vec = dtimes,
    max_ode_iter = 20,
    dur_tol = 1/60,
    mid_tol = 1/60,
    epoch_length_min = 3,
    min_observed_hours = 18
  ), regexp = ".*Try increasing num_ode_jumps.")

})

# test_that("Test for odeBisect() speed up if passing updated starting values", {
#
#   ## set up times and light entrainment profile
#   times = seq(0, 24*7, by = .05) # 3-minute intervals
#   light <- lightCycle(times, l1=1000, l2=5) # generate standard light profile
#
#   start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
#   dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps
#
#   # create a list for deSolve::ode arguments #
#   desolve_list <- list(
#     # initial values, arbitrary
#     y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
#     times = times,
#     func = "derivsc_p",
#     parms = unlist(hclParms(mu = 16.5, tau = 24.5)), # set desired mu and tau values
#     dllname = "rHCL",
#     initforc = "forcc_p",
#     forcings = cbind(times, light),
#     fcontrol = list(method = "linear", rule=2, f=0),
#     initfunc = "parmsc_p",
#     nout = 0,
#     events = list(func = "eventc_p", root = TRUE),
#     rootfun = "rootc_p",
#     nroot = 1
#   )
#
#   # extract sleep duration summary of synthetic data #
#   syn_sol <- odeIter(desolve_args = desolve_list, dtime_vec = dtimes,
#                      max_iter = 20, dur_tol = 1/60, mid_tol = 1/60,
#                      epoch_length_min = 3, min_observed_hours = 18)
#
#   # check for convergence on synthetic data
#   if(syn_sol$converge==FALSE){
#     stop("Synthetic syn_sol data did not converge")
#   }
#   # extract results
#   syn_mid <- syn_sol$sleep_sum$sleep_midpoint[nrow(syn_sol$sleep_sum)] # sleep midpoint
#   syn_dur <- syn_sol$sleep_sum$sleep_duration[nrow(syn_sol$sleep_sum)] # sleep duration
#
#   ## estimate mu ##
#   min_mu <- desolve_list[["parms"]][["Hzero"]] + desolve_list[["parms"]][["ca_par"]] + desolve_list[["parms"]][["delta"]]*.05 # constrain mu to be greater than this
#   desolve_list[["parms"]][["tau_c"]] <- 24.2 # for estimating mu, tau_c can be fixed to 24 to aid convergence
#
#   # not updating starting state #
#   test1_start <- Sys.time()
#   test1 <- odeBisect(
#     param_lower = 13,
#     param_upper = 30,
#     observed_param = syn_dur,
#     root_stop = 1e4,
#     max_iter = 100,
#     abs_tol = 1e-8,
#     method = "mu",
#     num_ode_jumps = 10,
#     desolve_args = desolve_list,
#     dtime_vec = dtimes,
#     max_ode_iter = 20,
#     dur_tol = 1/60,
#     mid_tol = 1/60,
#     epoch_length_min = 3,
#     min_observed_hours = 18,
#     update_y0 = FALSE
#   )
#   test1_end <- Sys.time()
#
#   # updating starting state #
#   test2_start <- Sys.time()
#   test2 <- odeBisect(
#     param_lower = 13,
#     param_upper = 30,
#     observed_param = syn_dur,
#     root_stop = 1e4,
#     max_iter = 100,
#     abs_tol = 1e-8,
#     method = "mu",
#     num_ode_jumps = 10,
#     desolve_args = desolve_list,
#     dtime_vec = dtimes,
#     max_ode_iter = 20,
#     dur_tol = 1/60,
#     mid_tol = 1/60,
#     epoch_length_min = 3,
#     min_observed_hours = 18,
#     update_y0 = TRUE
#   )
#   test2_end <- Sys.time()
#
#
#   (test2_end - test2_start) < (test1_end - test1_start)
#
#   expect_equal(abs(test2$minimum - 16.5) < .05, TRUE)
#   expect_equal(test2$objective < .03, TRUE)
#
#
#
#
# })
