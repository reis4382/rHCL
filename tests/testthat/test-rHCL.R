
# Optimization control functions ------------------------------------------

test_that("durationOptControl() correctly formats arguments", {

  expect_equal(durationOptControl(maximum=TRUE)$maximum, TRUE)
  test <- durationOptControl()
  expect_equal(names(test), c("param_lower", "param_upper", "bisect_root_stop",
                              "bisect_max_iter", "bisect_abs_tol", "bisect_max_jumps"))
  expect_equal(test$param_upper, 30)

  expect_equal(durationOptControl(bisect_max_jumps = 20)$bisect_max_jumps, 20)

  expect_error(durationOptControl(param_upper = "Hello"),
               regexp = "The following argument\\(s\\) need to be numeric:")
  expect_error(durationOptControl(param_lower = "Hello"),
               regexp = "The param_lower argument in durationOptControl")

})

test_that("durationOptControl() rejects 'f' or 'interval' as arguments", {

  expect_error(durationOptControl(f = "Hello"),
               regexp = "durationOptControl\\(\\) cannot accept arguments for 'f'")
  expect_error(durationOptControl(interval = "Hello"),
               regexp = "durationOptControl\\(\\) cannot accept arguments for 'f'")
})

test_that("midpiontOptControl() correctly formats arguments", {
  expect_equal(midpointOptControl(maximum=TRUE)$maximum, TRUE)
  test <- midpointOptControl()
  expect_equal(names(test), c("param_lower", "param_upper", "bisect_root_stop",
                              "bisect_max_iter", "bisect_abs_tol", "bisect_max_jumps"))
  expect_equal(test$param_upper, 25)

  expect_equal(midpointOptControl(bisect_max_jumps = 20)$bisect_max_jumps, 20)

  expect_error(midpointOptControl(param_upper = "Hello"),
               regexp = "The following argument\\(s\\) need to be numeric:")

})

test_that("midpointOptControl() rejects 'f' or 'interval' as arguments", {

  expect_error(midpointOptControl(f = "Hello"),
               regexp = "midpointOptControl\\(\\) cannot accept arguments for 'f'")
  expect_error(midpointOptControl(interval = "Hello"),
               regexp = "midpointOptControl\\(\\) cannot accept arguments for 'f'")
})


# rhcl() tests ------------------------------------------------------------

test_that("rhcl() correctly optimizes parameters", {

  ## set up times and light entrainment profile
  start_time <- as.POSIXct("2025-06-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  end_time <- as.POSIXct("2025-06-07 23:59:00", format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")
  times <- seq(start_time, end_time, by = "1 min") # dtime

  light_times = ctimeCalc(times) # convert for use with lightCycle()

  df <- data.frame(
    times = times,
    light = lightCycle(light_times, l1=1000, l2=5) # generate standard light profile
  )


  ## set up a synthetic sleep wake cycle ##

  # create a list for deSolve::ode arguments #
  desolve_list <- list(
    # initial values, arbitrary
    y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
    times = light_times,
    func = "derivsc_p",
    parms = unlist(hclParms(mu = 16.5, tau = 24.5)), # set desired mu and tau values
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(light_times, df$light),
    fcontrol = list(method = "linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func = "eventc_p", root = TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  # extract sleep duration summary of synthetic data #
  syn_sol <- odeIter(desolve_args = desolve_list, dtime_vec = df$times,
                     max_iter = 20, dur_tol = 1/60, mid_tol = 1/60,
                     epoch_length_min = 1, min_observed_hours = 18)

  # check for convergence on synthetic data
  if(syn_sol$converge==FALSE){
    stop("Synthetic syn_sol data did not converge")
  }
  # extract results
  syn_mid <- syn_sol$sleep_sum$sleep_midpoint[nrow(syn_sol$sleep_sum)] # sleep midpoint
  syn_dur <- syn_sol$sleep_sum$sleep_duration[nrow(syn_sol$sleep_sum)] # sleep duration
  df$sleep <- syn_sol$ode_res$S # add sleep states to synthetic data

  ### Optimize both parameters ###
  min_mu <- desolve_list[["parms"]][["Hzero"]] + desolve_list[["parms"]][["ca_par"]] + desolve_list[["parms"]][["delta"]]*.05 # constrain mu to be greater than this
  if(min_mu > 16.5){
    stop("Mu for synthetic data is less than the allowed minimum")
  }


  ## Restrict range of optimization and reduce optimize/bisect tolerances to
  # speed up tests

  # start1 <- Sys.time()
  res <- rhcl(df = df,
              time_var = "times",
              sleep_var = "sleep",
              light_var = "light",
              epoch_length_min = 1,
              y0 = NULL,
              ode_parms = hclParms(),
              sleep_mid = syn_mid,
              sleep_dur = syn_dur,
              min_observed_hours = 18,
              max_ode_iter = 20,
              dur_tol = 1/60,
              mid_tol = 1/60,
              compiled = TRUE,
              opt_method = "bisect",
              duration_opt_control = durationOptControl(param_lower = 16.47, param_upper = 16.55,
                                                        bisect_root_stop = .1),
              midpoint_opt_control = midpointOptControl(param_lower = 24.47, param_upper = 24.55,
                                                        bisect_root_stop = .1)
  )
  # print(Sys.time()-start1)

  # start2 <- Sys.time()
  # reduce optimize tolerance to speed up test
  res2 <- rhcl(df = df,
               time_var = "times",
               sleep_var = "sleep",
               light_var = "light",
               epoch_length_min = 1,
               y0 = NULL,
               ode_parms = hclParms(),
               sleep_mid = syn_mid,
               sleep_dur = syn_dur,
               min_observed_hours = 18,
               max_ode_iter = 20,
               dur_tol = 1/60,
               mid_tol = 1/60,
               compiled = TRUE,
               opt_method = "optimize",
               duration_opt_control = durationOptControl(param_lower = 16.49, param_upper = 16.5,
                                                         tol = .1),
               midpoint_opt_control = midpointOptControl(param_lower = 24.49, param_upper = 24.5,
                                                         tol = .1)
  )
  # print(Sys.time() - start2)

  # R code is too slow
  # res2 <- rhcl(df = df,
  #              time_var = "times",
  #              sleep_var = "sleep",
  #              light_var = "light",
  #              y0 = NULL,
  #              ode_parms = hclParms(),
  #              sleep_mid = syn_mid,
  #              sleep_dur = syn_dur,
  #              max_iter = 20,
  #              dur_tol = 1/60,
  #              mid_tol = 1/60,
  #              compiled = FALSE
  # )

  ### Test that sleep information doesn't need to be provided ###
  res3 <- rhcl(df = df,
               time_var = "times",
               sleep_var = "sleep",
               light_var = "light",
               epoch_length_min = 1,
               y0 = NULL,
               ode_parms = hclParms(),
               sleep_mid = NULL,
               sleep_dur = NULL,
               min_observed_hours = 18,
               max_ode_iter = 20,
               dur_tol = 1/60,
               mid_tol = 1/60,
               compiled = TRUE,
               opt_method = "optimize",
               duration_opt_control = durationOptControl(param_lower = 16.49, param_upper = 16.5,
                                                         tol = .1),
               midpoint_opt_control = midpointOptControl(param_lower = 24.49, param_upper = 24.5,
                                                         tol = .1)
  )

  ## test bisection method ##
  expect_equal(abs(res$opt_results$value[1] - 16.5) < .05, TRUE)
  expect_equal(abs(res$opt_results$value[2] - 24.5) < .05, TRUE)
  expect_equal(res$opt_convergence_status, 1)
  expect_equal(res$ode_converge, TRUE)

  ## test optimize method ##
  expect_equal(abs(res2$opt_results$value[1] - 16.5) < .05, TRUE)
  expect_equal(abs(res2$opt_results$value[2] - 24.5) < .05, TRUE)
  expect_equal(res2$opt_convergence_status, 1)
  expect_equal(res2$ode_converge, TRUE)

  ## Test that sleep values are correctly used if summaries are not provided
  expect_equal(res2, res3)

})

