test_that("new_rhcl_mod() makes a class", {

  res <- new_rhcl_mod()

  expect_equal(methods::is(res, "rhcl_mod"), TRUE)
  expect_equal(methods::is(unclass(res), "list"), TRUE)

})

test_that("summary.rhcl_mod method works", {

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
               opt_method = "optimize",
               duration_opt_control = durationOptControl(param_lower = 16.49, param_upper = 16.5,
                                                         tol = .1),
               midpoint_opt_control = midpointOptControl(param_lower = 24.49, param_upper = 24.5,
                                                         tol = .1)
  )

  res_sum <- summary(res)

  expect_equal(res_sum$opt$value, res$opt_results$value)
})
