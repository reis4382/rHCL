
# Midpoint optimization tests ---------------------------------------------

test_that("odeOptim_midpoint() returns 13^2 on non-convergence", {

  ## using tryCatch() so that dyn.load is always removed,
  ## even if the test errors out
  tryCatch(
    {
      ## check if .dll is loaded, load if needed (will unload after test)
      ## TODO - is there a better way of loading c code functions for testing?
      if(!"rHCL" %in% names(getLoadedDLLs())){
        # using here package to find root of rstudio project directory b/c
        # when running test suite the working directory switches to test folder
        dyn.load(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }

      ## set up times and light entrainment profile
      times <- seq(0, 24*30, by = .2) # 12-minute intervals
      light <- lightCycle(times) # default light profile
      light2 <- rep(0, length(times))

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
                                max_iter = 20, dur_tol = 1/60, mid_tol = 1/60) # if testing excessively large tau that won't converge
      res2 <- odeOptim_midpoint(tau_c = 24.2, sleep_mid = 5.25, desolve_args = desolve_list2,
                                max_iter = 20, dur_tol = 1/60, mid_tol = 1/60) # if model won't converge because of insufficient light


    },
    finally = {
      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )

  expect_equal(res1, 13^2)
  expect_equal(res2, 13^2)

})



test_that("odeOptim_midpoint() works", {


  ## using tryCatch() so that dyn.load is always removed,
  ## even if the test errors out
  tryCatch(
    {
      ## check if .dll is loaded, load if needed (will unload after test)
      ## TODO - is there a better way of loading c code functions for testing?
      if(!"rHCL" %in% names(getLoadedDLLs())){
        # using here package to find root of rstudio project directory b/c
        # when running test suite the working directory switches to test folder
        dyn.load(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }

      ## set up times and light entrainment profile
      times <- seq(0, 24*7, by = .2) # 12-minute intervals
      light <- lightCycle(times) # generate standard light profile

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
                         max_iter = 20, dur_tol = 1/60, mid_tol = 1/60)
      syn_mid <- syn_sol$sleep_sum$sleep_midpoint[nrow(syn_sol$sleep_sum)]

      # optimize - optimize if for 1D optimization.
      opt_midpoint <- optimize(f = odeOptim_midpoint, interval = c(22, 26), sleep_mid = syn_mid,
                      desolve_args = desolve_list, max_iter = 20,
                      dur_tol = 1/60, mid_tol = 1/60)


      ## alternative tau_c value ##
      desolve_list2 <- desolve_list
      desolve_list2[["forcings"]] <- cbind(times, lightCycle(times, l2=10)) # decrease night light exposure to aid more extreme entrainment
      desolve_list2$parms[["tau_c"]] <- 24.4 # alternative tau_c


      # generate synthetic results and extract sleep midpoint
      syn_sol2 <- odeIter(desolve_args = desolve_list2,
                          max_iter = 20, dur_tol = 1/60, mid_tol = 1/60)
      # check for convergence on synthetic data
      if(syn_sol2$converge==FALSE){
        stop("Synthetic syn_sol2 data did not converge")
      }
      syn_mid2 <- syn_sol2$sleep_sum$sleep_midpoint[nrow(syn_sol2$sleep_sum)]

      # optimize
      opt_midpoint2 <- optimize(f = odeOptim_midpoint, interval = c(22, 26), sleep_mid = syn_mid2,
                               desolve_args = desolve_list2, max_iter = 20,
                               dur_tol = 1/60, mid_tol = 1/60)


    },
    finally = {

      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )

  ## check if optimizer got close to values when rounded ##
  expect_equal(round(opt_midpoint$minimum, 2), 23.99)
  expect_equal(round(opt_midpoint2$minimum, 2), 24.40)

  # alternative checks allowing for close values that may round differently
  # expect_equal(abs(opt_midpoint$minimum - 23.99) <= .01, TRUE) # test that estimated minimum for first data is close to 23.99
  # expect_equal(abs(opt_midpoint2$minimum - 24.4) <= .01, TRUE) # test that estimated minimum for second data is close to 24.4

  ## check that deviances (squared) for minimum values are below .03 ##
  expect_equal(opt_midpoint$objective < .03, TRUE)
  expect_equal(opt_midpoint2$objective < .03, TRUE)

})


# Duration optimization tests ---------------------------------------------

test_that("odeOptim_duration() returns 24^2 upon non-convergence", {


  ## using tryCatch() so that dyn.load is always removed,
  ## even if the test errors out
  tryCatch(
    {
      ## check if .dll is loaded, load if needed (will unload after test)
      ## TODO - is there a better way of loading c code functions for testing?
      if(!"rHCL" %in% names(getLoadedDLLs())){
        # using here package to find root of rstudio project directory b/c
        # when running test suite the working directory switches to test folder
        dyn.load(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }

      ## set up times and light entrainment profile
      times <- seq(0, 24*30, by = .2) # 12-minute intervals
      light <- lightCycle(times) # default light profile
      light2 <- rep(0, length(times))

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
                                max_iter = 20, dur_tol = 1/60, mid_tol = 1/60) # if testing excessively large tau that won't converge
      res2 <- odeOptim_duration(mu = 17.87, sleep_dur = 7.4, desolve_args = desolve_list2,
                                max_iter = 20, dur_tol = 1/60, mid_tol = 1/60) # if model won't converge because of insufficient light


    },
    finally = {
      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )

  expect_equal(res1, 24^2)
  expect_equal(res2, 24^2)
})



test_that("odeOptim_duration() works", {


  ## using tryCatch() so that dyn.load is always removed,
  ## even if the test errors out
  tryCatch(
    {
      ## check if .dll is loaded, load if needed (will unload after test)
      ## TODO - is there a better way of loading c code functions for testing?
      if(!"rHCL" %in% names(getLoadedDLLs())){
        # using here package to find root of rstudio project directory b/c
        # when running test suite the working directory switches to test folder
        dyn.load(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }

      ## set up times and light entrainment profile
      times <- seq(0, 24*7, by = .2) # 12-minute intervals
      light <- lightCycle(times) # generate standard light profile

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
                         max_iter = 20, dur_tol = 1/60, mid_tol = 1/60)
      syn_dur <- syn_sol$sleep_sum$sleep_duration[nrow(syn_sol$sleep_sum)]

      # optimize - optimize if for 1D optimization.
      min_mu <- desolve_list[["parms"]][["Hzero"]] + desolve_list[["parms"]][["ca_par"]] + desolve_list[["parms"]][["delta"]]*.05 # constrain mu to be greater than this

      opt_duration <- optimize(f = odeOptim_duration, interval = c(min_mu, 30), sleep_dur = syn_dur,
                               desolve_args = desolve_list, max_iter = 20,
                               dur_tol = 1/60, mid_tol = 1/60)


      ## alternative mu value ##
      desolve_list2 <- desolve_list
      desolve_list2[["forcings"]] <- cbind(times, lightCycle(times, l2=10)) # decrease night light exposure to aid more extreme entrainment
      desolve_list2$parms[["mu"]] <- 21.04 # alternative mu


      # generate synthetic results and extract sleep midpoint
      syn_sol2 <- odeIter(desolve_args = desolve_list2,
                          max_iter = 20, dur_tol = 1/60, mid_tol = 1/60)
      # check for convergence on synthetic data
      if(syn_sol2$converge==FALSE){
        stop("Synthetic syn_sol2 data did not converge")
      }
      syn_dur2 <- syn_sol2$sleep_sum$sleep_duration[nrow(syn_sol2$sleep_sum)]

      # optimize
      min_mu2 <- desolve_list2[["parms"]][["Hzero"]] + desolve_list2[["parms"]][["ca_par"]] + desolve_list2[["parms"]][["delta"]]*.05 # constrain mu to be greater than this

      opt_duration2 <- optimize(f = odeOptim_duration, interval = c(min_mu2 , 30), sleep_dur = syn_dur2,
                                desolve_args = desolve_list2, max_iter = 20,
                                dur_tol = 1/60, mid_tol = 1/60)


    },
    finally = {

      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )

  # ## check if optimizer fairly close to actual values ##
  expect_equal(abs(opt_duration$minimum - 16.5) < .15, TRUE)
  expect_equal(abs(opt_duration2$minimum - 21.04) < .15, TRUE)

  ## check that deviances (squared) for minimum values are below .03 ##
  expect_equal(opt_duration$objective^2 < .03, TRUE)
  expect_equal(opt_duration2$objective^2 < .03, TRUE)

})


# Joint parameter optimization --------------------------------------------

test_that("odeOptim_tauAndMu() returns 13^2 + 24^2 upon non-convergence", {


  ## using tryCatch() so that dyn.load is always removed,
  ## even if the test errors out
  tryCatch(
    {
      ## check if .dll is loaded, load if needed (will unload after test)
      ## TODO - is there a better way of loading c code functions for testing?
      if(!"rHCL" %in% names(getLoadedDLLs())){
        # using here package to find root of rstudio project directory b/c
        # when running test suite the working directory switches to test folder
        dyn.load(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }

      ## set up times and light entrainment profile
      times <- seq(0, 24*30, by = .2) # 12-minute intervals
      light <- lightCycle(times) # default light profile
      light2 <- rep(0, length(times))

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

      # run through optim_tauAndMu function #
      res1 <- odeOptim_tauAndMu(pars=c(24.2, 1000), sleep_mid = 3.5, sleep_dur = 7.4, desolve_args = desolve_list,
                                max_iter = 20, dur_tol = 1/60, mid_tol = 1/60) # if testing excessively large mu that won't converge
      res2 <- odeOptim_tauAndMu(pars=c(35, 17), sleep_mid = 3.5, sleep_dur = 7.4, desolve_args = desolve_list,
                                max_iter = 20, dur_tol = 1/60, mid_tol = 1/60) # if testing excessively large tau that won't converge
      res3 <- odeOptim_tauAndMu(pars=c(24.2, 17.87), sleep_dur = 7.4, desolve_args = desolve_list2,
                                max_iter = 20, dur_tol = 1/60, mid_tol = 1/60) # if model won't converge because of insufficient light


    },
    finally = {
      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )

  expect_equal(res1, 13^2 + 24^2)
  expect_equal(res2, 13^2 + 24^2)
  expect_equal(res3, 13^2 + 24^2)
})



test_that("odeOptim_tauAndMu() works", {


  ## using tryCatch() so that dyn.load is always removed,
  ## even if the test errors out
  tryCatch(
    {
      ## check if .dll is loaded, load if needed (will unload after test)
      ## TODO - is there a better way of loading c code functions for testing?
      if(!"rHCL" %in% names(getLoadedDLLs())){
        # using here package to find root of rstudio project directory b/c
        # when running test suite the working directory switches to test folder
        dyn.load(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }

      ## set up times and light entrainment profile
      times <- seq(0, 24*7, by = .2) # 12-minute intervals
      light <- lightCycle(times, l2=10) # generate standard light profile

      ## set up a synthetic sleep wake cycle ##

      # create a list for deSolve::ode arguments #
      desolve_list <- list(
        # initial values, should correspond to a 4:24 am ymin time and ~7.4 hours of nightly sleep
        y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
        times = times,
        func = "derivsc_p",
        parms = unlist(hclParms(tau_c = 24.4, mu = 16.5)), # set desired mu value
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

      # extract sleep outcomes of synthetic data #
      syn_sol <- odeIter(desolve_args = desolve_list,
                         max_iter = 20, dur_tol = 1/60, mid_tol = 1/60)
      if(syn_sol$converge==FALSE){
        stop("Synthetic syn_sol data did not converge")
      }

      syn_mid <- syn_sol$sleep_sum$sleep_midpoint[nrow(syn_sol$sleep_sum)]
      syn_dur <- syn_sol$sleep_sum$sleep_duration[nrow(syn_sol$sleep_sum)]

      # optimize for multiple parameters
      opt_optim1 <- optim(par = c(24.2, 17), fn = odeOptim_tauAndMu,
                          sleep_mid = syn_mid, sleep_dur = syn_dur,
                          desolve_args = desolve_list, max_iter = 20,
                          dur_tol = 1/60, mid_tol = 1/60,
                          control = list(maxit = 1000))

      # # try bobyqa for gradient free optimization
      # min_mu <- desolve_list[["parms"]][["Hzero"]] + desolve_list[["parms"]][["ca_par"]] + desolve_list[["parms"]][["delta"]]*.05 # constrain mu to be greater than this
      #
      # opt_bobyqa1 <- nloptr::bobyqa(x0 = c(24.2, 17.87), fn = odeOptim_tauAndMu,
      #                               lower = c(22, min_mu),
      #                               upper = c(26, 30),
      #                               sleep_mid = syn_mid, sleep_dur = syn_dur,
      #                               desolve_args = desolve_list, max_iter = 20,
      #                               dur_tol = 1/60, mid_tol = 1/60,
      #                               control = list(xtol_rel = 1e-8))


      ## alternative mu value ##
      desolve_list2 <- desolve_list
      desolve_list2[["forcings"]] <- cbind(times, lightCycle(times, l2=10)) # decrease night light exposure to aid more extreme entrainment
      desolve_list2$parms[["tau_c"]] <- 23.8
      desolve_list2$parms[["mu"]] <- 19.9 # alternative mu


      # generate synthetic results and extract sleep midpoint
      syn_sol2 <- odeIter(desolve_args = desolve_list2,
                          max_iter = 20, dur_tol = 1/60, mid_tol = 1/60)
      # check for convergence on synthetic data
      if(syn_sol2$converge==FALSE){
        stop("Synthetic syn_sol2 data did not converge")
      }
      syn_mid2 <- syn_sol2$sleep_sum$sleep_midpoint[nrow(syn_sol2$sleep_sum)]
      syn_dur2 <- syn_sol2$sleep_sum$sleep_duration[nrow(syn_sol2$sleep_sum)]

      # optimize
      opt_optim2 <- optim(par = c(24.2, 18.5), fn = odeOptim_tauAndMu,
                          sleep_mid = syn_mid2, sleep_dur = syn_dur2,
                          desolve_args = desolve_list2, max_iter = 20,
                          dur_tol = 1/60, mid_tol = 1/60)
#
#       # try bobyqa for gradient free optimization
#       min_mu2 <- desolve_list2[["parms"]][["Hzero"]] + desolve_list2[["parms"]][["ca_par"]] + desolve_list2[["parms"]][["delta"]]*.05 # constrain mu to be greater than this
#
#       opt_bobyqa2 <- nloptr::bobyqa(x0 = c(24.2, 18.5), fn = odeOptim_tauAndMu,
#                                     lower = c(22, min_mu2),
#                                     upper = c(26, 30),
#                                     sleep_mid = syn_mid2, sleep_dur = syn_dur2,
#                                     desolve_args = desolve_list2, max_iter = 20,
#                                     dur_tol = 1/60, mid_tol = 1/60,
#                                     control = list(xtol_rel = 1e-8))

    },
    finally = {

      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )

  # ## check if optimizer close to actual values - Not the most accurate ##
  expect_equal(abs(opt_optim1$par[1] - 24.4) < .15, TRUE)
  expect_equal(abs(opt_optim1$par[2] - 16.5) < .2, TRUE)

  expect_equal(abs(opt_optim2$par[1] - 23.8) < .15, TRUE)
  expect_equal(abs(opt_optim2$par[2] - 19.9) < .15, TRUE)

  ## check that deviances (squared) for minimum values are below .03 ##
  expect_equal(opt_optim1$value < .03, TRUE)
  expect_equal(opt_optim2$value < .03, TRUE)
})


# residualCheck function --------------------------------------------------

test_that("residualCheck() works", {
  expect_equal(residualCheck(1, .05), FALSE)
  expect_equal(residualCheck(-1, .05), FALSE)
  expect_equal(residualCheck(.0001, .05), TRUE)
  expect_equal(residualCheck(-.03, -.03), TRUE)
})
