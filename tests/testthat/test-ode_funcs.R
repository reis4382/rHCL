# Tests of hclParms(), which generates ODE parameters -----------------------

test_that("hclParms() returns correct format and values", {
  expect_equal(class(hclParms()), "list") # test that a list is returned
  expect_equal(length(hclParms()), 16) # test length of returned list
  expect_equal(names(hclParms()), c("mu", "chi", "Hzero", "delta", "ca_par",
                                    "tau_c", "f_par", "G_par", "p_par", "k_par",
                                    "little_b_par", "gamma", "alpha_zero",
                                    "beta", "Izero", "kappa")) # correct names are returned
  # manually check all terms
  expect_equal(as.numeric(hclParms()),
               c(17.87, 45, 13, 1, 1.72, 24.2, 0.99669, 19.9, 0.6, 0.55, 0.4, 0.23, 0.16, 0.013, 9500, 12/pi))
})

test_that("hclParms() correctly accepts argument changes", {
  expect_equal(hclParms(tau_c = 24)[["tau_c"]], 24)
  expect_equal(as.numeric(hclParms(mu = 15, gamma = .5)[c("mu", "gamma")]), c(15, 0.5))
})

test_that("hclParms() correctly accepts time_scale arguments", {
  expect_equal(as.numeric(hclParms(time_scale="hours")[c("chi", "kappa")]), c(45, 12/pi)) # scale to hours
  expect_equal(as.numeric(hclParms(time_scale="mins")[c("chi", "kappa")]), c(45*60, 12/pi*60)) # scale to minutes
  expect_equal(as.numeric(hclParms(time_scale="secs")[c("chi", "kappa")]), c(45*60*60, 12/pi*60*60)) # scale to seconds

  expect_equal(as.numeric(hclParms(time_scale=1)[c("chi", "kappa")]), c(45, 12/pi)) # scale to hours
  expect_equal(as.numeric(hclParms(time_scale=1/60)[c("chi", "kappa")]), c(45*60, 12/pi*60)) # scale to minutes
  expect_equal(as.numeric(hclParms(time_scale=1/60/60)[c("chi", "kappa")]), c(45*60*60, 12/pi*60*60)) # scale to seconds

  expect_error(hclParms(time_scale = "days"), regexp = "time_scale argument must be either a numeric.*")
})


test_that("hclParms() correctly rejects incorrect parameter names", {
  expect_error(hclParms(whoops = 5), regexp = "unused argument \\(whoops = 5\\).*")
  expect_error(hclParms(arg = 3, whoops = 5), regexp = "unused arguments \\(arg = 3.*")
  expect_error(hclParms(k_bar = 5), regexp = "unused argument \\(k_bar = 5\\).*")
})

test_that("hclParms() rejects non-numeric arguments (other than time_scale)", {
  expect_error(hclParms(mu = "nope", p_par = "hi"), regexp = "The following arguments need to be numeric: mu, p_par")

})


# Compare C and R ODE code ------------------------------------------------

test_that("R and C code for derivatives return same results", {
  ## using tryCatch() so that light.int() and dyn.load are always removed,
  ## even if the test errors out
  tryCatch(
    {
      ## create times
      times = c(0, 12, 24, 24.2, 48, 48.4, 72, 72.6, 96, 96.8, 120, 121, 144, 145.2)

      ## Create light interpolation function for R code ##
      assign("light.int",
             approxfun(x=times, y=rep(0, length(times)), method="linear", rule=2),
             envir = .GlobalEnv) # create interpolation function; Note this is creating a global environment and needs to be cleaned up
      ## TODO - Is there a way to make light.int() available for testing but not place it in the global env?

      ## check if .dll is loaded, load if needed (will unload after test)
      ## TODO - is there a better way of loading c code functions for testing?
      if(!"rHCL" %in% names(getLoadedDLLs())){
        # using here package to find root of rstudio project directory b/c
        # when running test suite the working directory switches to test folder
        dyn.load(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }

      # first test - NO LIGHT
      sol_r <- deSolve::ode(y = c(h = 13, n = 0, x = 1, y = 0, S = 0),
                            func = dHCL,
                            times = times,
                            parms = hclParms(),
                            events = list(func = dEventFunc, root = TRUE),
                            rootfun = dRootFunc)

      sol_c <- deSolve::ode(
        y = c(h = 13, n = 0, x = 1, y = 0, S = 0),
        times = times,
        func = "derivsc_p",
        parms = unlist(hclParms()),
        dllname = "rHCL",
        initforc = "forcc_p",
        forcings = cbind(times, rep(0, length(times))),
        fcontrol = list(method="linear", rule=2, f=0),
        initfunc = "parmsc_p",
        nout = 0,
        events = list(func="eventc_p", root=TRUE),
        rootfun = "rootc_p",
        nroot = 1
      )

      ### second test - Light exposure (no roots) ###
      times <- seq(from=0, to=24.2, by = .1) # start times
      light <- rep(0, length(times)) # start light vector
      light[times > 8 & times < 22] <- 5000 # light exposure during "day"
      # re-assign light func
      assign("light.int",
             approxfun(x=times, y=light, method="linear", rule=2),
             envir = .GlobalEnv) # create interpolation function; Note this is creating a global environment and needs to be cleaned up

      sol_r2 <- deSolve::ode(y = c(h = 13, n = 0, x = 1, y = 0, S = 0),
                             func = dHCL,
                             times = times,
                             parms = hclParms())

      sol_c2 <- deSolve::ode(
        y = c(h = 13, n = 0, x = 1, y = 0, S = 0),
        times = times,
        func = "derivsc_p",
        parms = unlist(hclParms()),
        dllname = "rHCL",
        initforc = "forcc_p",
        forcings = cbind(times, light),
        fcontrol = list(method="linear", rule=2, f=0),
        initfunc = "parmsc_p",
        nout = 0,
        # events = list(func="eventc_p", root=TRUE),
        # rootfun = "rootc_p",
        # nroot = 1
      )

      ### Third test - Light exposure w/ roots ###
      sol_r3 <- deSolve::ode(y = c(h = 13, n = 0, x = 1, y = 0, S = 0),
                             func = dHCL,
                             times = times,
                             parms = hclParms(),
                             events = list(func = dEventFunc, root = TRUE),
                             rootfun = dRootFunc)

      sol_c3 <- deSolve::ode(
        y = c(h = 13, n = 0, x = 1, y = 0, S = 0),
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

    },
    finally = {
      ## Clean up the function I added to the global environment ##
      if(exists("light.int", where = .GlobalEnv)){
        rm(light.int, envir = .GlobalEnv)
      } # remove light approxfun

      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )

  expect_equal(as.numeric(sol_r), as.numeric(sol_c))
  expect_equal(as.numeric(sol_r2), as.numeric(sol_c2))
  expect_equal(as.numeric(sol_r3), as.numeric(sol_c3))

})


# Test that ODEs return expected values -----------------------------------

test_that("ODE output in the absence of light", {

  ## using tryCatch() so that light.int() and dyn.load are always removed,
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

      ## Run up to 30 days of lightlessness (no roots)
      # tau 24.2
      ## create times
      # tau_test <- c(24.18, 24.19, 24.2, 24.21, 24.22)
      # times <- c(0, tau_test, tau_test*4, tau_test*30)
      times <- seq(0, 91, by = .2)

      sol_tau_24pt2 <- deSolve::ode(
        y = c(h = 13, n = 0, x = 0, y = -1, S = 0),
        times = times,
        func = "derivsc_p",
        parms = unlist(hclParms(tau_c=18)),
        dllname = "rHCL",
        initforc = "forcc_p",
        forcings = cbind(times, rep(0, length(times))),
        fcontrol = list(method="linear", rule=2, f=0),
        initfunc = "parmsc_p",
        nout = 0,
        # events = list(func="eventc_p", root=TRUE),
        # rootfun = "rootc_p",
        # nroot = 1
      )

      # # tau 18
      # ## create times
      # tau_test <- c(17.98, 17.99, 18, 18.01, 18.02)
      # times <- c(0, tau_test, tau_test*4, tau_test*30)
      # sol_tau_18 <- deSolve::ode(
      #   y = c(h = 13, n = 0, x = 1, y = 0, S = 0),
      #   times = times,
      #   func = "derivsc_p",
      #   parms = unlist(hclParms(tau_c=18)),
      #   dllname = "rHCL",
      #   initforc = "forcc_p",
      #   forcings = cbind(times, rep(0, length(times))),
      #   fcontrol = list(method="linear", rule=2, f=0),
      #   initfunc = "parmsc_p",
      #   nout = 0,
      #   # events = list(func="eventc_p", root=TRUE),
      #   # rootfun = "rootc_p",
      #   # nroot = 1
      # )
    },
    finally = {
      ## Clean up the function I added to the global environment ##
      if(exists("light.int", where = .GlobalEnv)){
        rm(light.int, envir = .GlobalEnv)
      } # remove light approxfun

      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )
  browser()

})


### TODO - NEED TO TEST THAT THE ODEs ARE RETURNING VALUES THAT ARE CORRECT



# # Test C code for differential equations ----------------------------------
# test_that("ode c code correctly works (ignoring roots)", {
#
#   # check if .dll is loaded, load if needed (will unload after iterations)
#   if(!is.loaded("src/rHCL.dll")){
#     dyn.load("src/rHCL.dll")
#   }
#
#   # times <- seq(from=0, to=24, by = 1) # time points to evaluate
#   times = c(0, 12, 24, 24.2, 48, 48.4, 72, 72.6, 96, 96.8, 120, 121, 144, 145.2)
#
#
#   # using tryCatch to make sure dyn.unload() runs even if there is an error
#   tryCatch(
#     sol <- deSolve::ode(
#       y = c(h = 13, n = 0, x = 1, y = 0, S = 0),
#       times = times,
#       func = "derivsc_p",
#       parms = unlist(hclParms()),
#       dllname = "rHCL",
#       initforc = "forcc_p",
#       forcings = cbind(times, rep(0, length(times))),
#       fcontrol = list(method="linear", rule=2, f=0),
#       initfunc = "parmsc_p",
#       nout = 0,
#       # events = list(func="eventc_p", root=TRUE),
#       # rootfun = "rootc_p",
#       # nroot = 1
#     ),
#     finally = {
#       # unload .dll
#       dyn.unload("src/rHCL.dll")
#     }
#   )
#
#   ttimes <- atan2(sol[,"y"], sol[,"x"]) # calcuate radians
#   ttimes[ttimes < 0] <- (2*pi) + ttimes[ttimes<0] # adjust negative values
#   sol2 <- cbind(sol, ttimes*12/pi) # convert to time
#
#   browser()
#
#
#     ## TODO - test that sleep pressure builds up as expected based on equation 1 (in the absence of switching)
#
# })

## TODO Test that different time scales lead to same results
