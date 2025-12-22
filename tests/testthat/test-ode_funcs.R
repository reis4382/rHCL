# Tests of hclParms(), which generates ODE parameters -----------------------

test_that("hclParms() returns correct format and values", {
  expect_equal(class(hclParms()), "list") # test that a list is returned
  expect_equal(length(hclParms()), 16) # test length of returned list
  expect_equal(names(hclParms()), c("mu", "chi", "Hzero", "delta", "ca_par",
                                    "tau_c", "f_par", "G_par", "p_par", "k_par",
                                    "little_b_par", "gamma", "alpha_zero",
                                    "beta", "Izero", "time_scale")) # correct names are returned
  # manually check all terms
  expect_equal(as.numeric(hclParms()),
               c(17.87, 45, 13, 1, 1.72, 24.2, 0.99669, 19.9, 0.6, 0.55, 0.4, 0.23, 0.16, 0.013, 9500, 1))
})

test_that("hclParms() correctly accepts argument changes", {
  expect_equal(hclParms(tau_c = 24)[["tau_c"]], 24)
  expect_equal(as.numeric(hclParms(mu = 15, gamma = .5)[c("mu", "gamma")]), c(15, 0.5))
})

test_that("hclParms() correctly accepts time_scale arguments", {
  expect_equal(as.numeric(hclParms(time_scale="hours")[c("chi", "time_scale")]), c(45, 1)) # scale to hours
  expect_equal(as.numeric(hclParms(time_scale="mins")[c("chi", "time_scale")]), c(45*60, 60)) # scale to minutes
  expect_equal(as.numeric(hclParms(time_scale="secs")[c("chi", "time_scale")]), c(45*60*60, 60*60)) # scale to seconds

  expect_equal(as.numeric(hclParms(time_scale=1)[c("chi", "time_scale")]), c(45, 1)) # scale to hours
  expect_equal(as.numeric(hclParms(time_scale=1/60)[c("chi", "time_scale")]), c(45*60, 60)) # scale to minutes
  expect_equal(as.numeric(hclParms(time_scale=1/60/60)[c("chi", "time_scale")]), c(45*60*60, 60*60)) # scale to seconds

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

      # tau 24.2
      ## create times
      times_24pt2 <- seq(0, 49, by = .2)

      sol_tau_24pt2 <- deSolve::ode(
        y = c(h = 5, n = 0, x = 0, y = -1, S = 0),
        times = times_24pt2,
        func = "derivsc_p",
        parms = unlist(hclParms()),
        dllname = "rHCL",
        initforc = "forcc_p",
        forcings = cbind(times_24pt2, rep(0, length(times_24pt2))),
        fcontrol = list(method="linear", rule=2, f=0),
        initfunc = "parmsc_p",
        nout = 0,
        # events = list(func="eventc_p", root=TRUE),
        # rootfun = "rootc_p",
        # nroot = 1
      )

      # # tau 18
      # ## create times
      times_18 <- seq(0, 55, by = .2)

      sol_tau_18 <- deSolve::ode(
        y = c(h = 5, n = 0, x = 0, y = -1, S = 0),
        times = times_18,
        func = "derivsc_p",
        parms = unlist(hclParms(tau_c=18)),
        dllname = "rHCL",
        initforc = "forcc_p",
        forcings = cbind(times_18, rep(0, length(times_18))),
        fcontrol = list(method="linear", rule=2, f=0),
        initfunc = "parmsc_p",
        nout = 0,
        # events = list(func="eventc_p", root=TRUE),
        # rootfun = "rootc_p",
        # nroot = 1
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

  ### Tests ###
  ## 24.2 tau mins - test that yminimums occur on tau period lengths ##
  tau24pt2_mins <- minMaxFinder(sol_tau_24pt2[,'y'])
  # allow tolerance - seq() has rounding error
  expect_equal(tau24pt2_mins, which((times_24pt2 %% 24.2) > -1e-10 & (times_24pt2 %% 24.2) < 1e-10)[-1])

  ## 18 tau mins - test that yminimums occur on tau period lengths##
  tau18_mins <- minMaxFinder(sol_tau_18[,"y"])
  expect_equal(tau18_mins, which((times_18 %% 18) > -1e-10 & (times_18 %% 18) < 1e-10)[-1])

  ## test homeostatic sleep pressure accumulation in absence of sleep ##
  tau_24pt2_pressure <- sleepHomeostasis(
    hclParms()[["mu"]], tswitch=0, h_tswitch=5, time=times_24pt2,
    chi = hclParms()[["chi"]], s=0
    )

  tau_18_pressure <- sleepHomeostasis(
    hclParms()[["mu"]], tswitch=0, h_tswitch=5, time=times_18,
    chi = hclParms()[["chi"]], s=0
  )

  expect_equal(sol_tau_24pt2[,"h"], tau_24pt2_pressure)
  expect_equal(sol_tau_18[,"h"], tau_18_pressure)


  ## test that photoreceptor fraction never changes ##
  expect_equal(sum(sol_tau_24pt2[,"n"]==0), nrow(sol_tau_24pt2))
  expect_equal(sum(sol_tau_18[,"n"]==0), nrow(sol_tau_18))
})


# Test ODE output when exposed to light -----------------------------------

test_that("Regular light leads to expected 24-hour period once entrained",{

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
      light <- rep(0, length(times)) # light vector
      light[(times %% 24) > 8 & (times %%24) < 22] <- 1000 # 1000 lux exposure from 8 am - 10 pm

      sol <- deSolve::ode(
        # initial values, should correspond to a 4:24 am ymin time and ~7.4 hours of nightly sleep
        # NOTE: changing the model would likely change this, so tests will be based on post-entrainment estimates
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

  ### tests ###
  ## Allow entrainment in beginning, so extract ~final week
  sol2 <- sol[sol[,"time"] > 24*23, ]
  ymin_times <- times[minMaxFinder(sol2[,"y"])] %% 24 # daily times of ymin
  expect_equal(var(ymin_times), 0) # times of day should be the same (no variance)

  ## check that sleep pressure is changing as expected ##
  # separate sleep and wake runs, starting first period at beginning of data
  sleep_runs <- as.data.frame(unclass(rle(sol[,"S"]))) # create run data.frame
  sleep_runs$ends <- cumsum(sleep_runs$lengths) # ends of runs
  sleep_runs$starts <- sleep_runs$ends - sleep_runs$lengths + 1 # starts of runs

  # prepare data frame for calculates
  h_df <- as.data.frame(sol) # convert to df to help format
  h_df$tswitch <- NA # initialize
  h_df$h_tswitch <- NA # initialize

  # iterate through runs and
  for(i in 1:nrow(sleep_runs)){
    # extract the time of the sleep state start
    h_df[sleep_runs$starts[i]:sleep_runs$ends[i], "tswitch"] <- sol[sleep_runs$starts[i], "time"]
    # extract the sleep pressure of the sleep state start
    h_df[sleep_runs$starts[i]:sleep_runs$ends[i], "h_tswitch"] <- sol[sleep_runs$starts[i], "h"]
  }

  # calculate sleep pressures (might be a bit inefficient) #
  new_h <- apply(h_df, 1, function(x){
    # calculate sleep pressure (row-by-row, a little inefficient)
    new_h <- sleepHomeostasis(mu = hclParms()[["mu"]],
                              tswitch = x[["tswitch"]],
                              h_tswitch = x[["h_tswitch"]],
                              time = x[["time"]],
                              chi = hclParms()[["chi"]],
                              s = x[["S"]])
    return(new_h)
  })

  # browser()
  expect_equal(h_df$h, new_h)
})



# Test that different time scales lead to same results --------------------
test_that("Changing time scales does not affect results", {

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

      ## Create time variables ##
      times_secs <- seq(0, 1440*60, by = 6*60) # 6-minute intervals in seconds
      times_mins <- times_secs / 60
      times_hours <- times_mins / 60
      light_vec <- rep(0, length(times_secs))
      light_vec[times_hours > 8 & times_hours < 22] <- 5000 # set daytime light

      sol_hours <- deSolve::ode(
        y = c(h = 13, n = 0, x = 0, y = -1, S = 0),
        times = times_hours,
        func = "derivsc_p",
        parms = unlist(hclParms(time_scale = "hours")),
        dllname = "rHCL",
        initforc = "forcc_p",
        forcings = cbind(times_hours, light_vec),
        fcontrol = list(method="linear", rule=2, f=0),
        initfunc = "parmsc_p",
        nout = 0,
        events = list(func="eventc_p", root=TRUE),
        rootfun = "rootc_p",
        nroot = 1
      )


      # times_mins <- seq(0, 1250, by = 6)
      # light_vec <- light_vec[1:length(times_mins)]

      sol_mins <- deSolve::ode(
        y = c(h = 13, n = 0, x = 0, y = -1, S = 0),
        times = times_mins,
        func = "derivsc_p",
        parms = unlist(hclParms(time_scale = "mins")),
        dllname = "rHCL",
        initforc = "forcc_p",
        forcings = cbind(times_mins, light_vec),
        fcontrol = list(method="linear", rule=2, f=0),
        initfunc = "parmsc_p",
        nout = 0,
        events = list(func="eventc_p", root=TRUE),
        rootfun = "rootc_p",
        nroot = 1
      )
      # browser()

      sol_secs <- deSolve::ode(
        y = c(h = 13, n = 0, x = 0, y = -1, S = 0),
        times = times_secs,
        func = "derivsc_p",
        parms = unlist(hclParms(time_scale = "secs")),
        dllname = "rHCL",
        initforc = "forcc_p",
        forcings = cbind(times_secs, light_vec),
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

  # browser()

  ### Tests ###
  expect_equal(sol_hours[,!colnames(sol_hours) %in% "time"], sol_mins[,!colnames(sol_mins) %in% "time"])
  expect_equal(sol_mins[,!colnames(sol_mins) %in% "time"], sol_secs[,!colnames(sol_secs) %in% "time"])

})


# TODO create tests for forced wake addition ------------------------------
test_that("Forced wake functions operate correctly with forced wake input", {
  ## using tryCatch() so that light.int() and dyn.load are always removed,
  ## even if the test errors out

  ## I can't get forced wake to work in compiled code. Seems to be a possible but
  ## in the root function - the roots are being weird about interpolating the
  ## forced wake values. They seem to be failing to reset it if around an event?
  ## R code seems to work though.

  tryCatch(
    {
      ## create times
      times = seq(from = 0, to = 60, by = .1)
      light <- rep(0, length(times)) # light vector
      light[(times %% 24) > 8 & (times %% 24) < 22] <- 1000 # 1000 lux exposure from 8 am - 10 pm
      f_wake <- rep(0, length(times)) # forced wake vector
      f_wake[(times %% 24) > 2 & (times %%24) < 3] <- 1 # force wake between 2 and 3 am

      ## Create light interpolation function for R code ##
      # NOTE: C code will require constant interpolation for both (can't have difference methods)
      assign("light.int",
             approxfun(x=times, y=light, method="constant", rule=2),
             envir = .GlobalEnv) # create interpolation function; Note this is creating a global environment and needs to be cleaned up

      ## create force wake interpolation function ##
      assign("force.wake",
             approxfun(x=times, y=f_wake, method="constant", rule=2),
             envir = .GlobalEnv) # create interpolation function; Note this is creating a global environment and needs to be cleaned up

      ## TODO - Is there a way to make light.int() available for testing but not place it in the global env?


      ## first test - R code ##
      # Function w force wake
      sol_r <- deSolve::ode(y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
                             func = dHCL,
                             times = times,
                             parms = hclParms(),
                             events = list(func = dEventFunc, root = TRUE),
                             rootfun = dRootFunc_FW)



    },
    finally = {
      ## Clean up the function I added to the global environment ##
      if(exists("light.int", where = .GlobalEnv)){
        rm(light.int, envir = .GlobalEnv)
      } # remove light approxfun

      if(exists("force.wake", where = .GlobalEnv)){
        rm(force.wake, envir = .GlobalEnv)
      } # remove force wake approxfun
    }
  )

  # test that no enforced wake times are sleeping
  expect_equal(sum(sol_r[(sol_r[,"time"] %% 24) > 2 & (sol_r[,"time"] %% 24) <3 ,"S"] !=0), 0)

})
