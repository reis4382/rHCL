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

  ## create times
  times = c(0, 12, 24, 24.2, 48, 48.4, 72, 72.6, 96, 96.8, 120, 121, 144, 145.2)

  the$light_int <- stats::approxfun(x=times, y=rep(0, length(times)), method="linear", rule=2) # update package custom ("the") environment

  # ## Create light interpolation function for R code ##
  # assign("light.int",
  #        stats::approxfun(x=times, y=rep(0, length(times)), method="linear", rule=2),
  #        envir = .GlobalEnv) # create interpolation function; Note this is creating a global environment and needs to be cleaned up
  # ## TODO - Is there a way to make light.int() available for testing but not place it in the global env?
  #
  ## check if .dll is loaded, load if needed (will unload after test)
  ## TODO - is there a better way of loading c code functions for testing?
  # if(!"rHCL" %in% names(getLoadedDLLs())){
  #   # using here package to find root of rstudio project directory b/c
  #   # when running test suite the working directory switches to test folder
  #   dyn.load(paste(here::here(), "src/rHCL.dll", sep = "/"))
  # }


  # first test - NO LIGHT
  sol_r <- deSolve::ode(y = c(h = 13, n = 0, x = 1, y = 0, S = 0),
                        func = dHCL,
                        times = times,
                        parms = hclParms(),
                        events = list(func = dEventFunc, root = TRUE),
                        rootfun =  dRootFunc)


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
  the$light_int <- stats::approxfun(x=times, y=light, method="linear", rule=2)
  #
  # assign("light.int",
  #        stats::approxfun(x=times, y=light, method="linear", rule=2),
  #        envir = .GlobalEnv) # create interpolation function; Note this is creating a global environment and needs to be cleaned up

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



  expect_equal(as.numeric(sol_r), as.numeric(sol_c))
  expect_equal(as.numeric(sol_r2), as.numeric(sol_c2))
  expect_equal(as.numeric(sol_r3), as.numeric(sol_c3))

})


# Test that ODEs return expected values -----------------------------------
test_that("ODE output in the absence of light", {

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

  ### Tests ###
  expect_equal(sol_hours[,!colnames(sol_hours) %in% "time"], sol_mins[,!colnames(sol_mins) %in% "time"])
  expect_equal(sol_mins[,!colnames(sol_mins) %in% "time"], sol_secs[,!colnames(sol_secs) %in% "time"])

})


# TODO create tests for forced wake addition ------------------------------
test_that("Forced wake functions operate correctly with forced wake input", {

  ## create times
  times = seq(from = 0, to = 60, by = .1)
  light <- rep(0, length(times)) # light vector
  light[(times %% 24) > 8 & (times %% 24) < 22] <- 1000 # 1000 lux exposure from 8 am - 10 pm
  f_wake <- rep(0, length(times)) # forced wake vector
  f_wake[(times %% 24) >= 2 & (times %%24) < 3] <- 1 # force wake between 2 and 3 am

  ## Create light interpolation function for R code ##
  # NOTE: C code will require constant interpolation for both (can't have different methods)
  the$light_int <- stats::approxfun(x=times, y=light, method="constant", rule=2) # update package custom ("the") environment
  the$force_wake <- stats::approxfun(x=times, y=f_wake, method="constant", rule=2) # update package custom ("the") environment

  ## first test - R code ##
  # Function w force wake
  sol_r <- deSolve::ode(y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
                        func = dHCL,
                        times = times,
                        parms = hclParms(),
                        events = list(func = dEventFunc, root = TRUE),
                        rootfun = dRootFunc_FW)

  # test that no enforced wake times are sleeping
  expect_equal(sum(sol_r[(sol_r[,"time"] %% 24) >= 2 & (sol_r[,"time"] %% 24) <3 ,"S"] !=0), 0)

  ## C code ##
  sol_c <- deSolve::ode(y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0, fwake = 0),
                        times = times,
                        func = "derivsc_p_fw",
                        parms = unlist(hclParms()),
                        dllname = "rHCL",
                        initforc = "forcc_p_fw",
                        forcings = list(cbind(times, light),
                                        cbind(times, f_wake)),
                        fcontrol = list(method = "constant", rule=2, f=0),
                        initfunc = "parmsc_p",
                        nout = 0,
                        events = list(func = "eventc_p", root = TRUE),
                        rootfun = "rootc_p_fw",
                        nroot = 1)

  sol_c <- as.data.frame(sol_c)

  # test that no enforced wake times are sleeping

  # NOTE: due to the hacky way fwake is implemented in C code, the first epoch with
  # fwake == 1 won't be sufficient to trigger the increased circadian threshold.
  # For example, have fwake == 1 start at 3 am (during normal sleep) won't
  # result in spontaneous forced wake until the first epoch after 3 am.
  # A hacky way to solve this may be to cause the epoch before the first actual
  # fwake to also fwake == 1. Similarly, may need to end fwake an epoch early
  # to have the threshold lower back down.

  # therefore, only times
  expect_equal(sum(sol_c[(sol_c[,"time"] %% 24) > 2 & (sol_c[,"time"] %% 24) <3 ,"S"] !=0), 0)

  ## To get sol_c to match sol_r (which is what I want), let's try adjusting the
  # force wake epochs #
  f_wake2 <- f_wake
  f_wake_inds <- which((times%%24)==2) # indices in f_wake equal to 2 am
  f_wake2[f_wake_inds-1] <- 1 # assign preceding index to be 1
  f_wake_inds2 <- which((times%%24)==3)
  f_wake2[f_wake_inds2-1] <- 0 # assign preceding index to be 0

  sol_c2 <- deSolve::ode(y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0, fwake = 0),
                        times = times,
                        func = "derivsc_p_fw",
                        parms = unlist(hclParms()),
                        dllname = "rHCL",
                        initforc = "forcc_p_fw",
                        forcings = list(cbind(times, light),
                                        cbind(times, f_wake2)),
                        fcontrol = list(method = "constant", rule=2, f=0),
                        initfunc = "parmsc_p",
                        nout = 0,
                        events = list(func = "eventc_p", root = TRUE),
                        rootfun = "rootc_p_fw",
                        nroot = 1)

  sol_c2 <- as.data.frame(sol_c2)
  sol_c2 <- sol_c2[,!names(sol_c2) %in% "fwake"] # this derivative column won't be that useful

  # due to the hacky approach for C, sol_r and sol_c2 will not be identical w/ respect to
  # states, as the transition to sleep for sol_c2 is likely taking place betwen observed
  # times. Hopefully sleep states are the same though, at least for the most part.
  sol_r <- as.data.frame(sol_r)
  expect_equal(sol_r$S, sol_c2$S)

  h_diff <- mean(sol_r$h - sol_c2$h)
  expect_equal(abs(h_diff) < 1e-2, TRUE) # h_diff is the most different
  # because differences in transitions to sleep/wake directly impact
  # sleep pressure build up. The effect of transition differences on
  # x, y, and n will be less, because they are only affected by the
  # sleep gating of light.

  x_diff <- mean(sol_r$x - sol_c2$x)
  expect_equal(abs(x_diff) < 1e-3, TRUE)

  y_diff <- mean(sol_r$y - sol_c2$y)
  expect_equal(abs(y_diff) < 1e-3, TRUE)

  n_diff <- mean(sol_r$n - sol_c2$n)
  expect_equal(abs(n_diff) < 1e-3, TRUE)

})

test_that("forced wake code returns same results if no forced wake", {

  ## create times
  times = seq(from = 0, to = 60, by = .1)
  light <- rep(0, length(times)) # light vector
  light[(times %% 24) > 8 & (times %% 24) < 22] <- 1000 # 1000 lux exposure from 8 am - 10 pm
  f_wake <- rep(0, length(times)) # forced wake vector

  ## Create light interpolation function for R code ##
  # NOTE: C code will require constant interpolation for both (can't have different methods)
  the$light_int <- stats::approxfun(x=times, y=light, method="constant", rule=2) # update package custom ("the") environment
  the$force_wake <- stats::approxfun(x=times, y=f_wake, method="constant", rule=2) # update package custom ("the") environment

  ## first test - R code ##
  # Function w force wake
  sol_r1 <- deSolve::ode(y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
                        func = dHCL,
                        times = times,
                        parms = hclParms(),
                        events = list(func = dEventFunc, root = TRUE),
                        rootfun = dRootFunc_FW)

  sol_r2 <- deSolve::ode(y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
                         func = dHCL,
                         times = times,
                         parms = hclParms(),
                         events = list(func = dEventFunc, root = TRUE),
                         rootfun = dRootFunc)

  expect_equal(sol_r1, sol_r2)

  ## C code ##
  sol_c1 <- deSolve::ode(y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0, fwake = 0),
                        times = times,
                        func = "derivsc_p_fw",
                        parms = unlist(hclParms()),
                        dllname = "rHCL",
                        initforc = "forcc_p_fw",
                        forcings = list(cbind(times, light),
                                        cbind(times, f_wake)),
                        fcontrol = list(method = "constant", rule=2, f=0),
                        initfunc = "parmsc_p",
                        nout = 0,
                        events = list(func = "eventc_p", root = TRUE),
                        rootfun = "rootc_p_fw",
                        nroot = 1)

  sol_c2 <- deSolve::ode(y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
                         times = times,
                         func = "derivsc_p",
                         parms = unlist(hclParms()),
                         dllname = "rHCL",
                         initforc = "forcc_p",
                         forcings = list(cbind(times, light)),
                         fcontrol = list(method = "constant", rule=2, f=0),
                         initfunc = "parmsc_p",
                         nout = 0,
                         events = list(func = "eventc_p", root = TRUE),
                         rootfun = "rootc_p",
                         nroot = 1)

  # compare after dropping extra "fwake" column
  sol_c1 <- as.data.frame(sol_c1[,-ncol(sol_c1)])
  sol_c2 <- as.data.frame(sol_c2)

  expect_equal(sol_c1, sol_c2)

  # compare c and r results when no forced wake
  sol_r1 <- as.data.frame(sol_r1)

  expect_equal(sol_r1, sol_c1)

})

test_that("rapid forced wake transitions work for r and c code", {

  times <- c(0, .1, .2, .3, .4, .5, .6)
  light <- rep(0, length(times))
  f_wake <- c(0, 1, 0, 1, 0, 1, 0)

  ## Create light interpolation function for R code ##
  # NOTE: C code will require constant interpolation for both (can't have different methods)
  the$light_int <- stats::approxfun(x=times, y=light, method="constant", rule=2) # update package custom ("the") environment
  the$force_wake <- stats::approxfun(x=times, y=f_wake, method="constant", rule=2) # update package custom ("the") environment

  ## first test - R code ##
  # Function w force wake
  sol_r <- deSolve::ode(y = c(h = 14, n = .152, x = -0.966, y = -0.558, S = 1),
                         func = dHCL,
                         times = times,
                         parms = hclParms(),
                         events = list(func = dEventFunc, root = TRUE),
                         rootfun = dRootFunc_FW)
  sol_r <- as.data.frame(sol_r)

  expect_sleep <- c(1, 0, 1, 0, 1, 0, 1) # expected sleep pattern
  expect_equal(sol_r$S, expect_sleep)

  ## C code ##
  c_f_wake <- c(1, 0, 1, 0, 1, 0, 0) # need to shift epochs forward - note
  # that involves altering the starting state of fwake if needed

  sol_c <- deSolve::ode(y = c(h = 14, n = .152, x = -0.966, y = -0.558, S = 1, fwake = 0),
                        times = times,
                        func = "derivsc_p_fw",
                        parms = unlist(hclParms()),
                        dllname = "rHCL",
                        initforc = "forcc_p_fw",
                        forcings = list(cbind(times, light),
                                        cbind(times, c_f_wake)),
                        fcontrol = list(method = "constant", rule=2, f=0),
                        initfunc = "parmsc_p",
                        nout = 0,
                        events = list(func = "eventc_p", root = TRUE),
                        rootfun = "rootc_p_fw",
                        nroot = 1)
  sol_c <- as.data.frame(sol_c)

  expect_equal(sol_c$S, expect_sleep)


})

test_that("forced wake c code works with small time steps", {

  ## create times
  times = seq(from = 0, to = 60, by = .5/60)
  light <- rep(0, length(times)) # light vector
  light[(times %% 24) > 8 & (times %% 24) < 22] <- 1000 # 1000 lux exposure from 8 am - 10 pm
  f_wake <- rep(0, length(times)) # forced wake vector
  f_wake[(times %% 24) >= 2 & (times %%24) < 3] <- 1 # force wake between 2 and 3 am

  ## C code ##
  sol_c <- deSolve::ode(y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0, fwake = 0),
                        times = times,
                        func = "derivsc_p_fw",
                        parms = unlist(hclParms()),
                        dllname = "rHCL",
                        initforc = "forcc_p_fw",
                        forcings = list(cbind(times, light),
                                        cbind(times, f_wake)),
                        fcontrol = list(method = "constant", rule=2, f=0),
                        initfunc = "parmsc_p",
                        nout = 0,
                        events = list(func = "eventc_p", root = TRUE),
                        rootfun = "rootc_p_fw",
                        nroot = 1)

  sol_c <- as.data.frame(sol_c)

  # NOTE: due to the hacky way fwake is implemented in C code, the first epoch with
  # fwake == 1 won't be sufficient to trigger the increased circadian threshold.
  # For example, have fwake == 1 start at 3 am (during normal sleep) won't
  # result in spontaneous forced wake until the first epoch after 3 am.
  # A hacky way to solve this may be to cause the epoch before the first actual
  # fwake to also fwake == 1. Similarly, may need to end fwake an epoch early
  # to have the threshold lower back down.

  # therefore, only times
  expect_equal(sum(sol_c[(sol_c[,"time"] %% 24) > 2 & (sol_c[,"time"] %% 24) <3 ,"S"] !=0), 0)

})

# Test that Forger 1999 odes work -----------------------------------------

test_that("dForger ODE functions work", {

  ## create times
  times = c(0, 12, 24, 24.2, 48, 48.4, 72, 72.6, 96, 96.8, 120, 121, 144, 145.2)

  ## Create light interpolation function for R code ##
  the$light_int <- stats::approxfun(x=times, y=rep(0, length(times)), method="linear", rule=2) # update package custom ("the") environment

  # first test - NO LIGHT
  sol_r <- deSolve::ode(y = c(n = 0, x = 1, y = 0),
                        func = dForger,
                        times = times,
                        parms = hclParms())

  sol_c <- deSolve::ode(
    y = c(n = 0, x = 1, y = 0),
    times = times,
    func = "derivsc_forger",
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(times, rep(0, length(times))),
    fcontrol = list(method="linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0
  )

  ### second test - Light exposure (no roots) ###
  times <- seq(from=0, to=24.2, by = .1) # start times
  light <- rep(0, length(times)) # start light vector
  light[times > 8 & times < 22] <- 1000 # light exposure during "day"
  # re-assign light func
  the$light_int <- stats::approxfun(x=times, y=light, method="linear", rule=2) # update package custom ("the") environment

  sol_r2 <- deSolve::ode(y = c(n = 0, x = 1, y = 0),
                         func = dForger,
                         times = times,
                         parms = hclParms())

  sol_c2 <- deSolve::ode(
    y = c(n = 0, x = 1, y = 0),
    times = times,
    func = "derivsc_forger",
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(times, light),
    fcontrol = list(method="linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0
  )

  #### Third test - extended days to establish convergence ####
  times <- seq(from=0, to=24*50, by = .1) # start times
  light <- rep(0, length(times)) # start light vector
  light[(times%%24) > 6 & (times%%24) < 22] <- 1000 # light exposure during "day"

  # re-assign light func
  the$light_int <- stats::approxfun(x=times, y=light, method="linear", rule=2) # update package custom ("the") environment

  sol_c3 <- deSolve::ode(
    y = c(n = 0, x = 1, y = 0),
    times = times,
    func = "derivsc_forger",
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(times, light),
    fcontrol = list(method="linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0
  )

  ## test that r and c code are identical ##
  expect_equal(as.numeric(sol_r), as.numeric(sol_c))
  expect_equal(as.numeric(sol_r2), as.numeric(sol_c2))

  ## Test that results are as expected ##

  res_sub <- as.data.frame(sol_c3[(1+24*49*10):nrow(sol_c3), ]) # extract final day
  min_time <- res_sub$time[which.min(res_sub$y)] %% 24

  expect_equal(min_time, 3.5)

})

test_that("deSolve::ode() works with example data", {

  # been getting some weird, unreproducible deSolve errors, where the step size
  # is 0. Seems to happen after I run the package examples, after which any
  # call to ode() fails. Restarting R seems to help.
  check_df <- dfPrep(df = rhcl_df, time_var = "times", light_var = "light")

  expect_no_error(deSolve::ode(
    y = c(n = 0, x = 1, y = 0),
    times = check_df$ctime,
    func = "derivsc_forger",
    parms = unlist(hclParms()),
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(check_df$ctime, check_df$light),
    fcontrol = list(method="linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0
  ))

})

test_that("Ensure C code doesn't generate NaNs in edge case with negative light", {

  # edge case in real data, where it appeared that light was being treated as a
  # negative value (likely floating point issue) and somehow that was causing
  # n, x, and y to become NaN. This is because the pow() function in C cannot
  # take a negative base if the exponent is not an integer, and will return NaN.

  times <- c(126.3 - (2/60), 126.3 - (1/60), 126.3, 126.3 + (1/60)) # isolating an issue from real data
  light <- c(12.32, 9.24, 0, 0)

  sol <- deSolve::ode(
    # state right before NaN issue occurs
    # y = c(h = 11.1346787402, n = 0.0006959234, x = 0.9053542840, y = -0.6286116919, S = 1),
    # state two before NaN issue
    y = c(h = 11.1388, n = 0.000705029, x = 0.9027115, y = -0.6325563, S = 1),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms(mu = 14.0499525437)), # mu at which error occurs
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

  ## explicitly include negative light values - note: full set of functions will
  # prepare data and cause error if negative light values are provided.
  light2 <- c(12, -1, -1, -1)
  sol2 <- deSolve::ode(
    # state right before NaN issue occurs
    # y = c(h = 11.1346787402, n = 0.0006959234, x = 0.9053542840, y = -0.6286116919, S = 1),
    # state two before NaN issue
    y = c(h = 11.1388, n = 0.000705029, x = 0.9027115, y = -0.6325563, S = 1),
    times = times,
    func = "derivsc_p",
    parms = unlist(hclParms(mu = 14.0499525437)), # mu at which error occurs
    dllname = "rHCL",
    initforc = "forcc_p",
    forcings = cbind(times, light2),
    fcontrol = list(method="linear", rule=2, f=0),
    initfunc = "parmsc_p",
    nout = 0,
    events = list(func="eventc_p", root=TRUE),
    rootfun = "rootc_p",
    nroot = 1
  )

  expect_equal(sum(is.na(sol)), 0)
  expect_equal(sum(is.na(sol2)), 0)
})

