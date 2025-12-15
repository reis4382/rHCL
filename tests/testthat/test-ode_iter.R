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

      sol <- odeIter(desolve_args = desolve_list, max_iter = 20)

      ## No light exposure ##
      desolve_list2 <- desolve_list
      desolve_list2[["forcings"]] <- cbind(times, rep(0, length(times)))

      sol2 <- odeIter(desolve_args = desolve_list2, max_iter = 21)

    },
    finally = {

      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )

  ## Check - model converges when sufficient light ##
  expect_equal(sol$converge, TRUE)
  expect_match(sol$conv_message, regexp = "Convergence obtained after 3 iterations.")

  ## Check - model does not converge when no light ##
  expect_equal(sol2$converge, FALSE)
  expect_match(sol2$conv_message, regexp = "The model did not converge.")
  expect_equal(sol2$iterations, 21)
})
