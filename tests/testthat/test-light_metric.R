# Tests of novel light metric ---------------------------------------------

test_that("circLight() works", {

  tryCatch(
    {
      ## check if .dll is loaded, load if needed (will unload after test)
      ## TODO - is there a better way of loading c code functions for testing?
      if(!"rHCL" %in% names(getLoadedDLLs())){
        # using here package to find root of rstudio project directory b/c
        # when running test suite the working directory switches to test folder
        dyn.load(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }

      ## generate test profile - initial state conditions are those at midnight of someone
      ## with midsleep of 4:30 and sleep duration of 8.22 under typical Skeldon 2017 light profile
      # TODO cannot achieve convergence that matches those sleep outcomes under the default light profile.
      # probably need to reach out to authors.

      # says to use default parameters

      times <- seq(0, 24*30, by = 1/60) # 12-minute intervals
      light <- lightCycle(times) # generate light profile in skeldon 2017 paper (see function documentation for ref)

      desolve_list <- list(
        y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
        times = times,
        func = "derivsc_p",
        parms = unlist(hclParms(mu = 19.18, tau_c = 24.3271)),
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

      sol <- odeIter(desolve_args = desolve_list, max_iter = 40, dur_tol=1/60, mid_tol = 5/60)
      browser()
      sol$sleep_sum



    },
    finally = {

      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
  )


})
