
# Optimization control functions ------------------------------------------

test_that("durationOptControl() correctly formats arguments", {

  durationOptControl(maximum = TRUE) # test optional optimize() argument

})


# rhcl() tests ------------------------------------------------------------

test_that("rhcl() correctly optimizes parameters", {


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
      times = seq(0, 24*7, by = 1/60) # 1-minute intervals

      df <- data.frame(
        times = times,
        light = lightCycle(times, l1=1000, l2=5) # generate standard light profile
      )


      ## set up a synthetic sleep wake cycle ##

      # create a list for deSolve::ode arguments #
      desolve_list <- list(
        # initial values, arbitrary
        y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
        times = df$times,
        func = "derivsc_p",
        parms = unlist(hclParms(mu = 16.5, tau = 24.5)), # set desired mu and tau values
        dllname = "rHCL",
        initforc = "forcc_p",
        forcings = cbind(df$times, df$light),
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

      res <- rhcl(df = df,
                  time_var = "times",
                  sleep_var = "sleep",
                  light_var = "light",
                  y0 = NULL,
                  ode_parms = hclParms(),
                  sleep_mid = syn_mid,
                  sleep_dur = syn_dur,
                  max_ode_iter = 20,
                  dur_tol = 1/60,
                  mid_tol = 1/60,
                  compiled = TRUE,
                  opt_method = "bisect"
                  )

      res2 <- rhcl(df = df,
                  time_var = "times",
                  sleep_var = "sleep",
                  light_var = "light",
                  y0 = NULL,
                  ode_parms = hclParms(),
                  sleep_mid = syn_mid,
                  sleep_dur = syn_dur,
                  max_ode_iter = 20,
                  dur_tol = 1/60,
                  mid_tol = 1/60,
                  compiled = TRUE,
                  opt_method = "optimize"
      )

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


    },
    finally = {

      # unload .dll
      if("rHCL" %in% names(getLoadedDLLs())){
        dyn.unload(paste(here::here(), "src/rHCL.dll", sep = "/"))
      }
    }
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
})

