## code to prepare dataset goes here ##

## Skeldon 2017 paper default light values ##
## Per personal correspondence w/ Prof. Skeldon, midsleep timing should be 3:25 am
## when using s1 = 8, s2 = 17 or 2:55 am when using s1 = 7.5, s2 = 16.5

times <- seq(0, 24*14-(1/60), by = 1/60) # 1-minute intervals

start_dtime <- lubridate::ymd_hms("2025-01-01 00:00:00", tz = "America/Denver")
dtimes <- start_dtime + times * 60 * 60 # POSIXct format stamps
light <- lightCycle(times*60*60, l1=700, l2=40, c=.6, s1 = 8, s2 = 17, time_scale = "secs")

## Time and light data ##
rhcl_df <- data.frame(
  times = dtimes,
  light = light
)

## generate spontaneous sleep values for default light profile
desolve_list <- list(
  # initial values, arbitrary
  y = c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0),
  times = times, # desolve needs numbers
  func = "derivsc_p",
  parms = unlist(hclParms(mu = 17.1, tau = 24.22)), # set desired mu and tau values
  dllname = "rHCL",
  initforc = "forcc_p",
  forcings = cbind(times, rhcl_df$light),
  fcontrol = list(method = "linear", rule=2, f=0),
  initfunc = "parmsc_p",
  nout = 0,
  events = list(func = "eventc_p", root = TRUE),
  rootfun = "rootc_p",
  nroot = 1
)

# extract sleep duration summary of synthetic data #
syn_sol <- odeIter(desolve_args = desolve_list, dtime_vec = rhcl_df$times,
                   max_iter = 20, dur_tol = 1/60, mid_tol = 1/60,
                   epoch_length_min = 1, min_observed_hours = 18)

## merge sleep data ##
rhcl_df$sleep <- syn_sol$ode_res$S


usethis::use_data(rhcl_df, overwrite = TRUE) # write data
