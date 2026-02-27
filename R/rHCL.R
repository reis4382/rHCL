#' Control of sleep duration parameter (\eqn{\mu}) optimization
#'
#' Creates a list of arguments that will be used to control the optimization
#' of \eqn{\mu}.
#'
#' @param param_lower Lower boundary of parameter for estimation (numeric). Default is NULL,
#' which will use the minimum allowable value as described by Skeldon et al. 2023.
#' If using the bisection method, this will be the lower boundary of estimation. If
#' using [optimize()], this will be the lower end point of the "interval" argument.
#' @param param_upper Upper boundary of parameter for estimation (numeric).
#' If using the bisection method, this will be the upper boundary of estimation. If
#' using [optimize()], this will be the upper end point of the "interval" argument.
#' Default is 30.
#' @param bisect_root_stop Value for the squared residual that is considered sufficient
#' for stopping the search, if using the bisection method. Any parameter value that
#' produces a squared residual less than root_stop will be considered the root.
#' Default is 1e-4.
#' @param bisect_max_iter Maximum number of iterations to search for the root if
#' using the bisection method. Default is 100.
#' @param bisect_abs_tol Absolute tolerance for stopping the root search if using
#' the bisection method. The search will stop if half the difference between the
#' new lower and upper bounds is less than abs_tol. No warning is given if
#' search is stopped due to abs_tol in the absence of root_stop being achieved.
#' Absolute tolerance is chosen over relative tolerance because the "midpoint"
#' step c may not actually be the middle of a and b, due to jumps made for
#' odeIter() convergence. Default is 1e-8.
#' @param bisect_max_jumps Maximum number of jumps that will be made by the
#' bisection method in order to address non-convergence within iterations
#' of the ordinary differential equations. Default is 10.
#' @param ... Optional arguments that can be provided to [optimize()] if using the
#' optimize method. See documentation for [optimize()].
#'
#' @returns A list of argument values.
#' @export
#'
#' @examples
#'
durationOptControl <- function(
    param_lower = NULL,
    param_upper = 30,
    bisect_root_stop = 1e-4,
    bisect_max_iter = 100,
    bisect_abs_tol = 1e-8,
    bisect_max_jumps = 10,
    ...){

  ## Set up return list
  par_list <- list(
    param_lower = param_lower,
    param_upper = param_upper,
    bisect_root_stop = bisect_root_stop,
    bisect_max_iter = bisect_max_iter,
    bisect_abs_tol = bisect_abs_tol,
    bisect_max_jumps = bisect_max_jumps
  )

  ## ensure inputs are correctly numeric - handle param_lower separately as it can be NULL ##
  num_classes <- unlist(lapply(par_list[!names(par_list) %in% c("param_lower")], is, "numeric"))
  not_num <- names(num_classes[!num_classes])
  if(length(not_num) > 0){
    stop(paste("The following argument(s) need to be numeric:", paste(not_num, collapse = ", ")))
  }

  ## check that param_lower is either NULL or numeric ##
  param_lower_check <- is.null(param_lower) | is(param_lower, "numeric")
  if(!param_lower_check){
    stop("The param_lower argument in durationOptControl() must be either NULL or numeric.")
  }

  ## handle additional optimize arguments. No checks for these, will be passed along
  ## to optimize() as they are
  optimize_args <- list(...) # unpack optimize arguments

  par_list <- c(par_list, optimize_args) # merge lists

  return(par_list)
}

#' Control of sleep midpoint parameter (\eqn{\tau}) optimization
#'
#' Creates a list of arguments that will be used to control the optimization
#' of \eqn{\tau}.
#'
#' @param param_lower Lower boundary of parameter for estimation (numeric).
#' If using the bisection method, this will be the lower boundary of estimation. If
#' using [optimize()], this will be the lower end point of the "interval" argument.
#' Default is 23, based on research into the distribution of the intrinsic
#' circadian period length in humans (see details for reference).
#' @param param_upper Upper boundary of parameter for estimation (numeric).
#' If using the bisection method, this will be the upper boundary of estimation. If
#' using [optimize()], this will be the upper end point of the "interval" argument.
#' Default is 25, based on research into the distribution of the intrinsic
#' circadian period length in humans (see details for reference).
#' @param bisect_root_stop Value for the squared residual that is considered sufficient
#' for stopping the search, if using the bisection method. Any parameter value that
#' produces a squared residual less than root_stop will be considered the root.
#' Default is 1e-4.
#' @param bisect_max_iter Maximum number of iterations to search for the root if
#' using the bisection method. Default is 100.
#' @param bisect_abs_tol Absolute tolerance for stopping the root search if using
#' the bisection method. The search will stop if half the difference between the
#' new lower and upper bounds is less than abs_tol. No warning is given if
#' search is stopped due to abs_tol in the absence of root_stop being achieved.
#' Absolute tolerance is chosen over relative tolerance because the "midpoint"
#' step c may not actually be the middle of a and b, due to jumps made for
#' iterative convergence of the ordinary differential equations. Default is 1e-8.
#' @param bisect_max_jumps Maximum number of jumps that will be made by the
#' bisection method in order to address non-convergence within iterations
#' of the ordinary differential equations. Default is 10.
#' @param ... Optional arguments that can be provided to [optimize()] if using the
#' optimize method. See documentation for [optimize()]. Note that the function (f),
#' interval, and corresponding ... arguments are already provided, so including those
#' (or the lower or upper) arguments will likely cause an error.
#'
#' @returns A list of argument values.
#'
#' @details Reference for default values of param_lower and param_upper.
#'
#' Duffy JF, Cain SW, Chang AM, Phillips AJ, Münch MY, Gronfier C,
#' Wyatt JK, Dijk DJ, Wright KP Jr, Czeisler CA. Sex difference in the
#' near-24-hour intrinsic period of the human circadian timing system.
#' Proc Natl Acad Sci U S A. 2011 Sep 13;108 Suppl 3(Suppl 3):15602-8.
#' doi: 10.1073/pnas.1010666108. Epub 2011 May 2. PMID: 21536890;
#' PMCID: PMC3176605.
#'
#' @export
#'
#' @examples
#'
midpointOptControl <- function(
    param_lower = 23,
    param_upper = 25,
    bisect_root_stop = 1e-4,
    bisect_max_iter = 100,
    bisect_abs_tol = 1e-8,
    bisect_max_jumps = 10,
    ...){

  ## Set up return list
  par_list <- list(
    param_lower = param_lower,
    param_upper = param_upper,
    bisect_root_stop = bisect_root_stop,
    bisect_max_iter = bisect_max_iter,
    bisect_abs_tol = bisect_abs_tol,
    bisect_max_jumps = bisect_max_jumps
  )

  ## ensure inputs are correctly numeric ##
  num_classes <- unlist(lapply(par_list, is, "numeric"))
  not_num <- names(num_classes[!num_classes])
  if(length(not_num) > 0){
    stop(paste("The following argument(s) need to be numeric:", paste(not_num, collapse = ", ")))
  }

  ## handle additional optimize arguments. No checks for these, will be passed along
  ## to optimize() as they are
  optimize_args <- list(...) # unpack optimize arguments

  par_list <- c(par_list, optimize_args) # merge lists

  return(par_list)
}


#' Optimization of \eqn{\mu} and \eqn{\tau} parameters of the Homeostatic-Circadian-Light model
#'
#' Function to estimate \eqn{\mu} and \eqn{\tau} parameters of the Homeostatic-Circadian-Light (HCL)
#' model, as described in the paper by Skeldon et al. (2023). See reference.
#' Ordinary differential equations (ODEs) are handled by the deSolve package.
#' Parameters are optimized sequentially, with \eqn{\mu} optimized first to best
#' match the observed sleep duration, and \eqn{\tau} optimized second to best match
#' the observed sleep midpoint. Optimization options include use of a
#' bisection approach or the [optimize()] function (see details).
#'
#'
#' @param df A data frame representing epoch-level light exposure and sleep/wake
#' states. The data frame should have three columns (in any order): 1) A column representing to cumulative
#' time of each epoch, in 24-hour decimal format (e.g., 2:30 am on day 1 = 2.5,
#' 2:30 am on day 2 = 26.5, etc); 2) a column with binary values representing sleep
#' (1) and wake (0) states (optional if sleep outcomes are provided, see sleep_dur and
#' sleep_mid arguments); and 3) a column with light exposure in lux. Additional
#' columns will not be used. Missing values should be imputed prior to the use
#' of this function. Gaps in time should be fine, although note that deSolve will
#' perform linear interpolation on light. As such, large time gaps may result
#' in poor estimation.
#'
#' @param time_var A string representing the name of the time column in df.
#'
#' @param light_var A string representing the name of the light column in df.
#'
#' @param y0 A named vector representing the initial states of the variables
#' used in the HCL model. The vector must have 5 elements with the following
#' names, in that order: "h", "n", "x", "y", "S". Given that the true starting
#' state is likely to be unknown, particularly given that the HCL model optimizes
#' over different parameters, choice of y0 values is essentially arbitrary. The
#' default (NULL) will use an arbitrary starting point, which will be addressed
#' through the iterative ODE process. It is possible that changes to y0 may
#' speed up ODE convergence, but it is unlikely that any meaningful improvements
#' will be achieved given the number of different parameters that will be tested.
#'
#' @param ode_parms A list of named values for each parameter required by the
#' HCL model. The list must be generated using the [hclParms()] function.
#' See documentation for [hclParms()] for more detail.
#'
#' @param sleep_dur The observed value for sleep duration in hours. If NULL, this
#' value will be calculated from the sleep data in df, specifically as the average
#' sleep per noon-to-noon day.
#'
#' @param sleep_mid The observed value of the sleep midpoint in 24-hour decimal format.
#' If NULL, this value will be calculated from the sleep data in df, specifically
#' as a weighted (by duration) circular average of all sleep periods in a
#' noon-to-noon 24 hour day.
#'
#' @param sleep_var A string representing the name of the sleep/wake column in df.
#' This only needs to be provided if either sleep_dur or sleep_mid are NULL, as
#' it will be used to calculate the missing value using the observed data.
#'
#' @param max_ode_iter The maximum number of iterations permitted for each run of
#' the ODE models to establish convergence. The default is 20. Increasing this number
#' may help some cases where ODE models are not converging. However, models that
#' do not converge between 20-40 iterations probably won't be helped by further iterations.
#' For example, there may be insufficient light exposure to entrain at the specified
#' \eqn{\tau}. The function will explore a range of parameter values to identify
#' those that lead to convergence, which will then be compared against the observed
#' sleep outcomes.
#'
#' @param dur_tol Tolerance allowed for sleep duration (hours) to determine convergence
#' during ODE iteration. The default is 1/60 (i.e., 1 minute), meaning that
#' the average sleep duration for successive iterations of the ODE model must not differ
#' by more than a minute. Increasing this number will make ODE convergence easier
#' to obtain.
#'
#' @param mid_tol Tolerance allowed for sleep midpoint (hours) to determine convergence
#' during ODE iteration. The default is 1/60 (i.e., 1 minute), meaning that
#' the average sleep midpoint for successive iterations of the ODE model must not differ
#' by more than a minute. Increasing this number will make ODE convergence easier
#' to obtain. Given that ODE non-convergence is typically due to sleep midpoint
#' changes over iterations, relaxing this can be a good initial step to troubleshooting.
#'
#' @param compiled Boolean. If TRUE (default), deSolve will be called using complied C code
#' instead of R code, which is much faster. C and R code returns identical results,
#' so leaving this argument as TRUE is recommended.
#'
#' @param opt_method A string representing the desired method for optimizing
#' \eqn{\mu} and \eqn{\tau} parameters. Options are "bisect" or "optimize".
#' "bisect" will use a bisection approach that includes additional steps for
#' addressing non-convergence of the ODE iterations. "optimize" (default) will use the
#' [optimize()] function (see details for caveat). Based on extremely limited testing,
#' "optimize" appears slightly faster, although "bisect" is also fast provided it
#' does not need to spend excessive time jumping around at the lower or upper
#' boundaries to establish ODE convergence when estimating \eqn{\tau}.
#' This means more extreme values for param_lower or param_upper when
#' optimizing on sleep midpoint will slow down the bisection approach.
#'
#' @param duration_opt_control A list of named values for the control of \eqn{\mu}
#' optimization. Values must be provided using the [durationOptControl()] function.
#' See [durationOptControl()] function documentation for more details.
#'
#' @param midpoint_opt_control A list of named values for the control of \eqn{\tau}
#' optimization. Values must be provided using the [midpointOptControl()] function.
#' See [midpointOptControl()] function documentation for more details.
#'
#' @returns A list with the following elements:
#' \itemize{
#'  \item{"opt_results"}{A data frame with values and squared residuals for \eqn{\mu} and \eqn{\tau}.}
#'  \item{"opt_convergence_status"}{A number representing the convergence status
#'  of parameter optimization. 1 = converged, 0 = not converged.}
#'  \item{"opt_convergence_message"}{A message corresponding to the optimization convergence status.}
#'  \item{"ode_df"}{A data frame with the solved states of the ODE equations at each timepoint.}
#'  \item{"ode_sleep_sum"}{A data frame with the summaries of sleep for each ODE iteration.}
#'  \item{"ode_converge"}{A boolean indicating if the final ODE model converged (TRUE) or not (FALSE)}
#'  \item{"ode_converge_message"}{A message corresponding to the convergence status of the final ODE model.}
#'  \item{"ode_converge_df"}{A data frame with additional information on the convergence of the final ODE model.}
#'  \item{"ode_iterations"}{The number of iterations the final ODE model took to converge. If
#'  convergence did not occur, this will be the maximum number of iterations allowed when calling the function.}
#' }
#'
#' @reference Skeldon AC, Rodriguez Garcia T, Cleator SF, Della Monica C,
#' Ravindran KKG, Revell VL, Dijk DJ. Method to determine whether sleep
#' phenotypes are driven by endogenous circadian rhythms or environmental light
#' by combining longitudinal data and personalised mathematical models. PLoS
#' Comput Biol. 2023 Dec 22;19(12):e1011743. doi: 10.1371/journal.pcbi.1011743.
#' PMID: 38134229; PMCID: PMC10817199.
#'
#' @details Note that [optimize()] does not allow the return of NAs in the objective function.
#' However, the ODE models may not converge for certain parameter values,
#' particular values of \eqn{\tau} that are further away from 24. As such,
#' residuals are not available in these cases. This is handled by returning an
#' arbitrarily high value to [optimize()] at the point of ODE non-convergence.
#' This seems to work fine in practice, but it is unclear if this may cause
#' estimation problems in certain cases.
#'
#' @export
#'
#' @examples
#'
rhcl <- function(
    df,
    time_var,
    light_var,
    y0 = NULL,
    ode_parms = hclParms(),
    sleep_dur = NULL,
    sleep_mid = NULL,
    sleep_var = NULL,
    max_ode_iter = 20,
    dur_tol = 1/60,
    mid_tol = 1/60,
    compiled = TRUE,
    opt_method = c("bisect", "optimize"),
    duration_opt_control = durationOptControl(),
    midpoint_opt_control = midpointOptControl()
    ){

  ### TODO - build in checks ###
  ## TODO - consider a function that does some preprocessing of the data.frame
  # such as accepting POSIXct and converting, as well as running all the checks

  ### Pre-process data.frame and check data ###
  df <- dfPrep(df = df, time_var = time_var, light_var = light_var, sleep_var = sleep_var)

  ## Require sleep_var argument if either sleep_dur or sleep_mid is NULL ##
  if(is.null(sleep_dur) | is.null(sleep_mid)){
    if(is.null(sleep_var)){
      stop(paste("Name of sleep variable in data frame must be provided",
                 "as 'sleep_var' argument if either 'sleep_dur' or 'sleep_mid'",
                 "are NULL."))
    } else{
      ## Calculated observed sleep values if not provided ##
      sleep_sum <- sleepSummary(df=df, sleep_var = sleep_var, time_var = "ctime")

      # sleep duration if needed
      if(is.null(sleep_dur)){
        sleep_dur <- sleep_sum$summar_sleep_dur_noon_24hr
      }

      # sleep midpoint if needed
      if(is.null(sleep_mid)){
        sleep_mid <- sleep_sum$summary(sleep_mid)
      }
    }
  }

  ## check that y0 has appropriate length and names ##
  if(!is.null(y0)){
    if(length(y0)!=5 | names(y0) != c("h", "n", "x", "y", "S")){
      stop("y0 must be a named vector with 5 values and names = c('h', 'n', 'x', 'y', 'S')")
    }
  } else{
    # TODO - arbitrary values for now. Could possibly speed up convergence by basing values on
    # intial time of data
    y0 <- c(h = 13.15, n = .152, x = -0.966, y = -0.558, S = 0)
  }


  ## if opt_method left as default, use optimize (faster based on limited testing)
  if(identical(opt_method, c("bisect", "optimize"))){
    opt_method <- "optimize"
  }

  ### wrap in tryCatch() to ensure interpolation function for R code, which
  # is added to the global environment, is always removed even if an error occurs
  # There has to be a way for deSolve to access the interpolation function without it
  # being in the global environment, but I haven't figured out how.
  tryCatch(
    {
      ### Set up input for deSolve::ode() ###
      if(compiled){
        ## desolve list for compiled code ##
        desolve_list <- list(
          y = y0, # initial values
          times = df[["ctime"]], # times vector
          func = "derivsc_p", # c function to call for derivative equations
          parms = unlist(ode_parms), #parameters
          dllname = "rHCL", # c library for package
          initforc = "forcc_p", # c function for forcing variable initialization
          forcings = cbind(df[["ctime"]], df[[light_var]]), # matrix of forcing variables
          fcontrol = list(method = "linear", rule=2, f=0), # forcing control arguments
          initfunc = "parmsc_p", # c function for initializing parameters for deSolve
          nout = 0, # number of additional variables for deSolve to return
          events = list(func = "eventc_p", root = TRUE), # arguments for events
          rootfun = "rootc_p", # c function for roots
          nroot = 1 # number of roots for deSolve to track
        )

      } else{
        ## desolve list for R code ##
        desolve_list <- list(
          y = y0, # initial values
          func = dHCL, # R derivative function
          times = df[["ctime"]], # times vector
          parms = ode_parms, #parameters
          events = list(func = dEventFunc, root = TRUE),
          rootfun = dRootFunc
        )

        ## Create light interpolation function for R code ##
        assign("light.int",
               approxfun(x=df[["ctime"]], y=df[[light_var]], method="linear", rule=2),
               envir = .GlobalEnv) # create interpolation function; Note this is creating a function in the global environment and needs to be cleaned up
      }


      ### optimize sleep duration first ###
      ## set tau_c to 24.2 for convergence purposes ##
      desolve_list[["parms"]][["tau_c"]] <- 24.2

      ## check min mu compared to provided value ##
      min_mu <- desolve_list[["parms"]][["Hzero"]] + desolve_list[["parms"]][["ca_par"]] + desolve_list[["parms"]][["delta"]]*.05 # constrain mu to be greater than this

      # if lower bound for mu is NULL, use min_mu
      if(is.null(duration_opt_control[["param_lower"]])){
        duration_opt_control[["param_lower"]] <- min_mu
      }
      # else if lower bound for mu is too low, use min_mu
      else if(duration_opt_control[["param_lower"]] < min_mu){
        warning(paste("Minimum mu value provided in dur_control argument (param_lower) is",
                      "below the minimum allowed value. Replaced with", round(min_mu,0)))
        duration_opt_control[["param_lower"]] <- min_mu
      }


      ## bisect or optimize methods ##
      if(opt_method == "bisect"){
        opt_duration <- odeBisect(
          param_lower = duration_opt_control[["param_lower"]],
          param_upper = duration_opt_control[["param_upper"]],
          observed_param = sleep_dur,
          root_stop = duration_opt_control[["bisect_root_stop"]],
          max_iter = duration_opt_control[["bisect_max_iter"]],
          abs_tol = duration_opt_control[["bisect_abs_tol"]],
          method = "mu",
          num_ode_jumps = duration_opt_control[["bisect_max_jumps"]],
          desolve_args = desolve_list,
          max_ode_iter = max_ode_iter,
          dur_tol = dur_tol,
          mid_tol = mid_tol
        )

      } else if(opt_method == "optimize"){
        ## extract additional arguments being passed to optimize() ##
        optimize_args_dur <- duration_opt_control[!names(duration_opt_control) %in% c(
          "param_lower", "param_upper", "bisect_root_stop", "bisect_max_iter",
          "bisect_abs_tol", "bisect_max_jumps")] # remove additional arguments

        # prepare arguments for do.call
        optimize_args_dur <- c(list("f" = odeOptim_duration,
                                "interval" = c(duration_opt_control[["param_lower"]], duration_opt_control[["param_upper"]]),
                                "sleep_dur" = sleep_dur,
                                "desolve_args" = desolve_list,
                                "max_iter" = max_ode_iter,
                                "dur_tol" = dur_tol,
                                "mid_tol" = mid_tol),
                           optimize_args_dur)

        opt_duration <- do.call("optimize", optimize_args_dur) # optimize mu


      }

      ### optimize sleep midpoint ###
      ## switch mu to estimated value ##
      desolve_list[["parms"]][["mu"]] <- opt_duration$minimum

      ## bisect or optimize methods ##
      if(opt_method == "bisect"){
        opt_midpoint <- odeBisect(
          param_lower = midpoint_opt_control[["param_lower"]],
          param_upper = midpoint_opt_control[["param_upper"]],
          observed_param = sleep_mid,
          root_stop = midpoint_opt_control[["bisect_root_stop"]],
          max_iter = midpoint_opt_control[["bisect_max_iter"]],
          abs_tol = midpoint_opt_control[["bisect_abs_tol"]],
          method = "tau_c",
          num_ode_jumps = midpoint_opt_control[["bisect_max_jumps"]],
          desolve_args = desolve_list,
          max_ode_iter = max_ode_iter,
          dur_tol = dur_tol,
          mid_tol = mid_tol
        )

        final_res <- opt_midpoint[["ode_res"]] # extract final ODE results

      } else if(opt_method == "optimize"){
        ## extract additional arguments being passed to optimize() ##
        optimize_args_mid <- midpoint_opt_control[!names(midpoint_opt_control) %in% c(
          "param_lower", "param_upper", "bisect_root_stop", "bisect_max_iter",
          "bisect_abs_tol", "bisect_max_jumps")] # remove additional arguments

        # prepare arguments for do.call
        optimize_args_mid <- c(list("f" = odeOptim_midpoint,
                                    "interval" = c(midpoint_opt_control[["param_lower"]], midpoint_opt_control[["param_upper"]]),
                                    "sleep_mid" = sleep_mid,
                                    "desolve_args" = desolve_list,
                                    "max_iter" = max_ode_iter,
                                    "dur_tol" = dur_tol,
                                    "mid_tol" = mid_tol),
                               optimize_args_mid)

        opt_midpoint <- do.call("optimize", optimize_args_mid) # optimize tau_c

        ## obtain solved ODE, as optimize does not return it like the bisection method does
        desolve_list[["parms"]][["tau_c"]] <- opt_midpoint$minimum
        final_res <- odeIter(desolve_args=desolve_list, max_iter = max_ode_iter, dur_tol = dur_tol, mid_tol = mid_tol)
      }


    }, finally = {
      ## Clean up the function I added to the global environment ##
      if(exists("light.int", where = .GlobalEnv)){
        rm(light.int, envir = .GlobalEnv)
      } # remove light approxfun
    }
  )

  ### Check outcome convergence ###
  duration_converge <- opt_duration$objective < .03 # residual (already squared) < .03
  midpoint_converge <- opt_midpoint$objective < .03 # residual (already squared) < .03

  if(duration_converge & midpoint_converge){
    conv_status <- 1 # convergence status
    conv_message <- "Convergence (squared residual < .03) obtained for both sleep duration and midpoint."
  } else{
    conv_status <- 0
    conv_message <- "Convergence (squared residual < .03)  not obtained for both sleep duration and midpoint. Check squared residuals to diagnose."
  }

  ### prepare other results ###
  res <- list(opt_results =
                ## Results of optimization ##
                data.frame(parameter = c("mu", "tau_c"),
                           value = c(opt_duration$minimum, opt_midpoint$minimum),
                           resid_squared = c(opt_duration$objective, opt_midpoint$objective)),
              ## Status of optimization convergence ##
              opt_convergence_status = conv_status,
              ## Optimization message ##
              opt_convegence_message = conv_message,
              ## Results of ODEs using final estimated parameters ##
              ode_df = final_res$ode_res,
              ## Summary of sleep per iteration for ODE results ##
              ode_sleep_sum = final_res$sleep_sum,
              ## Convergence status for ODE run using final parameters ##
              ode_converge = final_res$converge,
              ## Convergence message for final ODE run ##
              ode_converge_message = final_res$conv_message,
              ## Dataframe showing results of ODE convergence over iterations ##
              ode_converge_df = final_res$converge_df,
              ## Number of iterations in final ODE ##
              ode_iterations = final_res$iterations
              )

  ### trigger warnings if convergence for either optimization or final ODE run is not obtained ###
  if(!res$opt_convergence_status){
    warning(res$opt_convergence_message)
  }

  if(!res$ode_converge){
    warning(res$ode_converge_message)
  }

  ### return results ###
  return(res)
}
