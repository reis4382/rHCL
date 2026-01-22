#' Differential equations for HCL model to be used with [deSolve::ode()]
#'
#' This implements the equations described by Skeldon et al. 2023 for
#' the HCL model. An edit is made for the derivative of the photoreceptor activation
#' to be consistent with Forger 1999. Specifically, time scaling is explicitly
#' placed into the equation to allow non-hour time scales. Requires a global
#' interpolation function (named "light.int") that provides the interpolated light
#' value at time t.
#'
#'
#' @param time A vector of time (must be properly scaled with time_scale parameter)
#' @param states A vector of named states and corresponding values for current state of system
#' @param parms A vector of named parameters for use in the equations
#'
#' @returns A list of the derivatives for each state.
#' @noRd
#' @keywords internal
#'
dHCL <- function(time, states, parms){
  # construct environment for accessing states and parms
  with(as.list(c(states, parms)), {
    Itilde = light.int(time) # note: light.int must be a function from approxfun() specified in the environment.

    ### Auxiliary values ###
    Itilde = (1 - S) * Itilde # Equation 5: set light to 0 if asleep
    beta_hat = G_par * alpha_zero * (Itilde / Izero)^p_par * (1 - n) # Equation 7
    B_par = (1 - little_b_par*x) * (1 - little_b_par*y) * beta_hat # Equation 10

    ### Derivatives ###
    ## sleep homeostatic pressure derivative ##
    dhdt = (-h + (1 - S)*mu) / chi # Equation 2

    ## photoreceptor derivative ##
    # Note addition of leading "60*", which is present in Forger 1999.
    # Skeldon 2023 appear to drop this and incorporate it into default parameters
    # (alpha_zero and beta), specifically scaled for time in seconds. However,
    # that does not allow for flexibility in time input. I also think it the
    # alterations of the default parameters may introduce an error in the calculation
    # of beta_hat. Modifying formula here
    dndt = (60*(alpha_zero * (Itilde / Izero)^p_par * (1 - n) - beta*n)) / time_scale # Equation 6

    ## derivative of x (I believe this is xc in forger 1999) ##
    dxdt = (gamma*(x - (4*x^3/3)) - y*((24 / (f_par * tau_c))^2 + k_par * B_par)) / (12/pi * time_scale) # Equation 8

    ## derivative of y (I believe this is x in forger 1999) ##
    dydt = (x + B_par) / (12/pi * time_scale) # Equation 9

    ## "derivative" for sleep (always 0 b/c it doesn't change dynamically, only during root function) ##
    # This is used to help the root switching functions #
    dsleepdt = 0

    ### Return derivatives - must match order of state variables ###
    return(list(c(dhdt = dhdt, dndt = dndt, dxdt = dxdt, dydt = dydt, dsleepdt = dsleepdt)))
  })
}

#' Function to calculate circadian value for homeostatic sleep switching
#'
#' @param x Circadian pacemaker auxiliary value at a given time
#' @param y Circadian pacemaker primary value at a given time
#'
#' @returns A numeric value indicating the circadian wake propensity at given time
#' @noRd
#' @keywords internal
#'
circFunction <- function(x, y){
  # Scalars taken from appendix of skeldon 2023 paper
  alpha_vals <- c(0.7896, -0.3912, 0.7583, -0.4442, 0.0250, -0.9647)

  # calculate circadian value
  circ_val <- (alpha_vals %*% c(1, x, y, x^2, x*y, y^2))[[1]] # matrix multiplication

  return(circ_val)
}

#' Function to identify roots during deSolve ODE calculations
#'
#' Roots are when the homeostatic sleep pressure crosses the appropriate
#' threshold.
#'
#' @param time Current time of the ODE equations
#' @param states Vector with named values representing current state of ODE system.
#' @param parms Parameter list used for ODE functions
#'
#' @returns Boolean if root is found at current time step.
#' @noRd
#' @keywords internal
#'
dRootFunc <- function(time, states, parms){
  # attach parameter and states
  with(as.list(c(states, parms)), {

    # calculate current circadian wake propensity
    circ_prop <- circFunction(x=x, y=y)

    # determine if appropriate threshold is crossed based on current sleep/wake
    if(S == 0){
      h_thresh = Hzero + 0.5 * delta + ca_par * circ_prop; # eq. 3; threshold for sleep if awake
    } else if (S == 1){
      h_thresh = Hzero - 0.5 * delta + ca_par * circ_prop; # eq. 4; threshold for wake if asleep
    }

    return(h - h_thresh) # triggers when difference equals 0

  })
}

#' Function to identify roots during deSolve ODE calculations if enforcing wake periods
#'
#' Roots are when the homeostatic sleep pressure crosses the appropriate
#' threshold. Requires an additional global interpolation function (named "force.wake") that carries
#' forward any enforced wake forcing variable.
#'
#' @param time Current time of the ODE equations
#' @param states Vector with named values representing current state of ODE system.
#' @param parms Parameter list used for ODE functions
#'
#' @returns Boolean if root is found at current time step.
#' @noRd
#' @keywords internal
#'
dRootFunc_FW <- function(time, states, parms){
  # attach parameter and states
  with(as.list(c(states, parms)), {

    # calculate current circadian wake propensity
    circ_prop <- circFunction(x=x, y=y)

    ## enforce wake
    wake_thresh = 0
    if(force.wake(time)==1){
      wake_thresh = 100
    }

    # determine if appropriate threshold is crossed based on current sleep/wake
    if(S == 0){
      h_thresh = Hzero + 0.5 * delta + ca_par * circ_prop + wake_thresh; # eq. 3; threshold for sleep if awake
    } else if (S == 1){
      h_thresh = Hzero - 0.5 * delta + ca_par * circ_prop + wake_thresh; # eq. 4; threshold for wake if asleep
    }

    return(h - h_thresh) # triggers when difference equals 0

  })
}

#' Event function to switch sleep state when a root is identified by ode
#'
#' @param time Current time of the ODE equations
#' @param states Vector with named values representing current state of ODE system.
#' @param parms Parameter list used for ODE functions
#'
#' @returns A named vector of updated state variables (only sleep state will change).
#' @noRd
#' @keywords internal
#'
dEventFunc <- function(time, states, parms){
  # attach states and parameters
  with(as.list(c(states, parms)), {

    # reverse sleep state if root is detected
    states["S"] <- c(1,0)[S+1] # index will equal 1 when S=0 and 2 when S=1

    # return updated states
    return(states)
  })

}

#' Function establishing default parameters for ODE system
#'
#' This function sets up the parameters used as input to ordinary differential equations (ODEs).
#' Values are taken from default values of Skeldon 2023 paper (see references). However, modifications
#' have been made to the option for setting the time scale (\eqn{\kappa} parameter
#' in paper), which also affects the scaling of \eqn{\chi}. Additionally,
#' \eqn{\alpha<sub>0} and \eqn{\beta} have been tweaked to remain consistent
#' with the modifications made to the time_scale parameter and edits to the
#' photoreceptor derivative equation that now follows the Forger 1999 equation (see references).
#'
#' New parameter values can be specified by setting their respective arguments.
#' Unused arguments will trigger an error.
#'
#' @param mu Upper asymptote for sleep pressure
#' @param chi Time constant for sleep pressure decay/rise. Default units assumes time values are in hours.
#' @param Hzero Mean level of wake propensity rhythm
#' @param delta Separation between thresholds
#' @param ca_par Circadian wake propensity amplitude
#' @param tau_c Circadian period
#' @param f_par Correction factor for oscillator
#' @param G_par Gain factor determining magnitude of light effect
#' @param p_par Light sensitivity
#' @param k_par Determines relative effect of light on oscillator variables
#' @param little_b_par Sensitivity modulation factor
#' @param gamma Stiffness of van der Pol oscillator
#' @param alpha_zero Magnitude of effect of light on fraction of activated photoreceptors
#' @param beta Decay rate of fraction of activated photoreceptors
#' @param Izero Scaling factor for light
#' @param time_scale Scale of the time variable in hours. Three character values are allowed:
#' "hours" (default), "mins", and "secs", which indicate the time variable is
#' scaled to hours, minutes, or seconds of the day respectively. Alternatively, a
#' numeric value can be entered, representing how many hours of the day the time
#' variable represents. For example, enter 1/60 if time is scaled in minutes,
#' or 1/(60*60) if time is scaled to seconds.
#'
#' @returns A named list of all parameters required by the ODE functions.
#'
#' @references Skeldon AC, Rodriguez Garcia T, Cleator SF, Della Monica C,
#' Ravindran KKG, Revell VL, Dijk DJ. Method to determine whether sleep phenotypes
#' are driven by endogenous circadian rhythms or environmental light by combining
#' longitudinal data and personalised mathematical models. PLoS Comput Biol. 2023
#' Dec 22;19(12):e1011743. doi: 10.1371/journal.pcbi.1011743. PMID: 38134229;
#' PMCID: PMC10817199.
#'
#' Forger DB, Jewett ME, Kronauer RE. A simpler model of the human circadian
#' pacemaker. J Biol Rhythms. 1999 Dec;14(6):532-7.
#' doi: 10.1177/074873099129000867. PMID: 10643750.
#'
#' @export
#'
#' @examples
hclParms <- function(mu = 17.87,
                     chi = 45,
                     Hzero = 13,
                     delta = 1,
                     ca_par = 1.72,
                     tau_c = 24.2,
                     f_par = 0.99669,
                     G_par = 19.9,
                     p_par = 0.6,
                     k_par = 0.55,
                     little_b_par = 0.4,
                     gamma = 0.23,
                     alpha_zero = 0.16,
                     beta = 0.013,
                     Izero = 9500,
                     time_scale = "hours"){

  ## Set up return list
  par_list <- list(
    mu = mu,
    chi = chi,
    Hzero = Hzero,
    delta = delta,
    ca_par = ca_par,
    tau_c = tau_c,
    f_par = f_par,
    G_par = G_par,
    p_par = p_par,
    k_par = k_par,
    little_b_par = little_b_par,
    gamma = gamma,
    alpha_zero = alpha_zero,
    beta = beta,
    Izero = Izero,
    time_scale = NA # time scale constant
  )

  ## ensure inputs are correctly numeric ##
  num_classes <- unlist(lapply(par_list[!names(par_list) %in% "time_scale"], is, "numeric"))
  not_num <- names(num_classes[!num_classes])
  if(length(not_num) > 0){
    stop(paste("The following arguments need to be numeric:", paste(not_num, collapse = ", ")))
  }

  ## extract time scale adjustment ##
  if(time_scale == "hours"){
    time_scale <- 1
  } else if(time_scale == "mins"){
    time_scale = 1/60
  } else if(time_scale == "secs"){
    time_scale = 1/60/60
  } else if(is(time_scale, "numeric")){
    time_scale = time_scale
  } else{
    stop("time_scale argument must be either a numeric value or one of the following: 'hours', 'mins', 'secs'")
  }

  ## Adjust chi and kappa for time scale ##
  par_list[["chi"]] <- par_list[["chi"]] / (time_scale)
  par_list[["time_scale"]] <- 1 / (time_scale)

  return(par_list)
}



