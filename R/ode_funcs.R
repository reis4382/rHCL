#' Differential equations for HCL model to be used with deSolve ode()
#'
#' @details
#'This implements the equations described by Skeldon et al. 2023 for the HCL model.
#'
#'
#' @param time A vector of time (must be properly scaled with \eqn{\kappa} parameter)
#' @param states A vector of named states and corresponding values for current state of system
#' @param parms A vector of named parameters for use in the equations
#'
#' @returns A list of the derivatives for each state.
#' @noMd
#' @keywords internal
#'
dHCL <- function(time, states, parms){
  # construct environment for accessing states and parms
  with(as.list(c(states, parms)), {
    Itilde = light.int(time) # note: light.int must be a function from approxfun() specified in the environment.

    ### Auxiliary values ###
    Itilde = (1 - S) * Itilde # Equation 5: set light to 0 if asleep
    beta_hat = G_par * alpha_zero * (Itilde / I_zero)^p * (1 - n) # Equation 7
    B_par = (1 - little_b_par*x) * (1 - little_b_par*y) * beta_hat # Equation 10

    ### Derivatives ###
    ## sleep homeostatic pressure derivative ##
    dhdt = (-h - (1 - S)*mu) / chi # Equation 2

    ## photoreceptor derivative ##
    dndt = alpha_zero * (Itilde / Izero)^p * (1 - n) - beta*n # Equation 6

    ## derivative of x (xc in other models?) ##
    dxdt = (gamma_par*(x - (4*x^3/3)) - y*((24 / f_par * tau_c)^2 + k_par * B_par)) / kappa_par # Equation 8

    ## derivative of y (x in other models?) ##
    dydt = (x + B_par) / (kappa_par) # Equation 9

    ## "derivative" for sleep (always 0 b/c it doesn't change dynamically) ##
    # This is used to help the root switching functions #
    dsleepdt = 0

    ### Return derivatives - must match order of state variables ###
    return(list(c(dhdt = dhdt, dndt = dndt, dxdt = dxdt, dydt = dydt, dsleepdt = dsleepdt)))
  })
}

### TESTING ADDING SOME CODE ###
testFunc <- function(){
  NULL
}
