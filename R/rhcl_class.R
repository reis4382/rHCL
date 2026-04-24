#' Create class for displaying HCL model results
#'
#' @param opt_results A data.frame with optimization results.
#' @param opt_convergence_status 0 = not converged, 1 = converged
#' @param opt_convergence_message Message corresponding to convergence status.
#' @param ode_df Results of final ODE model from [ode_iter()].
#' @param ode_sleep_sum Results of final sleep summary from [ode_iter()].
#' @param ode_converge TRUE if final ODE model converged.
#' @param ode_converge_message Message corresonding to ode_converge.
#' @param ode_converge_df Convergence data.frame for ODE model.
#' @param ode_iterations Number of iterations it took to achieve convergence.
#' Will be the maximum number allowed if not converged.
#'
#' @returns A list with the rhcl_mod class.
#' @noRd
#'
new_rhcl_mod <- function(
    opt_results = NULL,
    opt_convergence_status = NULL,
    opt_convergence_message = NULL,
    ode_df = NULL,
    ode_sleep_sum = NULL,
    ode_convergence_status = NULL,
    ode_convergence_message = NULL,
    ode_converge_df = NULL,
    ode_iterations = NULL,
    epoch_length_min = NULL,
    min_observed_hours = NULL
    ){

  # TODO - could add checks here, but this is only for internal use #

  # combine elements into a list
  x <- list(
    opt_results = opt_results,
    opt_convergence_status = opt_convergence_status,
    opt_convergence_message = opt_convergence_message,
    ode_df = ode_df,
    ode_sleep_sum = ode_sleep_sum,
    ode_convergence_status = ode_convergence_status,
    ode_convergence_message = ode_convergence_message,
    ode_converge_df = ode_converge_df,
    ode_iterations = ode_iterations,
    epoch_length_min = epoch_length_min,
    min_observed_hours = min_observed_hours
  )

  # create the new results class
  structure(
    x,
    class = "rhcl_mod"
    )
}


#' @returns \code{NULL}
#' @method print rhcl_mod
#' @export
print.rhcl_mod <- function(x, ...){
  # prep numeric table #
  opts <- x$opt_results
  opts[,2:3] <- apply(opts[,2:3], 2, function(y){sprintf("%.2f", y)})

  # cat can't handle df
  cat("Optimization results:")
  cat("\n")
  cat("\n")
  print(opts)

  part2 <- paste("Optimization convergence status:", x$opt_convergence_message)
  part3 <- paste("ODE convergence status:", x$ode_convergence_message)

  cat("\n", "\n", part2, "\n", "\n", part3, "\n")

  invisible(x)
}

#' @returns \code{NULL}
#' @method summary rhcl_mod
#' @export
summary.rhcl_mod <- function(object, ...){

  res <- list(
    opt = object$opt_results,
    opt_convergence_status = object$opt_convergence_status,
    opt_convergence_message = object$opt_convergence_message,

    ode_convergence_status = object$ode_convergence_status,
    ode_convergence_message = object$ode_convergence_message,
    ode_iterations = object$ode_iterations,

    ode_sleep = sleepSummary(df=object$ode_df, sleep_var = "S", time_var = "dtime",
                             epoch_length_min = object$epoch_length_min,
                             min_observed_hours = object$min_observed_hours)$summary

  )



  class(res) <- "summary_rhcl_mod" # create summary class

  return(res)
}

#' @returns \code{NULL}
#' @method print summary_rhcl_mod
#' @export
print.summary_rhcl_mod <- function(x, ...){

  # prepare convergences #
  opt_converge <- c("FALSE", "TRUE")[x[["opt_convergence_status"]]+1]
  ode_converge <- c("FALSE", "TRUE")[x[["ode_convergence_status"]]+1]

  # prepare opt results #
  # prep numeric table #
  opts <- x$opt
  opts[,2:3] <- apply(opts[,2:3], 2, function(y){sprintf("%.2f", y)})

  # print output #
  cat("OPTIMIZATION RESULTS")
  cat("\n", "\n")
  cat("    Optimization convergence:", opt_converge)
  cat("\n", "\n")
  print(opts) # cat can't handle df
  cat("\n")
  cat(x[["opt_convergence_message"]])
  cat("\n", "\n")

  cat("--------------------------------------")
  cat("\n", "\n")

  cat("ODE MODEL RESULTS")
  cat("\n", "\n")
  cat("    ODE convergence", ode_converge)
  cat("\n")
  cat("    ODE iterations", x[["ode_iterations"]])
  cat("\n", "\n")
  cat("ODE Sleep Summary")
  cat("\n")
  print(x[["ode_sleep"]])
  cat("\n")
  cat(x[["ode_convergence_message"]])
  cat("\n")

  invisible(x)
}
