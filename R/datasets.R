#' Example light and sleep data
#'
#' Synthetic data set of light exposure and sleep over 14 days measured at
#' 1-minute epochs.
#'
#' @format ## `rhcl_df`
#' A data farme with 20160 rows and 3 columns:
#' \describe{
#'  \item{times}{Date and time values in POSIXct format}
#'  \item{light}{Values of light in lux}
#'  \item{sleep}{Binary variable: 0 = wake and 1 = sleep}
#' }
#' @source Synthetic data
"rhcl_df"
