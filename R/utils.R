minFinder <- function(y, time, ..., smooth = FALSE){
  # find the local minimums of a vector

  if(smooth){
    # create spline function
    s_func <- splinefun(x = y, y = time, method = "natural")

  }


}
