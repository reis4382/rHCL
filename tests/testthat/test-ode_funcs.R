# Test differential equations ---------------------------------------------
test_that("dHCL() function works", {

  ## prepare input ##
  times <- 0:1439 / 24
  assign("light.int", approxfun(x=times, y=0, method="linear", rule=2)) # create interpolation function
  states <- c(h = 13, n = 0, y = 0, x = 1)
  ## TODO - store parameters in an rdata object?


})
