# Create a new test file at tests/testthat/test-par-compartment-1.R
library(testthat)
library(PKstat)

test_that("par_compartment_1 works correctly", {
  data <- data.frame(time = c(0, 1, 2, 3), conc = c(10, 8, 6, 4))
  result <- par_1c(data, time = "time", conc = "conc")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("parameter", "value"))
})
