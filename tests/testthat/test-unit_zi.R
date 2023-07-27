library(saezi)
library(testthat)

data(pop)
data(samp)


lin_formula <- DRYBIO_AG_TPA_live_ADJ ~ tcc16 + tmean + tri

result <- unit_zi(samp, pop, lin_formula, domain_level = "COUNTYFIPS", mse_est = TRUE, parallel = FALSE)

test_that("result is a list", {
  expect_type(result, "list")
})

test_that("result[[1]] is a df", {
  expect_s3_class(result[[1]], "data.frame")
})

test_that("mse column is numeric", {
  expect_type(result[[1]][[2]], "double")
})

test_that("est column is numeric", {
  expect_type(result[[1]][[3]], "double")
})

test_that("number of rows in final data frame matches domain levels", {
  expect_equal(nrow(result[[1]]), length(unique(pop$COUNTYFIPS)))
})

test_that("check for error when pop_dat is a vector or list", {
  expect_error(unit_zi(samp_dat = samp$tcc16, pop_dat = pop$DRYBIO_AG_TPA_live_ADJ_TONS, 
                       lin_formula = DRYBIO_AG_TPA_live_ADJ_TONS ~ tcc16, parallel = FALSE))
})
