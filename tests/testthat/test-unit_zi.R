library(saeczi)
data(pop)
data(samp)

lin_formula <- DRYBIO_AG_TPA_live_ADJ ~ tcc16 + elev

set.seed(5)
result <- unit_zi(samp,
                  pop, 
                  lin_formula,
                  domain_level = "COUNTYFIPS",
                  mse_est = TRUE,
                  B = 5,
                  parallel = FALSE)

test_that("printed result is as expected", {
  expect_snapshot(result)
})

test_that("result is as expected", {
  expect_equal(result$res$mse[2], 80.5, tolerance = 0.01)
  expect_equal(result$res$mse[12], 144.4, tolerance = 0.01)
  expect_equal(result$res$mse[36], 116.6, tolerance = 0.01)
  
  expect_equal(result$res$est[8], 104.4, tolerance = 0.01)
  expect_equal(result$res$est[17], 66.7, tolerance = 0.01)
})

test_that("result[[2]] is a df", {
  expect_s3_class(result[[2]], "data.frame")
})

test_that("mse column is numeric", {
  expect_type(result[[2]][[2]], "double")
})

test_that("est column is numeric", {
  expect_type(result[[2]][[3]], "double")
})

test_that("number of rows in final data frame matches domain levels", {
  expect_equal(nrow(result[[2]]), length(unique(pop$COUNTYFIPS)))
})

test_that("check for error when pop_dat is a vector or list", {
  expect_error(unit_zi(samp_dat = samp$tcc16, pop_dat = pop$DRYBIO_AG_TPA_live_ADJ_TONS, 
                       lin_formula = DRYBIO_AG_TPA_live_ADJ_TONS ~ tcc16, parallel = FALSE))
})
