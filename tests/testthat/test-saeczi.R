data(pop)
data(samp)

set.seed(5)
suppressWarnings(
  result <- saeczi(samp,
                   pop, 
                   lin_formula = DRYBIO_AG_TPA_live_ADJ ~ tcc16,
                   domain_level = "COUNTYFIPS",
                   mse_est = TRUE,
                   B = 10L,
                   parallel = FALSE)
)

test_that("result$res is a df", {
  expect_s3_class(result$res, "data.frame")
})

test_that("mse column exists and is not NA", {
  expect_contains(names(result$res), "mse")
  expect_equal(all(!is.na(result$res$mse)), T)
})

test_that("correct number of rows in result data.frame", {
  expect_equal(nrow(result$res), length(unique(pop$COUNTYFIPS)))
})

test_that("transform function checks work", {
  
  # no inv function supplied
  expect_error(
    result <- saeczi(samp,
                     pop, 
                     lin_formula = DRYBIO_AG_TPA_live_ADJ ~ tcc16,
                     domain_level = "COUNTYFIPS",
                     mse_est = TRUE,
                     B = 10L,
                     parallel = FALSE,
                     transform_fun = sqrt)
  )
  
  # incorrect inv function
  expect_warning(
    result <- saeczi(samp,
                     pop, 
                     lin_formula = DRYBIO_AG_TPA_live_ADJ ~ tcc16,
                     domain_level = "COUNTYFIPS",
                     mse_est = TRUE,
                     B = 10L,
                     parallel = FALSE,
                     transform_fun = sqrt,
                     inv_transform_fun = \(x) x^3)
  )
  
  # both properly supplied
  expect_no_error(
    result <- saeczi(samp,
                     pop, 
                     lin_formula = DRYBIO_AG_TPA_live_ADJ ~ tcc16,
                     domain_level = "COUNTYFIPS",
                     mse_est = TRUE,
                     B = 10L,
                     parallel = FALSE,
                     transform_fun = sqrt,
                     inv_transform_fun = \(x) x^2)
  )
  
})

