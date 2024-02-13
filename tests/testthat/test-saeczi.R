data(pop)
data(samp)

lin_formula <- DRYBIO_AG_TPA_live_ADJ ~ tcc16 + elev

set.seed(5)
result <- saeczi(samp,
                 pop, 
                 lin_formula,
                 domain_level = "COUNTYFIPS",
                 mse_est = TRUE,
                 B = 5,
                 parallel = FALSE)

test_that("result$res is a df", {
  expect_s3_class(result$res, "data.frame")
})



test_that("correct number of rows in result data.frame", {
  expect_equal(nrow(result$res), length(unique(pop$COUNTYFIPS)))
})

