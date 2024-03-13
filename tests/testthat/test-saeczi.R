data(pop)
data(samp)

lin_formula <- sqrtbio ~ tcc16 

samp <- samp |> 
  dplyr::mutate(sqrtbio = sqrt(DRYBIO_AG_TPA_live_ADJ))
  

set.seed(5)
result <- saeczi(samp,
                 pop, 
                 lin_formula,
                 domain_level = "COUNTYFIPS",
                 mse_est = TRUE,
                 B = 10L,
                 parallel = FALSE)

test_that("result$res is a df", {
  expect_s3_class(result$res, "data.frame")
})

test_that("correct number of rows in result data.frame", {
  expect_equal(nrow(result$res), length(unique(pop$COUNTYFIPS)))
})

