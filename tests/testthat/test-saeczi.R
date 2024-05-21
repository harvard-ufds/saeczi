data(pop)
data(samp)

lin_formula <- sqrtbio ~ tcc16 

samp <- samp |> 
  dplyr::mutate(sqrtbio = sqrt(DRYBIO_AG_TPA_live_ADJ))
  

set.seed(5)
suppressWarnings(
  result <- saeczi(samp,
                   pop, 
                   lin_formula,
                   domain_level = "COUNTYFIPS",
                   mse_est = TRUE,
                   B = 10L,
                   parallel = FALSE)
)

test_that("result$res is a df", {
  expect_s3_class(result$res, "data.frame")
})


# if (Sys.info()["sysname"] == "Darwin") {
#   test_that("result is correct", {
#     expect_snapshot(result$res$est)
#   }) 
# }

test_that("mse column is not NA", {
  expect_equal(all(!is.na(result$res$mse)), T)
})

test_that("mse column exists", {
  expect_contains(names(result$res), "mse")
})

test_that("correct number of rows in result data.frame", {
  expect_equal(nrow(result$res), length(unique(pop$COUNTYFIPS)))
})

