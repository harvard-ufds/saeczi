data(samp)

res <- boot_rep(samp,
                "COUNTYFIPS",
                DRYBIO_AG_TPA_live_ADJ ~ tcc16 + (1 | COUNTYFIPS),
                DRYBIO_AG_TPA_live_ADJ != 0 ~ tcc16 + (1 | COUNTYFIPS))

names_nz <- samp[samp$DRYBIO_AG_TPA_live_ADJ > 0, ]$COUNTYFIPS |> unique()
names_all <- samp$COUNTYFIPS |> unique()


test_that("return object includes proper values", {
  expect_equal(names(res), c("beta_lm", "beta_glm", "u_lm", "u_glm"))
})


test_that("u_lm and u_glm are properly named vectors", {
  expect_named(res$u_lm, names_nz)
  expect_named(res$u_glm, names_all)
})
