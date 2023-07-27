library(tidyverse)
library(dplyr)
library(saezi)
library(testthat)
future::plan("multisession", workers = 4)

# plots <- readRDS("../zero-inflation-mse/Oregon/tsumdatp.rds")
# aux <- readRDS("../zero-inflation-mse/Oregon/pltassgn.rds")

dat <- readRDS("../zero-inflation-mse.nosync/data/ID.rds")

plots <- dat[[1]]
aux <- dat[[2]]



plots <- plots %>% 
  select("CN", "DRYBIO_AG_TPA_live_ADJ_TONS")

pop <- left_join(aux, plots, by = c("PLT_CN" ="CN"))

samp <- pop %>%
  group_by(COUNTYFIPS) %>%
  slice_sample(prop = 0.15, replace = FALSE)

lin_formula <- DRYBIO_AG_TPA_live_ADJ_TONS ~ tcc16 + tmean + tri

result <- unit_zi(samp, pop, lin_formula, domain_level = "COUNTYFIPS", B = 123, mse_est = TRUE, parallel = FALSE)

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
