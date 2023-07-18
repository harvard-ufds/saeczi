# let's test on NV and make sure we get the same results as the sim.

# run sim on MSU server, July 8, 2023

## variable selection with help from glmnet,
## setting B = 360 and reps = 500

future::plan("multisession", workers = 8)

source("../zero-inflation-mse/R/compareEstimators.R")

wa_dat <- readRDS("../zero-inflation-mse/data/new_states/WA.rds")
nv_dat <- readRDS("../zero-inflation-mse/data/new_states/NV.rds")

store <- list()
time_l <- list()

B = 360
num_samps = 500 

states <- c("NV", "WA")

t1 <- Sys.time()
store[[1]] <- compareEstimators(
  data_list = nv_dat,
  lin_formula = DRYBIO_AG_TPA_live_ADJ_TONS ~ tcc16 + wc3cl2 + wc3cl3 + elev,
  log_formula = DRYBIO_AG_TPA_live_ADJ_TONS ~ wc3cl2 + wc3cl3 + tmean + elev,
  sqrt_response = TRUE,
  B = B,
  num_samps = num_samps
) 
t2 <- Sys.time()
time_l[[1]] <- t2 - t1
saveRDS(store[[1]], "sim_results_nv.rds")

t1 <- Sys.time()
store[[2]] <-
  compareEstimators(
    data_list = wa_dat,
    lin_formula = DRYBIO_AG_TPA_live_ADJ_TONS ~ tcc16 + wc3cl2 + wc3cl3 + def,
    log_formula = DRYBIO_AG_TPA_live_ADJ_TONS ~ wc3cl2 + wc3cl3 + tcc16 + tmin01,
    sqrt_response = TRUE,
    B = B,
    num_samps = num_samps
  )
t2 <- Sys.time()
time_l[[3]] <- t2 - t1
saveRDS(store[[3]], "sim_results_wa.rds")


names(store) <- states

saveRDS(store, "sim_results.rds")