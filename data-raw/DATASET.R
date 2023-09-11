## code to prepare `DATASET` dataset goes here

plots <- readRDS("./data-raw/tsumdatp.rds")
plots <- plots %>% 
  select("CN", "DRYBIO_AG_TPA_live_ADJ")

aux <- readRDS("./data-raw/pltassgn.rds")
aux <- aux %>%
  select("CN", "COUNTYFIPS",6:13)

pop <- left_join(aux, plots, by ="CN")

samp <- pop %>%
  group_by(COUNTYFIPS) %>%
  slice_sample(prop = 0.15, replace = FALSE) %>% 
  ungroup()

pop <- pop %>%
  select(-"DRYBIO_AG_TPA_live_ADJ")


usethis::use_data(pop, overwrite = TRUE)
usethis::use_data(samp, overwrite = TRUE)
