data(samp)
data(pop)

res <- samp_by_grp(samp, pop, "COUNTYFIPS", 10)

test_that("plots per group are correct", {
  out <- vector(mode = "logical", length = 10)
  chk <- as.data.frame(table(samp$COUNTYFIPS))
  for (i in 1:10) {
    cmp <- as.data.frame(table(res[[i]]$COUNTYFIPS))
    mtch <- all.equal(chk$Freq, cmp$Freq)
    out[i] <- mtch
  }
  expect_true(all(out))
})

