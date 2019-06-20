context("compute hamming index")
set.seed(42)


asleep_sequence1 <- c(0, 1, 1, 1, 0, 1, 1, 1, 0, 0)
long_asleep_sequence <- c(0,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,0,1,1,1,0)


test_that("the hamming index is computed properly", {
  short_hi <- hamming_index(asleep_sequence1)
  long_hi <- hamming_index(long_asleep_sequence)
  expect_equal(short_hi, 0.1)
  expect_equal(long_hi, 0.1875)
})

test_that("if no sleep/doze state is observed, the function returns a as.double(NA)", {
  no_hi <- hamming_index(rep(F, 10))
  expect_true(is.na(no_hi), T)
  expect_true(is.double(no_hi), T)
})
