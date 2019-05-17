context("compute hamming index")
set.seed(42)


asleep_sequence1 <- c(0, 1, 1, 1, 0, 1, 1, 1, 0, 0)
asleep_sequence2 <- c(0, 1, 1, 0, 1, 0, 1, 0, 1, 1)

long_asleep_sequence_file <- system.file("extdata", "long_asleep_sequence.txt", package = "flyr")

long_asleep_sequence <- read.table(file = long_asleep_sequence_file)[,1]

test_that("the hamming index is computed properly", {
  short_hi <- hamming_index(asleep_sequence1)
  long_hi <- hamming_index(long_asleep_sequence)
  expect_equal(short_hi, 0.1)
  expect_equal(long_hi, 0.24)
})

test_that("if no sleep/doze state is observed, the function returns a as.double(NA)", {
  no_hi <- hamming_index(rep(F, 10))
  expect_true(is.na(no_hi), T)
  expect_true(is.double(no_hi), T)
})
