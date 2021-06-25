context("compute P(Doze|Wake) and P(Wake|Doze)")

asleep_sequence1 <- c(0, 1, 1, 1, 0, 1, 1, 1, 0, 0)
asleep_sequence2 <- c(0, 1, 1, 0, 1, 0, 1, 0, 1, 1)

test_that("P is computed properly", {
  expect_equal(p_doze(asleep_sequence1), 2/3)
  expect_equal(p_doze(asleep_sequence2), 4/4)
  expect_equal(p_wake(asleep_sequence1), 2/6)
  expect_equal(p_wake(asleep_sequence2), 3/5)
})
