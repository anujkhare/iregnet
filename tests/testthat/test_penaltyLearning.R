context("\npenalty learning data set")

data(penalty.learning)

fit <- with(penalty.learning, iregnet(X.mat, y.mat))
l1.norm.vec <- colSums(abs(fit$beta[-1,]))

test_that("only one model should have L1 arclength 0", {
  n.zero <- sum(l1.norm.vec == 0)
  expect_equal(n.zero, 1)
})
