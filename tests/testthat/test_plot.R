library(testthat)
library(iregnet)
context("\nplot")

data(prostate,package="ElemStatLearn")
pros <- subset(prostate,select=-train,train==TRUE)
y.mat <- cbind(pros$lpsa, pros$lpsa)
X.mat <- as.matrix(subset(pros, select=-lpsa))
fit <- iregnet(X.mat, y.mat)
gg <- plot(fit)
test_that("plot method returns ggplot", {
  expect_true(is.ggplot(gg))
})
