library(testthat)
library(iregnet)
context("\nCheck error messages produced")

x <- matrix(rnorm(25), 5, 5)
y <- rbind(  
  c(-Inf, 2),
  c(-Inf, 3),
  c(-Inf, 4),
  c(-Inf, 5),
  c(-Inf, 8))
test_that(
expect_error(iregnet(x, y, family = "gaussian"), "Target matrix completely left censored. Try adding more data")
)
y <- rbind(
  c(3, Inf),
  c(4, Inf),
  c(6, Inf),
  c(1, Inf),
  c(2, Inf))
test_that(
expect_error(iregnet(x, y, family = "gaussian"), "Target matrix completely right censored. Try adding more data")
)
# Custom data which works fine
y <- rbind(
  c(5, 10),
  c(6, Inf),
  c(-Inf, 2),
  c(1, 2),
  c(-Inf, 3))

# Custom data that produces NAN's error
y <- rbind(
  c(5, 10),
  c(6, Inf),
  c(9, Inf),
  c(3, 4),
  c(10, Inf))

y <- rbind(
  c(1, 2),
  c(-Inf, 5),
  c(3, Inf),
  c(-Inf, 3),
  c(-Inf, 4))


