context("Gradients of density functions")

glogisR <- function(x) {
### glogis in R
  res <- rep(0, length(x))
  isFinite <- !is.infinite(x)

  x <- x[isFinite]
  isNegative <- x < 0
  q <- exp(-abs(x))
  q <- 2*q^2*(1 + q)^-3 - q*(1 + q)^-2
  q[isNegative] <- -q[isNegative]
  res[isFinite] <- q
  res
}

gnormR <- function(x)
### gnorm in R
    -x * dnorm(x)

gcauchyR <- function(x)
### gcauchy(x) in R
    -2*x/pi*(1+x^2)^-2

test_that("gradient functions return the same as their R-equivalents", {
    x <- -10:10
    expect_equal(glogis(x), glogisR(x))
    expect_equal(gnorm(x), gnormR(x))
    expect_equal(gcauchy(x), gcauchyR(x))
})
