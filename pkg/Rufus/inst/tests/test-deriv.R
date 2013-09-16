context("Numerical derivatives")

test_that("gradient returns the right thing", {
### FIXME: What about names attributes?

    f <- function(x) sum(dnorm(x))
    df <- deriv(~ dnorm(x1) + dnorm(x2), c("x1", "x2"), func=TRUE)
    dg <- function(x1, x2) as.vector(attributes(df(x1, x2))[[1]])
    expect_equal(dg(1, 2), gradient(f, 1:2))
})

test_that("hessian returns the right thing", {

    f <- function(x) sum(dnorm(x))
    df <- deriv3(~ dnorm(x1) + dnorm(x2), c("x1", "x2"), func=TRUE)
    dg <- function(x1, x2) as.vector(attributes(df(x1, x2))[[1]])
    dH <- function(x1, x2)
        matrix(unlist(attributes(df(x1, x2))[[2]]), ncol=2)

    expect_equal(dH(1, 2), hessian(f, 1:2), tol=1e-6)
### FIXME: What about names attributes?
})

test_that("deriv12 returns the right thing", {

### FIXME: What about names attributes?

    f <- function(x) sum(dnorm(x))
    df <- deriv3(~ dnorm(x1) + dnorm(x2), c("x1", "x2"), func=TRUE)
    dg <- function(x1, x2) as.vector(attributes(df(x1, x2))[[1]])
    dH <- function(x1, x2)
        matrix(unlist(attributes(df(x1, x2))[[2]]), ncol=2)

    D12 <- deriv12(f, 1:2)
    expect_equal(dg(1, 2), D12$gradient, tol=1e-6)
    expect_equal(dH(1, 2), D12$Hessian, tol=1e-6)
})
