## Make a function:
f <- function(x) sum(dnorm(x))

## compute gradient and Hessian (first two derivatives):
D12 <- deriv12(fun=f, x=1:2)

## Compare with analytical Hessian:
df <- deriv3(~ dnorm(x1) + dnorm(x2), c("x1", "x2"), func=TRUE)
dg <- function(x1, x2) as.vector(attributes(df(x1, x2))[[1]])
dH <- function(x1, x2)
    matrix(unlist(attributes(df(x1, x2))[[2]]), ncol=2)

dg(1, 2) - D12$gradient
dH(1, 2) - D12$Hessian

stopifnot(
    isTRUE(all.equal(dg(1, 2), D12$gradient, tol=1e-6)),
    isTRUE(all.equal(dH(1, 2), D12$Hessian, tol=1e-6))
    )
