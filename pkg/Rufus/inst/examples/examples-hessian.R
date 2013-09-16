## Make a function:
f <- function(x) sum(dnorm(x))

## compute the Hessian (2nd derivative):
hessian(fun=f, x=1:2)

## Compare with analytical Hessian:
df <- deriv3(~ dnorm(x1) + dnorm(x2), c("x1", "x2"), func=TRUE)
dH <- function(x1, x2)
    matrix(unlist(attributes(df(x1, x2))[[2]]), ncol=2)

dH(1, 2) - hessian(f, 1:2)

stopifnot(
    isTRUE(all.equal(dH(1, 2), Rufus::hessian(f, 1:2), tol=1e-6))
    )

