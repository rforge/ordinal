## Make a function:
f <- function(x) sum(dnorm(x))

## compute the gradient:
gradient(fun=f, x=1:4)

## Compare with analytical gradient:
df <- deriv(~ dnorm(x1) + dnorm(x2), c("x1", "x2"), func=TRUE)
dg <- function(x1, x2) as.vector(attributes(df(x1, x2))[[1]])

gradient(f, 1:2) - dg(1, 2)

stopifnot(
    isTRUE(all.equal(dg(1, 2), gradient(f, 1:2)))
    )
