
## Make artificial data:
df <- data.frame(g = factor(letters[rep(1:2, c(4, 2))]),
                 f = factor(LETTERS[c(1,2,1,2,1,1)]))


## Make design matrix:
X <- model.matrix(~ g * f, data=df)
## Compare ncol before and after dropping redundant columns:
X
newX <- drop.coef(X)
newX
## Optionally suppres message (e.g. inside functions):
newX <- drop.coef(X, silent=TRUE)

## Essentially this is being computed:
qr.X <- qr(X, tol = 1e-7, LAPACK = FALSE)
newX <- X[, qr.X$pivot[1:qr.X$rank], drop = FALSE]
## is newX of full column rank?
ncol(newX) == qr(newX)$rank
## the number of columns being dropped:
ncol(X) - ncol(newX)

## A 'real' example:
## if(require(ordinal)) {
##     X <- model.matrix( ~ PRODID * DAY, data = soup)
##     ncol(X)
##     newX <- drop.coef(X)
##     ncol(newX)
## }
