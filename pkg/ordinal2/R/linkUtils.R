#################################
## pfun:

pgumbel <- function(q, location = 0, scale = 1, lower.tail = TRUE)
### CDF for the gumbel distribution
### Currently only unit length location and scale are supported.
    .C("pgumbel",
       q = as.double(q),
       length(q),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(lower.tail))$q

pgumbel2 <- function(q, location = 0, scale = 1, lower.tail = TRUE)
### CDF for the 'swapped' gumbel distribution
### Currently only unit length location and scale are supported.
    .C("pgumbel2",
       q = as.double(q),
       length(q),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(lower.tail))$q

pgumbelR <- function(q, location = 0, scale = 1, lower.tail = TRUE)
### R equivalent of pgumbel()
{
    q <- (q - location)/scale
    p <- exp(-exp(-q))
    if (!lower.tail) 1 - p else p
}

pgumbel2R <- function(q, location = 0, scale = 1, lower.tail = TRUE)
{
    q <- (-q - location)/scale
    p <- exp(-exp(-q))
    if (!lower.tail) p else 1 - p
}

pAOR <- function(q, lambda, lower.tail = TRUE) {
    if(lambda < 1e-6)
        stop("'lambda' has to be positive. lambda = ", lambda, " was supplied")
    p <- 1 - (lambda * exp(q) + 1)^(-1/lambda)
    if(!lower.tail) 1 - p else p
}

pAO <- function(q, lambda, lower.tail = TRUE)
    .C("pAO",
       q = as.double(q),
       length(q),
       as.double(lambda[1]),
       as.integer(lower.tail))$q


plgammaR <- function(eta, lambda, lower.tail = TRUE) {
    q <- lambda
    v <- q^(-2) * exp(q * eta)
    if(q < 0)
        p <- 1 - pgamma(v, q^(-2))
    if(q > 0)
        p <- pgamma(v, q^(-2))
    if(isTRUE(all.equal(0, q, tolerance = 1e-6)))
        p <- pnorm(eta)
    if(!lower.tail) 1 - p else p
}

plgamma <- function(eta, lambda, lower.tail = TRUE)
    .C("plgamma",
       eta = as.double(eta),
       length(eta),
       as.double(lambda[1]),
       as.integer(lower.tail[1]))$eta

#################################
## dfun:

dgumbel <- function(x, location = 0, scale = 1, log = FALSE)
### PDF for the gumbel distribution
### Currently only unit length location and scale are supported.
    .C("dgumbel",
       x = as.double(x),
       length(x),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(log))$x

dgumbel2 <- function(x, location = 0, scale = 1, log = FALSE)
### PDF for the 'swapped' gumbel distribution
### Currently only unit length location and scale are supported.
    .C("dgumbel2",
       x = as.double(x),
       length(x),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(log))$x

dgumbelR <- function(x, location = 0, scale = 1, log = FALSE)
### dgumbel in R
{
    q <- (x - location)/scale
    log.d <- -exp(-q) - q - log(scale)
    if (!log) exp(log.d) else log.d
}

dgumbel2R <- function(x, location = 0, scale = 1, log = FALSE)
{
    q <- (-x - location)/scale
    log.d <- -exp(-q) - q - log(scale)
    if (!log) exp(log.d) else log.d
}

dAOR <- function(eta, lambda, log = FALSE) {
### exp(eta) * (lambda * exp(eta) + 1)^(-1-1/lambda)
    if(lambda < 1e-6)
        stop("'lambda' has to be positive. lambda = ", lambda, " was supplied")
    log.d <- eta - (1 + 1/lambda) * log(lambda * exp(eta) + 1)
    if(!log) exp(log.d) else log.d
}

dAO <- function(eta, lambda, log = FALSE)
    .C("dAO",
       eta = as.double(eta),
       length(eta),
       as.double(lambda[1]),
       as.integer(log))$eta

dlgammaR <- function(x, lambda, log = FALSE) {
    q <- lambda
    q.2 <- q^(-2)
    qx <- q * x
    log.d <- log(abs(q)) + q.2 * log(q.2) -
        lgamma(q.2) + q.2 * (qx - exp(qx))
    if (!log) exp(log.d) else log.d
}

dlgamma <- function(x, lambda, log = FALSE)
    .C("dlgamma",
       x = as.double(x),
       length(x),
       as.double(lambda[1]),
       as.integer(log[1]))$x

#################################
## gfun:

ggumbel <- function(x)
### gradient of dgumbel(x) wrt. x
    .C("ggumbel",
       x = as.double(x),
       length(x))$x

ggumbel2 <- function(x)
### gradient of dgumbel(x) wrt. x
    .C("ggumbel2",
       x = as.double(x),
       length(x))$x

ggumbelR <- function(x){
### ggumbel in R
    q <- exp(-x)
    eq <- exp(-q)
    -eq*q + eq*q*q
}

ggumbel2R <- function(x) -ggumbelR(-x)

glogis <- function(x)
### gradient of dlogis
    .C("glogis",
       x = as.double(x),
       length(x))$x

gnorm <- function(x)
### gradient of dnorm(x) wrt. x
    .C("gnorm",
       x = as.double(x),
       length(x))$x

gcauchy <- function(x)
### gradient of dcauchy(x) wrt. x
    .C("gcauchy",
       x = as.double(x),
       length(x))$x

glogisR <- function(x) {
### glogis in R
    q <- exp(-x)
    2*q^2*(1 + q)^-3 - q*(1 + q)^-2
}

gnormR <- function(x)
### gnorm in R
    -x * dnorm(x)

gcauchyR <- function(x)
### gcauchy(x) in R
    -2*x/pi*(1+x^2)^-2

gAOR <- function(eta, lambda) {
    lex <- lambda * exp(eta)
    dAO(eta, lambda) * (1 - (1 + 1/lambda) * lex/(1 + lex))
}

gAO <- function(eta, lambda)
    .C("gAO",
       eta = as.double(eta),
       length(eta),
       as.double(lambda[1]))$eta

glgammaR <- function(x, lambda)
    (1 - exp(lambda * x))/lambda * dlgamma(x, lambda)

glgamma <- function(x, lambda)
    .C("glgamma",
       x = as.double(x),
       length(x),
       as.double(lambda[1]))$x

##################################################################
PFUN <- function(x, lambda = 1, link)
    .C("pfun",
       x = as.double(x),
       length(x),
       as.double(lambda),
       as.integer(link))$x

DFUN <- function(x, lambda = 1, link)
    .C("dfun",
       x = as.double(x),
       length(x),
       as.double(lambda),
       as.integer(link))$x

GFUN <- function(x, lambda = 1, link)
    .C("gfun",
       x = as.double(x),
       length(x),
       as.double(lambda),
       as.integer(link))$x
