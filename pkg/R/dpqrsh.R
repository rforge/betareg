## ----------------------------------------------------------------------------------
## beta distribution: regression specification
## (mean = mu, precision = phi)
## ----------------------------------------------------------------------------------

dbetar <- function(x, mu, phi, log = FALSE) {
  dbeta(x, shape1 = mu * phi, shape2 = (1 - mu) * phi, log = log)
}

pbetar <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  pbeta(q, shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lower.tail = lower.tail, log.p = log.p)
}

qbetar <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  qbeta(p, shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lower.tail = lower.tail, log.p = log.p)
}

rbetar <- function(n, mu, phi) {
  rbeta(n, shape1 = mu * phi, shape2 = (1 - mu) * phi)
}

sbetar <- function(x, mu, phi, parameter = c("mu", "phi"), drop = TRUE) {
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "phi")))
  xstar <- qlogis(x)
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
  s <- cbind(
    if("mu" %in% parameter) phi * (xstar - mustar),
    if("phi" %in% parameter) (mu * (xstar - mustar) + log(1 - x) - digamma((1 - mu) * phi) + digamma(phi))
  )
  colnames(s) <- c("mu", "phi")[c("mu", "phi") %in% parameter]
  if(drop) drop(s) else s
}

hbetar <- function(x, mu, phi, parameter = c("mu", "phi"), drop = TRUE) {
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "phi")))
  if(all(c("mu", "phi") %in% parameter)) parameter <- c(parameter, "mu:phi")
  n <- max(length(x), length(mu), length(phi))
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  psi1 <- trigamma(mu * phi)
  psi2 <- trigamma((1 - mu) * phi)
  a <- psi1 + psi2
  b <- psi1 * mu^2 + psi2 * (1 - mu)^2 - trigamma(phi)
  h <- cbind(
    if("mu" %in% parameter) phi^2 * (psi1 + psi2),
    if("phi" %in% parameter) psi1 * mu^2 + psi2 * (1 - mu)^2 - trigamma(phi),
    if("mu:phi" %in% parameter) phi * (mu * psi1 - (1 - mu) * psi2)
  )
  colnames(h) <- c("mu", "phi", "mu:phi")[c("mu", "phi", "mu:phi") %in% parameter]
  if(drop) drop(h) else h
}


## ----------------------------------------------------------------------------------
## (censored) four-parameter beta distribution in regression parametrization
## (mean = mu, precision = phi, support = (theta1, theta2))
## ----------------------------------------------------------------------------------

dbeta4 <- function(x, mu, phi, theta1, theta2 = 1 - theta1, log = FALSE, censored = FALSE) {
  out <- dbeta((x - theta1) / (theta2 - theta1),
    shape1 = mu * phi, shape2 = (1 - mu) * phi, log = log
  )
  out <- if(log) out - log(theta2 - theta1) else out/(theta2 - theta1)
  if(censored) {
    ## unify lengths of all variables
    n <- length(out)
    x <- rep_len(x, n)
    mu <- rep_len(mu, n)
    phi <- rep_len(phi, n)
    theta1 <- rep_len(theta1, n)
    theta2 <- rep_len(theta2, n)

    out[x <= 0] <- pbeta4(0, mu = mu[x <= 0], phi = phi[x <= 0], theta1 = theta1[x <= 0], theta2 = theta2[x <= 0], censored = FALSE, log.p = log)
    out[x >= 1] <- pbeta4(1, mu = mu[x >= 1], phi = phi[x >= 1], theta1 = theta1[x >= 1], theta2 = theta2[x >= 1], censored = FALSE, log.p = log, lower.tail = FALSE)
    out[x < 0 | x > 1] <- if(log) -Inf else 0
  }
  return(out)
}

pbeta4 <- function(q, mu, phi, theta1, theta2 = 1 - theta1, lower.tail = TRUE, log.p = FALSE, censored = FALSE) {
  out <- pbeta((q - theta1) / (theta2 - theta1),
    shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lower.tail = lower.tail, log.p = log.p
  )
  if(censored) {
    if(lower.tail) {
      out[q < 0] <- if(log.p) -Inf else 0
      out[q >= 1] <- if(log.p) 0 else 1
    } else {
      out[q < 0] <- if(log.p) 0 else 1
      out[q >= 1] <- if(log.p) -Inf else 0
    }
  }
  return(out)
}

qbeta4 <- function(p, mu, phi, theta1, theta2 = 1 - theta1, lower.tail = TRUE, log.p = FALSE, censored = FALSE) {
  q <- qbeta(p, shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lower.tail = lower.tail, log.p = log.p)
  q <- q * (theta2 - theta1) + theta1
  if(censored) {
    q[q < 0] <- 0
    q[q > 1] <- 1
  }
  return(q)
}

rbeta4 <- function(n, mu, phi, theta1, theta2 = 1 - theta1, censored = FALSE) {
  r <- theta1 + (theta2 - theta1) * rbeta(n, shape1 = mu * phi, shape2 = (1 - mu) * phi)
  if(censored) {
    r[r < 0] <- 0
    r[r > 1] <- 1
  }
  return(r)
}


## ----------------------------------------------------------------------------------
## (censored) exponential mixture beta distribution (see IWSM 2015 Proceedings)
## (mean = mu, precision = phi, contamination = nu)
## ----------------------------------------------------------------------------------

## auxiliary quadrature function
quadtable <- function(nquad = 20) {
  matrix(unlist(
    statmod::gauss.quad(n = nquad, kind = "laguerre", alpha = 0)
  ), nrow = nquad)
}

dbetax <- function(x, mu, phi, nu, log = FALSE, censored = FALSE, quad = 20) {
  ## standard beta distribution
  if(length(quad) == 1L) quad <- quadtable(quad)

  if(all(nu == 0)) return(dbeta(x, shape1 = mu * phi, shape2 = (1 - mu) * phi, log = log))

  if(!censored && any(x < 0 | x > 1)) {
    warning("computation of uncensored density outside the unit interval may be numerically unreliable")
  }

  ## unify lengths of all variables
  n <- max(length(x), length(mu), length(phi), length(nu))
  x <- rep_len(x, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  nu <- rep_len(nu, n)

  ## quadrature
  out <- apply(quad, 1, function(rule) {
    e <- rule[1] * nu
    rule[2] * dbeta((x + e)/(1 + 2 * e), shape1 = mu * phi, shape2 = (1 - mu) * phi)/(1 + 2 * e)
  })
  out <- if (is.null(dim(out))) sum(out) else rowSums(out)

  ## censoring and log transformation
  if(censored) {
    out[x <= 0] <- pbetax(0, mu = mu[x <= 0], phi = phi[x <= 0], nu = nu[x <= 0], censored = FALSE)
    out[x >= 1] <- pbetax(1, mu = mu[x >= 1], phi = phi[x >= 1], nu = nu[x >= 1], censored = FALSE, lower.tail = FALSE)
    out[x < 0 | x > 1] <- 0
  }
  if(log) out <- log(out)

  return(out)
}

pbetax <- function(q, mu, phi, nu, lower.tail = TRUE, log.p = FALSE, censored = FALSE, quad = 20) {
  if(length(quad) == 1L) quad <- quadtable(quad)
  ## standard beta distribution
  if(all(nu == 0)) return(pbeta(q, shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lower.tail = lower.tail, log.p = log.p))

  ## unify lengths of all variables
  n <- max(length(q), length(mu), length(phi), length(nu))
  q <- rep_len(q, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  nu <- rep_len(nu, n)

  ## quadrature
  out <- apply(quadrule, 1, function(rule) {
    e <- rule[1] * nu
    rule[2] * pbeta((q + e)/(1 + 2 * e), shape1 = mu * phi, shape2 = (1 - mu) * phi)
  })
  out <- if (is.null(dim(out))) sum(out) else rowSums(out)

  if(censored) {
    out[q < 0] <- 0
    out[q >= 1] <- 1
  }

  ## additional arguments
  if(!lower.tail) out <- 1 - out
  if(log.p) out <- log(out)
  return(out)
}

rbetax <- function(n, mu, phi, nu, censored = FALSE) {
  nu <- rexp(n, 1/nu)
  rbeta4(n, mu = mu, phi = phi, theta1 = -nu, theta2 = 1 + nu, censored = censored)
}

## ----------------------------------------------------------------------------------
