betamix <- function(formula, data, k, subset, na.action, weights, offset,
                    link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                    link.phi = "log", 
 		    control = betareg.control(...),
                    ID, nrep = 1, FLXcontrol = list(), cluster = NULL, ...)
{
  ## Determine model.frame as in betareg
  ## add variable cluster in m

  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  
  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  } else {
    if(length(formula)[2] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- if (missing(ID)) formula else as.Formula(formula(formula), ID) 
  
  ## evaluate model.frame
  mf[[1]] <- as.name("get_all_vars")
  mf <- na.omit(eval(mf, parent.frame()))

  n <- nrow(mf)
  
  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  
  ## offset
  offset <- model.offset(mf)
  if(!is.null(offset)) {
    if(length(offset) == 1) offset <- rep(offset, n)
    offset <- as.vector(offset)
  }

  fullformula <- formula(formula, rhs = 1, lhs = 1)
  precisionformula <- if (length(Formula(oformula))[2] > 1) 
    formula(formula, lhs = 0, rhs = 2) else FALSE
  if (!missing(ID)) fullformula <- formula(as.Formula(formula(fullformula), ID))

  if (!all.equal(as.integer(weights), as.numeric(weights)))
    stop("only integer weights allowed")

  rval <- flexmix:::stepFlexmix(fullformula, data = mf, k = k, weights = as.integer(weights),
                                model = FLXMRbeta(precision = precisionformula,
                                  offset = offset, link = link, link.phi = link.phi, control = control),
                                control = FLXcontrol, nrep = nrep, cluster = mf[["(cluster)"]])
  structure(list(flexmix = rval, call = cl), class = "betamix")
}

setOldClass("betamix")

## hand-crafted "Next()" to bridge to
## exported S4 classes "flexmix", argh!
print.betamix <- function(x, ...) {
  show(x$flexmix, ...)
  invisible(x)
}
logLik.betamix <- function(object, ...) logLik(object$flexmix, ...)
summary.betamix <- function(object, ...) summary(refit(object$flexmix, ...), ...)
posterior.betamix <- function(object, newdata, ...) {
  if (missing(newdata)) return(posterior(object$flexmix, ...))
  else return(posterior(object$flexmix, newdata = newdata, ...))
}
setMethod("posterior", "betamix", posterior.betamix)
## want to have
## sctest.betamix
coef.betamix <- function(object, model = c("full", "mean", "precision"), 
                         transpose = TRUE, ...) {
  model <- match.arg(model)
  if (model == "full") {
    COEFS <- parameters(object$flexmix, ...)
  }
  else {
    COEFS <- lapply(parameters(object$flexmix, simplify = FALSE, drop = FALSE, ...),
                    function(x) {
                      z <- sapply(x, "[[", model)
                      if (!is(z, "matrix")) z <- matrix(z, ncol = object$flexmix@k,
                                                        dimnames = list(model, names(x)))
                      z})[[1]]
  }
  t(COEFS)
}

setClass("FLXMRbeta",
         representation(precision="ANY",
                        link="character",
                        link.phi="character",
                        z="matrix",
                        precision.terms="ANY",
                        precision.xlevels="ANY",
                        precision.contrasts="ANY",
                        control="ANY"),
         contains = "FLXMR")


FLXMRbeta <- function(formula = .~., precision = FALSE, offset = NULL,
                      link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                      link.phi = "log", control = betareg.control())
{
  link <- match.arg(link)
  
  if (!(is(precision, "logical") | is(precision, "formula")))
    stop("precision needs to be either a logical or a formula")

  object <- new("FLXMRbeta", weighted=TRUE, formula=formula, precision=precision,
                link = link, link.phi = link.phi, control = control,
                name=paste("FLXMRbeta(link='", link, "', link.phi='", link.phi, "')", sep = ""))

  object@defineComponent <- expression({
    predict <- function(x, z, ...) {
      dotarg = list(...)
      if("offset" %in% names(dotarg))
        offset <- dotarg$offset
      p <- x%*%coef$mean
      if (!is.null(offset)) 
        p <-  p + offset
      q <- z%*%coef$precision
      list(mean = linkobjs$mean$linkinv(drop(p)),
           precision = linkobjs$precision$linkinv(drop(q)))
    }
    
    logLik <- function(x, y, z, ...) {
      pars <- predict(x, z, ...)
      dbeta(y, shape1 = pars$mean * pars$precision, shape2 = pars$precision * (1 - pars$mean), log=TRUE)
    }
    
    new("FLXcomponent",
        parameters=coef,
        logLik=logLik, predict=predict,
        df=df)
  })
  
  object@fit <- function(x, y, z, w){
    fit <- betareg.fit(x, as.vector(y), z, weights=w, offset=offset, link = link, link.phi = link.phi,
                       control = control)
    with(list(coef = fit$coefficients,
              df = ncol(x)+ncol(z),
              linkobjs = fit$link),
         eval(object@defineComponent))
  }
  object
}

setMethod("FLXgetModelmatrix", signature(model="FLXMRbeta"),
function(model, data, formula, lhs=TRUE, ...) {
  model <- callNextMethod(model, data, formula, lhs, ...)
  if (is(model@precision, "logical")) {
    if(model@precision) {
      model@z <- model@x
      model@precision.terms <- model@terms
      model@precision.contrasts <- model@contrasts
      model@precision.xlevels <- model@xlevels
    } else {
      model@z <- matrix(1, nrow = nrow(model@x), ncol = 1)
    }
  } else {
    mt1 <- if (is.null(model@precision.terms)) terms(model@precision, data=data) else model@precision.terms
    mf <- model.frame(delete.response(mt1), data=data, na.action = NULL)
    model@precision.terms<- attr(mf, "terms")
    model@z <- model.matrix(model@precision.terms, data=mf)
    model@precision.contrasts <- attr(model@z, "contrasts")
    model@precision.xlevels <- .getXlevels(model@precision.terms, mf)
  }
  model
})

setMethod("FLXmstep", signature(model = "FLXMRbeta"), function(model, weights, ...)
{
  apply(weights, 2, function(w) model@fit(model@x, model@y, model@z, w))
})

setMethod("fitted", signature(object="FLXMRbeta"),
function(object, components, ...)
{
    z <- list()
    for(n in 1:length(components)){
      z[[n]] <- list(components[[n]]@predict(object@x, object@z))
    }
    z
})

setMethod("predict", signature(object="FLXMRbeta"), function(object, newdata, components, ...)
{
  object <- FLXgetModelmatrix(object, newdata, formula = object@fullformula, lhs = FALSE)
  lapply(components, function(comp) comp@predict(object@x, object@z, ...))
})

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRbeta"), function(model, components, ...) {
  matrix(sapply(components, function(x)
                x@logLik(model@x, model@y, model@z)), nrow = nrow(model@y))
})

setMethod("FLXreplaceParameters", signature(object="FLXMRbeta"),
function(object, components, parms) {
  Design <- flexmix:::FLXgetDesign(object, components)
  lapply(seq_along(components), function(k) {
    Parameters <- list()
    parms_k <- parms[as.logical(Design[k,])]
    for (i in seq_along(components[[k]]@parameters)) {
      Parameters[[i]] <- parms_k[seq_along(components[[k]]@parameters[[i]])]
      attributes(Parameters[[i]]) <- attributes(components[[k]]@parameters[[i]])
      parms_k <- parms_k[-seq_along(components[[k]]@parameters[[i]])]
    }
    names(Parameters) <- names(components[[k]]@parameters)
    variables <- c("x", "y", "offset", "family")
    variables <- variables[variables %in% slotNames(object)]
    for (var in variables) 
      assign(var, slot(object, var))
    with(list(coef = Parameters, df = components[[k]]@df,
              linkobjs = get("linkobjs", environment(components[[k]]@predict))),
         eval(object@defineComponent))
  })
})

setMethod("refit_optim", signature(object = "FLXMRbeta"),
function(object, components, coef, se) {
  Design <- flexmix:::FLXgetDesign(object, components)
  x <- lapply(1:nrow(Design), function(k) {
    rval <- cbind(Estimate = coef[as.logical(Design[k,])],
                  "Std. Error" = se[as.logical(Design[k,])])
    mean <- rval[seq_along(components[[k]]@parameters$mean),,drop=FALSE]
    precision <- rval[-seq_along(components[[k]]@parameters$mean),,drop=FALSE]
    rownames(mean) <- names(components[[k]]@parameters$mean)
    rownames(precision) <- names(components[[k]]@parameters$precision)
    mean_zval <- mean[,1]/mean[,2]
    precision_zval <- precision[,1]/precision[,2]
    pars <- list(mean = new("Coefmat", cbind(mean, "z value" = mean_zval,
                   "Pr(>|z|)" = 2 * pnorm(abs(mean_zval), lower.tail = FALSE))))
    pars <- c(pars, list(precision = new("Coefmat", cbind(precision, "z value" = precision_zval,
                           "Pr(>|z|)" = 2 * pnorm(abs(precision_zval), lower.tail = FALSE)))))
    pars
  })
  names(x) <- paste("Comp", seq_along(x), sep = ".")
  x
})

setMethod("existGradient", signature(object = "FLXMRbeta"),
function(object) TRUE)

setMethod("FLXgradlogLikfun", signature(object="FLXMRbeta"),
function(object, components, weights, ...) {
  lapply(seq_along(components), function(k) {
    linkobjs <- get("linkobjs", environment(components[[k]]@predict))
    pp <- components[[k]]@predict(object@x, object@z)
    ystar <- qlogis(as.vector(object@y))
    mu <- pp$mean
    phi <- pp$precision
    mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
    res <- phi * (ystar - mustar) * linkobjs$mean$mu.eta(linkobjs$mean$linkfun(mu))
    Scores_x <- weights[,k] * res * object@x
    res <- (mu * (ystar - mustar) + log(1 - as.vector(object@y)) -
            digamma((1-mu) * phi) + digamma(phi)) * linkobjs$precision$mu.eta(linkobjs$precision$linkfun(phi))
    Scores_z <- weights[,k] * res * object@z
    cbind(Scores_x, Scores_z)
  })
})                       

