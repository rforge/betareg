plot.betareg <- function(x, which = 1:4,
  caption = c("Pearson residuals vs indices of obs.", "Deviance residuals vs indices of obs.",
    "Cook's distance plot", "Generalized leverage vs predicted values",
    "Half-normal plot of standardized residuals", "Half-normal plot of deviance residuals",
    "Pearson residuals vs linear predictor", "Deviance residuals vs linear predictor"),
    sub.caption = paste(deparse(x$call), collapse = "\n"), main = "", 
    ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
    ..., nsim = 100, level = 0.9)
{
  if(!is.numeric(which) || any(which < 1) || any(which > 8)) 
    stop("`which' must be in 1:8")
  rd <- residuals(x, type = "deviance")
  rp <- residuals(x, type = "pearson")
  n <- length(rd)
  k <- length(x$coefficients) - 1
  show <- rep(FALSE, 8)
  show[which] <- TRUE
  one.fig <- prod(par("mfcol")) == 1
  if(ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  if(show[1]) {
    plot(1:n, rp, xlab = "Obs. number", ylab = "Pearson residuals", main = main, ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[1], 3, 0.25)
    abline(h = 0, lty = 3, col = "gray")
  }
  if(show[2]) {
    plot(1:n, rd, xlab = "Indices of obs.", ylab = "Deviance residuals", main = main, ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[2], 3, 0.25)
    abline(h = 0, lty = 3, col = "gray")
  }
  if(show[3]) {
    plot(1:n, cooks.distance(x),
      xlab = "Obs. number", ylab = "Cook's distance", type = "h", main = main)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[3], 3, 0.25)
  }
  if(show[4]) {
    plot(fitted(x), gleverage(x),
      xlab = "Predicted values", ylab = "Generalized leverage", main = main, ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[4], 3, 0.25)
  }
  if(show[5]) {
    hn <- halfnormal.betareg(x, nsim = nsim, level = level, type = "pearson")
    plot(hn[,1], hn[,2], ylim = range(hn[,-1]), main = main,
      xlab = "Normal quantiles", ylab = "Absolute values of Pearson residuals", ...)
    lines(hn[,1], hn[,3],lty = 2)
    lines(hn[,1], hn[,4],lty = 1)
    lines(hn[,1], hn[,5],lty = 1)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[5], 3, 0.25)
  }
  if(show[6]) {
    hn <- halfnormal.betareg(x, nsim = nsim, level = level, type = "deviance")
    plot(hn[,1], hn[,2], ylim = range(hn[,-1]), main = main,
      xlab = "Normal quantiles", ylab = "Absolute values of deviance residuals", ...)
    lines(hn[,1], hn[,3],lty = 2)
    lines(hn[,1], hn[,4],lty = 1)
    lines(hn[,1], hn[,5],lty = 1)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[6], 3, 0.25)
  }
  if(show[7]) {
    plot(predict(x, type = "link"), rp,
      xlab = "Linear predictor", ylab = "Pearson residuals", main = main, ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[7], 3, 0.25)
    abline(h = 0, lty = 3, col = "gray")
  }
  if(show[8]) {
    plot(predict(x, type = "link"), rd,
      xlab = "Linear predictor", ylab = "Deviance residuals", main = main, ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[8], 3, 0.25)
    abline(h = 0, lty = 3, col = "gray")
  }

  if(!one.fig && par("oma")[3] >= 1) mtext(sub.caption, outer = TRUE, cex = 1.25)
  invisible()
}

halfnormal.betareg <- function(model, nsim = 100, level = 0.90, type = c("pearson", "deviance"))
{
  type <- match.arg(type)

  ## extract response y and regressors X
  y <- if(is.null(model$y)) model.response(model.frame(model)) else model$y
  x <- if(is.null(model$x)) model.matrix(model) else model$x
  offset <- if(is.null(model$offset)) rep(0, NROW(x)) else model$offset
  wts <- weights(model)

  n <- NROW(x)
  alpha <- (1 - level)/2
  mu <- fitted(model)
  phi <- tail(model$coefficients, 1)    
  res <- residuals(model, type = type)

  e <- matrix(0, n, nsim)
  e1 <- numeric(n)
  e2 <- numeric(n)
  
  for(i in 1:nsim) {
    ysim <- rbeta(n, mu * phi, (1 - mu) * phi)
    fit <- betareg(ysim ~ 0 + x, weights = wts, offset = offset)
    e[,i] <- sort(abs(residuals(fit, type = type)))
  }
  
  for(i in 1:n) {
    eo <- sort(e[i,])
    e1[i] <- quantile(eo, alpha)
    e2[i] <- quantile(eo, 1 - alpha)
  }
  
  e0 <- apply(e, 1, median)
  qq <- qnorm((n + 1:n + 0.5)/(2 * n + 1.125))
  
  cbind(qq, sort(abs(res)), e0, e1, e2)  
}
