\name{residuals.betareg}
\alias{residuals.betareg}

\title{Residuals Method for betareg Objects}

\description{
  Extract various types of residuals from beta regression models:
  raw response residuals (observed - fitted), Pearson residuals (raw residuals scaled by
  square root of variance function), deviance residuals (scaled log-likelihood contributions).
}

\usage{
\method{residuals}{betareg}(object, type = c("pearson", "deviance", "response"), \dots)
}

\arguments{
  \item{object}{fitted model object of class \code{"betareg"}.}
  \item{type}{character indicating type of residuals.}
  \item{\dots}{currently not used.}
}

\references{
Ferrari, S.L.P., and Cribari-Neto, F. (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
}

\seealso{\code{\link{betareg}}}

\examples{
data("GasolineYield", package = "betareg")

gy <- betareg(yield ~ gravity + pressure + temp10 + temp, data = GasolineYield)

cbind(
  residuals(gy, type = "pearson"),
  residuals(gy, type = "deviance"),
  residuals(gy, type = "response")
)
}

\keyword{regression}
