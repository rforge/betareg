\name{betamix}
\alias{betamix}
\alias{fitted,FLXMRbeta-method}
\alias{fitted,betamix-method}
\alias{posterior,betamix,ANY-method}
\alias{clusters,betamix,ANY-method}
\alias{predict,FLXMRbeta-method}
\alias{predict,betamix-method}
\title{Finite Mixtures of Beta Regression for Rates and Proportions}
\description{
  Fit finite mixtures of beta regression models for rates and
  proportions via maximum likelihood with the EM algorithm using a
  parametrization with mean (depending through a link function on the
  covariates) and precision parameter (called phi).
}
\usage{
betamix(formula, data, k, subset, na.action, weights, offset,
        link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
        link.phi = "log", control = betareg.control(...), ID, nstart = 1,
        FLXcontrol = list(), cluster = NULL, which = "BIC",...)
}
\arguments{
  \item{formula}{symbolic description of the model (of type \code{y ~ x}
    or \code{y ~ x | z}; for details see \code{\link{betareg}}).}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{k}{a vector of integers indicating the number of components of
    the finite mixture; passed in turn to the \code{k} argument
    of \code{\link[flexmix]{stepFlexmix}}.}
  \item{weights}{optional numeric vector of weights.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor.}
  \item{link}{character specification of the link function in
    the mean model (mu). Currently, \code{"logit"}, \code{"probit"},
    \code{"cloglog"}, \code{"cauchit"}, \code{"log"}, \code{"loglog"} are supported.
    Alternatively, an object of class \code{"link-glm"} can be supplied.}
  \item{link.phi}{character specification of the link function in
    the precision model (phi). Currently, \code{"identity"},
    \code{"log"}, \code{"sqrt"} are supported. The default is \code{"log"}
    unless \code{formula} is of type \code{y ~ x} where the default is
    \code{"identity"} (for backward compatibility).
    Alternatively, an object of class \code{"link-glm"} can be supplied.}
  \item{control}{a list of control arguments specified via
    \code{\link{betareg.control}}.}
  \item{ID}{grouping variable indicating if observations are from the same
    individual, i.e. the component membership is restricted to be the
    same for these observations.}
  \item{nstart}{for each value of \code{k} run
    \code{\link[flexmix]{stepFlexmix}} \code{nstart} times and keep only
    the solution with maximum likelihood.}
  \item{FLXcontrol}{object of class \code{FLXcontrol} or a named list;
    controls the EM algorithm and passed in turn to the \code{control}
    argument of \code{\link[flexmix]{flexmix}}.}
  \item{cluster}{Either a matrix with \code{k} columns of initial
    cluster membership probabilities for each observation; or a factor
    or integer vector with the initial cluster
    assignments of observations at the start of the EM
    algorithm. Default is random assignment into \code{k}
    clusters.}
  \item{which}{Number of model to get if \code{k} is a vector of
    integers longer than one. If character, interpreted as
    number of components or name of an information criterion.}
  \item{\dots}{arguments passed to \code{\link{betareg.control}}.}
}
\details{
  The arguments and the model specification are similar to
  \code{\link{betareg}}. Internally \code{\link[flexmix]{stepFlexmix}}
  is called with suitable arguments to fit the finite mixture model with
  the EM algorithm. 
}
\value{
  An object of class \code{"flexmix"} containing the best model with
  respect to the log likelihood or the one selected according to
  \code{which} if \code{k} is a vector of integers longer than 1.
}
\author{
  Bettina Gruen and Achim Zeileis
}
\seealso{
  \code{\link{betareg}}, \code{\link[flexmix]{flexmix}},
  \code{\link[flexmix]{stepFlexmix}}
}
\examples{
data("ReadingSkills", package = "betareg")
mix_beta <- betamix(accuracy ~ iq, k = 2, data = ReadingSkills,
                    cluster = ReadingSkills$dyslexia,
                    FLXcontrol = list(verbose = 10))
}
\keyword{regression}
\keyword{cluster}
