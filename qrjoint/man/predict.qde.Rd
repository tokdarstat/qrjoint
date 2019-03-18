\name{predict.qde}
\Rdversion{2.0}
\alias{predict.qde}
\title{Posterior predictive summary for quantile-based density estimation}
\description{Extract posterior predictive density estimate for \code{\link{qde}}}
\usage{ 
\method{predict}{qde}(object, burn.perc = 0.5, nmc = 200, yRange = range(object$y), yLength = 401, ...)
} 
\arguments{
    \item{object}{a fitted model of the class 'qde'.}
    \item{burn.perc}{a positive fraction indicating what fraction of the saved draws are to be discarded as burn-in}
    \item{nmc}{integer giving the number of samples, post burn-in, to be used in Monte Carlo averaging}
    \item{yRange}{range of values over which posterior predictive density is to be evaluated}
    \item{yLength}{number of grid points spanning yRange for posterior predictive density evaluation}
    \item{...}{currently no additional arguments are allowed}

}
\value{
Returns a list with three items:
\item{y}{vector giving the grid over which the posterior predictive density is evaluated.}
\item{fsamp}{a matrix with \code{yLength} many rows and \code{nmc} many columns. Each column corresponds to a draw of the response density from the posterior predictive.}
\item{fest}{summary of the posterior predictive density given by point-wise median, 2.5th and 97.5th percentiles.}
}

\seealso{\code{\link{qde}} and \code{\link{summary.qde}}.}
\examples{ 
# Plasma data analysis

data(plasma)
Y <- plasma$BetaPlasma
Y <- Y + 0.1 * rnorm(length(Y)) ## remove atomicity

# model fitting with 50 posterior samples from 100 iterations (thin = 2)
fit.qde <- qde(Y, 50, 2)
pred <- predict(fit.qde)
hist(Y, freq = FALSE, col = "gray", border = "white", ylim = c(0, max(pred$fest)))
matplot(pred$y, pred$fest, type="l", col=1, lty=c(2,1,2), ylab="Density", xlab="y")
}
\keyword{programming}
