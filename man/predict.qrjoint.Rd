\name{predict.qrjoint}
\Rdversion{2.0}
\alias{predict.qrjoint}
\title{Posterior predictive summary for quantile estimation}
\description{Extract quantile functions for \code{\link{qrjoint}}}
\usage{ 
\method{predict}{qrjoint}(object, newdata=NULL, summarize=TRUE, burn.perc = 0.5, nmc = 200, ...)
} 
\arguments{
    \item{object}{a fitted model of the class 'qrjoint'.}
    \item{newdata}{an optional data frame in which to look for variables with which to predict. If omitted, the fitted data are used.}
    \item{burn.perc}{a positive fraction indicating what fraction of the saved draws are to be discarded as burn-in}
    \item{nmc}{integer giving the number of samples, post burn-in, to be used in Monte Carlo averaging}
    \item{...}{currently no additional arguments are allowed}
}
\value{Either returns
\item{}{a matrix of posterior quantile-function medians if \code{summarize=TRUE}. Dimensions are \code{n} (number of rows in predicted data) x \code{L} (length of regularized tau.grid); or}
\item{}{a three dimensional array of posterior quantile-function samples if \code{summarize=FALSE}. Dimensions are \code{n} (number of rows in predicted data) x \code{L} (length of regularized tau.grid)  x \code{nmc} (retained posterior draws).}

}

\note{
When supplying \code{newdata}, new covariate values should lie within the convex hull of the original observed covariate space; otherwise, it is possible that extrapolated quantile functions will not obey the quantile monotonicity constraint.

Additional functionality is planned in future release to provide density function estimates given a set of input x.
}

\seealso{\code{\link{qrjoint}} and \code{\link{summary.qrjoint}}.}
\examples{ 
## Plasma data analysis

# recoding variables
data(plasma)
plasma$Sex <- as.factor(plasma$Sex)
plasma$SmokStat <- as.factor(plasma$SmokStat)
plasma$VitUse <- 3 - plasma$VitUse
plasma$VitUse <- as.factor(plasma$VitUse)

# Model fitting with 40 posterior samples from 80 iterations (thin = 2) is for
# illustration only. For practical model fitting, increase iterations, 
# e.g. nsamp = 500, thin = 20
fit.qrj <- qrjoint(BetaPlasma ~ Age + Sex + SmokStat + Quetelet + VitUse + Calories + 
        Fat + Fiber + Alcohol + Cholesterol + BetaDiet, plasma, nsamp = 40, thin = 2)
summary(fit.qrj, more = TRUE)

quants <- predict(fit.qrj)
matplot(fit.qrj$tau.g[fit.qrj$reg.ix], t(quants), type="l", xlab="p", ylab="Quantile Function", col="lightgray", lty=1)

}
\keyword{programming}
