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
    \item{newdata}{an optional data frame containing variables on which to predict. If omitted, the fitted data are used.}
    \item{summarize}{a logical indicating whether the quantile functions should be summarized across posterior draws into a single estimate (TRUE) or be left as individual samples (FALSE)}
    \item{burn.perc}{a positive fraction indicating what fraction of the saved draws are to be discarded as burn-in}
    \item{nmc}{integer giving the number of samples, post burn-in, to be used in Monte Carlo averaging}
    \item{...}{currently no additional arguments are allowed}
}
\value{Either returns
\item{}{a matrix of posterior quantile-function estimates if \code{summarize=TRUE}. Dimensions are \code{n} (number of rows in predicted data) x \code{L} (length of regularized tau.grid); or}
\item{}{a three dimensional array of posterior quantile-function samples if \code{summarize=FALSE}. Dimensions are \code{n} (number of rows in predicted data) x \code{L} (length of regularized tau.grid)  x \code{nmc} (retained posterior draws).}

}

\note{
When supplying \code{newdata}, new covariate values should lie within the convex hull of the original fit's covariate space; otherwise, it is possible that extrapolated quantile functions will not obey the quantile monotonicity constraint. For information on expanding the convex hull of the original fit see Detail section of \code{\link{qrjoint}}.

Additional functionality is planned in future release to provide density function estimates given a set of input covariates.
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
\dontrun{
fit.qrj <- qrjoint(BetaPlasma ~ Age + Sex + SmokStat + Quetelet + VitUse + Calories + 
        Fat + Fiber + Alcohol + Cholesterol + BetaDiet, plasma, nsamp = 40, thin = 2)

quants <- predict(fit.qrj)
matplot(fit.qrj$tau.g[fit.qrj$reg.ix], t(quants), type="l", xlab="p",
ylab="Quantile Function", col="lightgray", lty=1)
}
}
\keyword{programming}
