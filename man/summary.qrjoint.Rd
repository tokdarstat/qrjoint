\name{summary.qrjoint}
\Rdversion{0.1-1}
\alias{summary.qrjoint}
\title{Summary Method for qrjoint Model Fit}
\description{Summarize model fit, including MCMC details, for \code{\link{qrjoint}}} 
\usage{ 
\method{summary}{qrjoint}(object, ntrace = 1000, plot.dev = TRUE, more.details = FALSE, ...)
} 
\arguments{
\item{object}{a fitted model of the class 'qrjoint'.}
\item{ntrace}{number of draws to be included in trace plots}
\item{plot.dev}{logical indicator of whether to show trace plot of deviance}
\item{more.details}{logical indicating whether other details from MCMC are to be plotted}
\item{...}{a limited number of plotting controls that are passed onto the deviance plot}
}
\value{
Displays the trace of the deviance statistic. More details include trace plots of of the proximity parameter of each GP, a plot of Geweke p-values for (from \code{\link{geweke.diag}}) convergence of each model parameter and an image plot of parameter correlation. Also prints two versions of Watanabe AIC.

The following quantities are returned invisibly.
\item{deviance}{vector deviance statistic of the samples parameter draws} 
\item{pg}{a matrix with \code{nsamp} number of columns, each columns could be coerced into a matrix of dimension \code{ngrid * (p+1)}, where the columns gives the conditional posterior weights on the lambda grid values for the corresponding GP function.}       
\item{prox}{posterior draws of proximity in the form of a \code{(p+1)*nsamp} matrix.}
\item{ll}{a matrix of \code{n*nsamp} containing observation level log-likelihood contributions. Used to calculate \var{waic}, and could be used for other AIC calculations.}       
\item{waic}{Two versions of Watanabe AIC from Gelman, Hwang and Vehtari (2014).}
}

\references{Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive information criterion for Bayesian models. \emph{Stat Comput}, 24, 997-1016.}
\seealso{\code{\link{qrjoint}} and \code{\link{coef.qrjoint}}.}
\examples{ 
# Plasma data analysis

# recoding variables
data(plasma)
plasma$Sex <- as.factor(plasma$Sex)
plasma$SmokStat <- as.factor(plasma$SmokStat)
plasma$VitUse <- 3 - plasma$VitUse
plasma$VitUse <- as.factor(plasma$VitUse)

# creating predictors and response (beta carotene concentration in the plasma)
X <- model.matrix(BetaPlasma ~ Age + Sex + SmokStat + Quetelet + VitUse + Calories + 
        Fat + Fiber + Alcohol + Cholesterol + BetaDiet, data = plasma)[,-1]
Y <- plasma$BetaPlasma

# model fitting with 50 posterior samples from 100 iterations (thin = 2)
fit.qrj <- qrjoint(X, Y, 50, 2)
summary(fit.qrj, more = TRUE)
}
\keyword{programming}