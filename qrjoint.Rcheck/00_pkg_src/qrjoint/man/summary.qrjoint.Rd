\name{summary.qrjoint}
\Rdversion{2.0}
\alias{summary.qrjoint}
\title{Summary Method for qrjoint Model Fit}
\description{Summarize model fit, including MCMC details, for \code{\link{qrjoint}}.} 
\usage{ 
\method{summary}{qrjoint}(object, ntrace = 1000, burn.perc = 0.5,
        plot.dev = TRUE, more.details = FALSE, ...)
} 
\arguments{
\item{object}{a fitted model of the class 'qrjoint'.}
\item{ntrace}{number of draws to be included in trace plots}
\item{burn.perc}{fraction of MCMC draws to be discarded as burn-in.}
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
\item{ll}{a matrix of \code{n*nsamp} containing observation-level log-likelihood contributions. Used to calculate \var{waic}, and could be used for other AIC calculations.}
\item{rp}{a matrix of \code{n*nsamp} containing observation-level estimated response proportions (i.e. t such that \eqn{y = Q(t|x)}) at each posterior draw. These may be used in lieu of residuals to assess model fit and assumption of linearity.}
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

# Model fitting with 40 posterior samples from 80 iterations (thin = 2) is for
# illustration only. For practical model fitting, increase iterations, 
# e.g. nsamp = 500, thin = 20
fit.qrj <- qrjoint(BetaPlasma ~ Age + Sex + SmokStat + Quetelet + VitUse + Calories + 
        Fat + Fiber + Alcohol + Cholesterol + BetaDiet, plasma, nsamp = 40, thin = 2)
summ <- summary(fit.qrj, more = TRUE)

\dontrun{
# Visually assess uniformity of response proportions with histogram and qqplot
# Notes: Can assess across all MCMC draws (as below) or for single iteration;
# adjustments to response proportions will be needed for censored observations
hist(summ$rp, breaks=40, freq=F)
curve(dunif(x),add=T)
qqplot(summ$rp, qunif(ppoints(length(summ$rp))),xlab="actual", ylab="theoretical")
abline(0,1)

# Visually assess linearity assumption using response proportions
# Notes: Can assess across all MCMC draws or for single iteration (as below)

# Loess gives visual of center of response proportions across covariate;
# trend line should be near 0.5
library(ggplot2)
use <- sample(1:ncol(summ$rp),1)
plasma$rpsamp <- summ$rp[,use]
ggplot(data=plasma, aes(x=Age, y=rpsamp)) + geom_point() + geom_smooth(se=F, 
method="loess")

# Violin plot allows for assessment of entire distribution across covariate;
# densities within decile bins should be blocky-uniform 
cut_dec <- function(x) factor(cut(x, quantile(x,0:10/10),inc=TRUE),labels=1:10)
ggplot(data=plasma, aes(x=cut_dec(Age), y=rpsamp)) + geom_violin() + 
xlab("Age Decile Bins")
}

}
\keyword{programming}
