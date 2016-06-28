\name{coef.qrjoint}
\Rdversion{0.1-1}
\alias{coef.qrjoint}
\title{Regression Coefficient Extraction from qrjoint Model Fit}
\description{Post process MCMC output from \code{\link{qrjoint}} to create summaries of intercept and slope function estimates} 
\usage{ 
\method{coef}{qrjoint}(object, burn.perc = 0.5, nmc = 200, plot = FALSE, show.intercept = TRUE, 
        reduce = TRUE, ...)
} 
\arguments{
\item{object}{a fitted model of the class \code{qrjoint}.}
\item{burn.perc}{a positive fraction indicating what fraction of the saved draws are to be discarded as burn-in}
\item{nmc}{integer giving the number of samples, post burn-in, to be used in Monte Carlo averaging}
\item{plot}{logical indicating if plots are to be made}
\item{show.intercept}{whether to plot the intercept curve when \code{plot = TRUE}}
\item{reduce}{logical indicating if the tail-expanded grid of tau values is to be reduced to the regular increment grid}
\item{...}{limited plotting parameters that are passed onto the call of \code{\link{getBands}}}
}
\value{
Extracts posterior draws of intercept and slope functions from saved draws of model parameters. A plot may be obtained if \code{plot = TRUE}. A list is returned invisibly with two fields.

\item{beta.samp}{a matrix with \code{nmc} many columns and \code{(p+1)*length(tau.grid)} many rows.}
\item{beta.est}{a list of length (p+1), \eqn{j}{j}-th element giving a 3-column matrix of median, 2.5th and 97.5th percentiles of the posterior distribution of \eqn{\beta_j}{{beta}[j]}}
}


\seealso{\code{\link{qrjoint}} and \code{\link{summary.qrjoint}} for model fitting under qrjoint. Also see \code{\link{getBands}} for plotting credible bands for coefficients.}


\examples{ 
## Plasma data analysis

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

\dontrun{
betas <- coef(fit.qrj) ## no plots
}
betas <- coef(fit.qrj, plot = TRUE, col = "darkgreen") ## estimates are plotted
}
\keyword{programming}