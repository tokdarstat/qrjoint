\name{coef.qrjoint}
\Rdversion{2.0}
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
Extracts posterior draws of intercept and slope functions from saved draws of model parameters. A plot may be obtained if \code{plot = TRUE}. A list is returned invisibly with three fields.

\item{beta.samp}{a three-dimensional array of posterior samples of \eqn{\beta_j}{{beta}[j]}. Its dimensions are \code{L} (length of tau.grid)  x \code{p+1} (intercept + slopes) x \code{nmc} (retained posterior draws).}
\item{beta.est}{a three-dimensional array containing posterior summaries (2.5th,  50th, and 97.5th percentiles) of \eqn{\beta_j}{{beta}[j]}. Its dimensions are \code{L} (length of tau.grid)  x \code{p+1} (intercept + slopes) x 3 (posterior summaries).} 
\item{parametric}{a matrix with 3 columns giving the estimate (posterior median) and the lower and upper end points of the 95\% posterior credible interval for \eqn{\gamma_0}{gamma0}, \eqn{\gamma}{gamma}, \eqn{\sigma}{sigma}, and, \eqn{\nu}{nu}.}
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

# model fitting with 50 posterior samples from 100 iterations (thin = 2)
fit.qrj <- qrjoint(BetaPlasma ~ Age + Sex + SmokStat + Quetelet + VitUse + Calories + 
        Fat + Fiber + Alcohol + Cholesterol + BetaDiet, plasma, nsamp = 40, thin = 2)

\dontrun{
betas <- coef(fit.qrj) ## no plots
betas <- coef(fit.qrj, plot = TRUE) ## estimates are plotted
}
}
\keyword{programming}
