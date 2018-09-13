\name{waic}
\Rdversion{2.0}
\alias{waic}
\title{Watanabe Information Criterion}
\description{Calculates two versions of the Watanabe information criteria from MCMC draws.} 
\usage{ 
waic(logliks, print = TRUE)
} 
\arguments{
\item{logliks}{a matrix of observation level log-likelihood values, the columns are MCMC iterations and the rows are observations in the data}
\item{print}{logical whether to print the results}
}
\value{
Returns the two version of the WAIC
}

\seealso{\code{\link{summary.qrjoint}}}

\references{Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive information criterion for Bayesian models. \emph{Stat Comput}, 24, 997-1016.}

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
summary(fit.qrj, more = TRUE)

# the call to summary already shows the waic for the fitted model, it also returns 
# the observation level log-likelihood vales. To calculate waic from last 20 draws 
# we can use:

\dontrun{
summary(fit.qrj, more = TRUE)
ll <- sm$ll
waic(ll[,21:40])
}
}



\keyword{programming}
