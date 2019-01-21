\name{qrjoint}
\Rdversion{2.0}
\alias{qrjoint}
\alias{update.qrjoint}
\title{Joint Estimation of Linear Quantile Planes}
\description{ Estimate intercept and slope functions within a joint linear regression model of the quantiles, with possible right or left censoring of the response.
} 
\usage{ 
qrjoint(formula, data, nsamp = 1e3, thin = 10, cens = NULL, 
    wt = NULL, incr = 0.01, par = "prior", nknots = 6,
    hyper = list(sig = c(.1,.1), lam = c(6,4), kap = c(.1,.1,1)), 
    shrink = FALSE, prox.range = c(.2,.95), acpt.target = 0.15, 
    ref.size = 3, blocking = "std5", temp = 1, expo = 2, 
    blocks.mu, blocks.S, fix.nu=FALSE, fbase = c("t", "logistic", "unif"), verbose = TRUE)
    
\method{update}{qrjoint}(object, nadd, append = TRUE, ...)
} 
\arguments{
\item{formula}{an object of class "formula": a symbolic description of the model to be fitted. It must include at least one predictor.}
\item{data}{a data frame containing variables in the model.}
\item{nsamp}{number of posterior samples to be saved; defaults to 1000.}
\item{thin}{thinning rate for the Markov chain sampler -- one posterior sample is saved per \code{thin} iterations. Defaults to 10. The Markov chain sampler runs for a total of \code{nsamp * thin} many iterations.}
\item{cens}{censoring status of response. Must be a vector of length NROW(data), with 0 indicating no censoring, 1 indicating right censoring, and 2 indicating left censoring. If not supplied, defaults to all zeros.}
\item{wt}{weights attached to the observation units, expected to be non-negative numbers, and defaults to a vector of ones if not otherwise supplied.}
\item{incr}{tau grid increment. Defaults to 0.01.} 
\item{par}{character string indicating how the sampler is to be initialized. Three options are currently supported: "prior" to initialize at a random draw from the prior; "RQ" to initialize at a model space approximation of the estimates from \code{\link{rq}}; and, "noX" to initialize at a model with all slope functions being equal to zero, and the intercept function obtained by fitting \code{qde} to the response data alone.}
\item{nknots}{number of knots to be used for low rank approximation of the Gaussian process priors. Defaults to 6.}
\item{hyper}{hyperparameters of the prior distribution. Must be a list with some of all of the following fields: \code{sig}: a two vector giving the parameters of the inverse-gamma distribution on sigma-square that is used when \code{shrink}=TRUE, \code{lam}: a two vector giving the parameters of the beta distribution on proximity = \eqn{\exp(-0.01* \lambda^2)}{exp(-0.01 * lambda^2)}, and \code{kap}: a vector to be coerced into a \code{3 * nkap} matrix, with \code{nkap} being the number of components in the mixture of gamma prior on \code{kappa}, and each column of the matrix gives the shape, rate and mixing weight of a component.}
\item{shrink}{for applying shrinkage to gamma[0] and gamma. Defaults to FALSE, in which case a right Haar prior is used on (gamma[0], gamma, sigma2). If TRUE then a horseshoe type prior is used.}
\item{prox.range}{for specifying the range of length-scale parameter of the Gaussian process prior.}
\item{acpt.target}{target acceptance rate of the adaptive Metropolis sampler; defaults to 0.15}
\item{ref.size}{adaptation rate of the adaptive Metropolis sampler. The proposal density is updated once every \code{ref.size} iterations. Could be a single number or a vector of length same as the number of blocks.}
\item{blocking}{type of blocking to be applied. Either a character string specifying one to be chosen from the supplied menu (see Details), or a list giving user specified blocks. In the latter case, each element of the list is a logical vector of length equal to the total number of model parameters, which equals \code{(nknots+1)*(ncol(X)+1) + 2} indicating which model parameters belong to the block.}
\item{temp}{temperature of the log-likelihood function. The log-likelihood function is raised to the power of \code{temp}. Defaults to 1.}
\item{expo}{the exponent to be used in the covariance kernel of the Gaussian process priors. Defaults to 2, giving the standard squared-exponential covariance kernel.}
\item{blocks.mu}{initial block specific means in the form of  a list. If left unspecified then will be automatically generated as a list of vectors of zeros of appropriate lengths matching the corresponding block sizes.}
\item{blocks.S}{initial block specific covariance matrices in the form of a list. If left unspecified then will be automatically generated as a list of identity matrices of appropriate dimensions matching the corresponding block sizes. When \code{blocking} is chosen as one of the menu items of the form "std*", known prior covariance information and estimated variance matrices from \code{\link{rq}} are used.}
\item{fix.nu}{either the logical FALSE indicating that nu should be learned, or a positive real number giving the fixed value of nu, which is then excluded from MCMC updates}
\item{fbase}{either "t" (default), "logistic" or "unif" to indicate what base distribution is to be used.}
\item{verbose}{logical indicating whether MCMC progress should be printed, defaults to TRUE}
\item{object}{a fitted model of the class 'qrjoint'.}
\item{nadd}{number of additional MCMC samples.}
\item{append}{logical indicating whether new samples should be appended to old ones. If FALSE then old samples are discarded.}
\item{...}{no additional arguments are allowed}
}
\value{
\code{qrjoint(x, y, ...)} returns a `qrjoint' class object to be used by \code{\link{update.qrjoint}}, \code{\link{coef.qrjoint}} and \code{\link{summary.qrjoint}}. 

\code{update(object,...)} runs additional Markov chain iterations and appends thinned draws to an existing `qrjoint' object \code{object}. All relevant details are inherited from \code{object}.

Returned object is a list containing the following variables.

\item{par}{latest draw of the parameter vector}
\item{x}{scaled and centered design matrix}
\item{y}{response vector}    
\item{cens}{censoring status vector, 0=uncensored, 1=right censored, 2=left censored}
\item{wt}{vector of observation weights}
\item{shrink}{shrinkage indicator}
\item{hyper}{completed list of hyper-parameters}
\item{dim}{model dimension vector of the form c(n, p, length of tau grid, position of \eqn{\tau_0}{tau0} on the grid, nknots, length of lambda grid, nkap, total number of MCMC iterations, thin, nsamp)}
\item{gridmats}{details of covariance matrix factors etc, intended for internal use.}
\item{tau.g}{the tau grid}
\item{muV}{list of means for parameter blocks}
\item{SV}{list of covariance matrices for parameter blocks}
\item{blocks}{list of blocks}
\item{blocks.size}{vector of block lengths}
\item{dmcmcpar}{numeric vector containing details of adaptive MCMC runs, equals c(temp, decay rate of adaptation, vector of target acceptance rates for the blocks, vector of increment scales used in adaptation). Intended strictly for internal use.}
\item{imcmcpar}{numeric vector containing details of adaptive MCMC runs, equals c(number of parameter blocks, ref.size, indicator on whether details are to be printed during MCMC progress, rate of details printing, a vector of counters needed for printing). Intended strictly for internal use.}
\item{parsamp}{a long vector containing the parameter draws. Could be coerced into a matrix of dim \code{npar * nsamp}. Intended primarily for use by \code{\link{summary.qrjoint}} and \code{\link{coef.qrjoint}}.}
\item{acptsamp}{a long vector containing rates of acceptance statistics for parameter blocks. Could be coerced into a matrix of dim \code{nblocks * nsamp}. Not very informative, because thinning times and adaptation times may not be exactly synced.}
\item{lpsamp}{vector of log posterior values for the saved MCMC draws.}
\item{fbase.choice}{integer 1 for "t", 2 for "logistic" and 3 "unif" base.}
\item{prox}{vector of proximity (exp(-0.01*lambda^2)) grid values}
\item{reg.ix}{positions of the regular tau grid on the expanded tail-appended grid}
\item{runtime}{run time of the MCMC}
\item{call}{original model call}
\item{terms}{terms included in model frame}
}

\details{

A formula has an implied intercept term. This model requires that the intercept term be included; therefore, it cannot be explicitely removed via (y ~ 0 + x) or (y ~ -1 + x) constructs.

The model assumes each conditional quantile of the response is a hyper-plane. The intercept and slope functions (as functons of the quantile level) are estimated under the constraint that the resulting quantile planes are non-crossing over some closed, convex predictor domain. The domain is equated, by default, to the convex hull of the observed predictor vectors. The input argument \code{wt} provides more flexibility in the domain specification. All observation units are used in calculating the convex hull, but only those with non-zero weights are used in the likelihood evaluation. By including pseudo-points with zero weight (e.g. covariates from a test dataframe), the boundaries of the predictor domain can be expanded.

In running the MCMC, the following menu choices are available for blocking the parameter vector. Below, \code{p = ncol(X)}.

\code{"single"}: a single block containing all parameters

\code{"single2"}: one block containing all parameters and an additional block containing only (gamma[0], gamma, sigma, nu)

\code{"single3"}: like \code{"single2"}, but the second block is split into two further blocks, one with \eqn{(\gamma_0, \gamma)}{(gamma[0], gamma)}, the other with \eqn{(\sigma, \nu)}{(sigma, nu)}

\code{"std0"}: \code{p+1} blocks, \eqn{(j+1)}{(j+1)}-th contains \eqn{(W_{*j}, \gamma_j, \sigma, \nu)}{(W*[j], gamma[j], sigma, nu)}, \eqn{j= 0,\ldots ,p}{j = 0, ..., p}.

\code{"std1"}: total \code{p+2} blocks. First \code{p+1} blocks same as \code{"std0"} and one additional block of \eqn{(\gamma_0, \gamma, \sigma, \nu)}{(gamma[0], gamma, sigma, nu)}.

\code{"std2"}: total \code{p+3} blocks. First \code{p+1} blocks same as \code{"std0"} and two additional blocks of \eqn{(\gamma_0, \gamma)}{(gamma[0], gamma)} and \eqn{(\sigma, \nu)}{(sigma, nu)}

\code{"std3"}: total \code{p+3} blocks. First \code{p+1} blocks are \eqn{W_{*j}}{W*[j]}, \eqn{j = 0, \ldots, p}{j=0,...,p}, last two are \eqn{(\gamma_0, \gamma)}{(gamma[0], gamma)} and \eqn{(\sigma, \nu)}{(sigma, nu)}

\code{"std4"}: total \code{p+3} blocks. First \code{p+1} blocks are \eqn{(W_{*j}, \gamma_j)}{(W*[j], gamma[j])}, \eqn{j = 0, \ldots, p}{j=0,...,p}, last two are \eqn{(\gamma_0, \gamma)}{(gamma[0], gamma)} and \eqn{(\sigma, \nu)}{(sigma, nu)}

\code{"std5"}: total \code{p+4} blocks. First \code{p+3} are same as \code{"std4"} and one additional block containing all parameters.
}
\references{Yang, Y. and Tokdar, S.T., 2017. Joint estimation of quantile planes over arbitrary predictor spaces. Journal of the American Statistical Association, 112(519), pp.1107-1120.}
\seealso{\code{\link{summary.qrjoint}} and \code{\link{coef.qrjoint}}.}
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

\dontrun{
# additional MCMC runs to get 10 more samples (20 additional iterations)
fit.qrj <- update(fit.qrj, 10)
summary(fit.qrj, more = TRUE)
}

\dontrun{
### UIS data analysis (with right censoring)
data(uis)
uis.qrj <- qrjoint(Y ~ TREAT + NDT + IV3 + BECK + FRAC + 
                       RACE + AGE + SITE , data=uis, cens = (1 - uis$CENSOR), 
                     nsamp = 50, thin = 2, fix.nu = 1e5)
summary(uis.qrj, more = TRUE)

betas <- coef(uis.qrj, plot = TRUE, col = "darkgreen")
tau.grid <- uis.qrj$tau.g[uis.qrj$reg.ix]
L <- length(tau.grid)
beta.samp <- betas$beta.samp

# survival curve estimation for k randomly chosen subjects
n <- nrow(uis)
k <- 9
ix.sel <- sort(sample(n, k))
Qsel.gp <- predict(uis.qrj, newdata=uis[ix.sel,], summarize=FALSE)
  
colRGB <- col2rgb("darkgreen")/255
colTrans <- rgb(colRGB[1], colRGB[2], colRGB[3], alpha = 0.05)
par(mfrow = c(3,3), mar = c(4,3,2,1) + .1)
for(i in 1:k){
  plot(exp(apply(Qsel.gp[i,,],1,mean)), 1 - tau.grid, ty = "n", ann = FALSE, 
        bty = "n", xlim = exp(c(2, 8)), ylim = c(0,1), lty = 2, log = "x")
  for(j in 1:dim(beta.samp)[3])
      lines(exp(Qsel.gp[i,,j]), 1 - tau.grid, col = colTrans, lwd = 1)
  title(xlab = "Return time (days)", ylab = "Survival function", line = 2)
  title(main = bquote(Obs.Id == .(ix.sel[i])))
  grid()  
}
}
}
\keyword{programming}
