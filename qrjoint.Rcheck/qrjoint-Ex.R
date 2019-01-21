pkgname <- "qrjoint"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "qrjoint-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('qrjoint')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("chull.center")
### * chull.center

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: chull.center
### Title: Fast Interior Point Center of Multivariate Data
### Aliases: chull.center
### Keywords: programming

### ** Examples
 
p <- 9
n <- 200
u <- runif(n)
require(splines)
x <- bs(u, df = p)
chull.center(x, plot = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("chull.center", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("coef.qde")
### * coef.qde

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: coef.qde
### Title: Coefficient Extraction from qde Model Fit
### Aliases: coef.qde
### Keywords: programming

### ** Examples

## Plasma data analysis
data(plasma)
Y <- plasma$BetaPlasma
Y <- Y + 0.1 * rnorm(length(Y)) ## remove atomicity

# model fitting with 50 posterior samples from 100 iterations (thin = 2)
fit.qde <- qde(Y, 50, 2)
betas <- coef(fit.qde)
signif(betas$parametric, 3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("coef.qde", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("coef.qrjoint")
### * coef.qrjoint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: coef.qrjoint
### Title: Regression Coefficient Extraction from qrjoint Model Fit
### Aliases: coef.qrjoint
### Keywords: programming

### ** Examples
 
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

## Not run: 
##D betas <- coef(fit.qrj) ## no plots
##D betas <- coef(fit.qrj, plot = TRUE) ## estimates are plotted
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("coef.qrjoint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getBands")
### * getBands

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getBands
### Title: Posterior Credible Bands
### Aliases: getBands
### Keywords: programming

### ** Examples
 
## toy example

x <- 2*pi*seq(0,1,.01)
fsamp <- replicate(100, rnorm(1,0,0.1) + rnorm(1,1,0.2) * cos(x))
getBands(fsamp)
getBands(fsamp, x = x, col = 3, remove.edges = FALSE, xlab = "x", ylab = "f", bty = "n")
getBands(fsamp, x = x, col = "darkgreen", remove.edges = FALSE, xlab = "x", ylab = "f")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getBands", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plasma")
### * plasma

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plasma
### Title: Plasma Concentration of Beta-Carotene and Retinol
### Aliases: plasma
### Keywords: datasets

### ** Examples

data(plasma)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plasma", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.qde")
### * predict.qde

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.qde
### Title: Posterior predictive summary for quantile-based density
###   estimation
### Aliases: predict.qde
### Keywords: programming

### ** Examples
 
# Plasma data analysis

data(plasma)
Y <- plasma$BetaPlasma
Y <- Y + 0.1 * rnorm(length(Y)) ## remove atomicity

# model fitting with 50 posterior samples from 100 iterations (thin = 2)
fit.qde <- qde(Y, 50, 2)
pred <- predict(fit.qde)
hist(Y, freq = FALSE, col = "gray", border = "white", ylim = c(0, max(pred$fest)))
matplot(pred$y, pred$fest, type="l", col=1, lty=c(2,1,2), ylab="Density", xlab="y")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.qde", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.qrjoint")
### * predict.qrjoint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.qrjoint
### Title: Posterior predictive summary for quantile estimation
### Aliases: predict.qrjoint
### Keywords: programming

### ** Examples
 
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
## Not run: 
##D fit.qrj <- qrjoint(BetaPlasma ~ Age + Sex + SmokStat + Quetelet + VitUse + Calories + 
##D         Fat + Fiber + Alcohol + Cholesterol + BetaDiet, plasma, nsamp = 40, thin = 2)
##D 
##D quants <- predict(fit.qrj)
##D matplot(fit.qrj$tau.g[fit.qrj$reg.ix], t(quants), type="l", xlab="p",
##D ylab="Quantile Function", col="lightgray", lty=1)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.qrjoint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("qde")
### * qde

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: qde
### Title: Quantiles based Density Estimation
### Aliases: qde update.qde
### Keywords: programming

### ** Examples

## Plasma data analysis

data(plasma)
Y <- plasma$BetaPlasma

# model fitting with 100 posterior samples from 200 iterations (thin = 2)
# this is of course for illustration, for practical model fitting you
# ought to try at least something like nsamp = 500, thin = 20
fit.qde <- qde(Y, nsamp = 100, thin = 2)
summary(fit.qde, more = TRUE)
pred <- predict(fit.qde)
hist(Y, freq = FALSE, col = "gray", border = "white", ylim = c(0, max(pred$fest)))
lines(pred$y, pred$fest[,2])
lines(pred$y, pred$fest[,1], lty = 2)
lines(pred$y, pred$fest[,3], lty = 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("qde", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("qrjoint")
### * qrjoint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: qrjoint
### Title: Joint Estimation of Linear Quantile Planes
### Aliases: qrjoint update.qrjoint
### Keywords: programming

### ** Examples
 
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

## Not run: 
##D # additional MCMC runs to get 10 more samples (20 additional iterations)
##D fit.qrj <- update(fit.qrj, 10)
##D summary(fit.qrj, more = TRUE)
## End(Not run)

## Not run: 
##D ### UIS data analysis (with right censoring)
##D data(uis)
##D uis.qrj <- qrjoint(Y ~ TREAT + NDT + IV3 + BECK + FRAC + 
##D                        RACE + AGE + SITE , data=uis, cens = (1 - uis$CENSOR), 
##D                      nsamp = 50, thin = 2, fix.nu = 1e5)
##D summary(uis.qrj, more = TRUE)
##D 
##D betas <- coef(uis.qrj, plot = TRUE, col = "darkgreen")
##D tau.grid <- uis.qrj$tau.g[uis.qrj$reg.ix]
##D L <- length(tau.grid)
##D beta.samp <- betas$beta.samp
##D 
##D # survival curve estimation for k randomly chosen subjects
##D n <- nrow(uis)
##D k <- 9
##D ix.sel <- sort(sample(n, k))
##D Qsel.gp <- predict(uis.qrj, newdata=uis[ix.sel,], summarize=FALSE)
##D   
##D colRGB <- col2rgb("darkgreen")/255
##D colTrans <- rgb(colRGB[1], colRGB[2], colRGB[3], alpha = 0.05)
##D par(mfrow = c(3,3), mar = c(4,3,2,1) + .1)
##D for(i in 1:k){
##D   plot(exp(apply(Qsel.gp[i,,],1,mean)), 1 - tau.grid, ty = "n", ann = FALSE, 
##D         bty = "n", xlim = exp(c(2, 8)), ylim = c(0,1), lty = 2, log = "x")
##D   for(j in 1:dim(beta.samp)[3])
##D       lines(exp(Qsel.gp[i,,j]), 1 - tau.grid, col = colTrans, lwd = 1)
##D   title(xlab = "Return time (days)", ylab = "Survival function", line = 2)
##D   title(main = bquote(Obs.Id == .(ix.sel[i])))
##D   grid()  
##D }
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("qrjoint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("redmaple")
### * redmaple

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: redmaple
### Title: Basal Areas of Red Maple Trees
### Aliases: redmaple
### Keywords: datasets

### ** Examples

data(redmaple)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("redmaple", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.qde")
### * summary.qde

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.qde
### Title: Summary Method for Quantile based Density Estimation
### Aliases: summary.qde
### Keywords: programming

### ** Examples
 
# Plasma data analysis

data(plasma)
Y <- plasma$BetaPlasma
Y <- Y + 0.1 * rnorm(length(Y)) ## remove atomicity

# model fitting with 50 posterior samples from 100 iterations (thin = 2)
fit.qde <- qde(Y, 50, 2)
summary(fit.qde, more = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.qde", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.qrjoint")
### * summary.qrjoint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.qrjoint
### Title: Summary Method for qrjoint Model Fit
### Aliases: summary.qrjoint
### Keywords: programming

### ** Examples
 
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

## Not run: 
##D # Visually assess uniformity of quantile levels with histogram and qqplot
##D # Notes: Can assess across all MCMC draws (as below) or for single iteration;
##D # adjustments to quantile levels will be needed for censored observations
##D hist(summ$ql, breaks=40, freq=F)
##D curve(dunif(x),add=T)
##D qqplot(summ$ql, qunif(ppoints(length(summ$ql))),xlab="actual", ylab="theoretical")
##D abline(0,1)
##D 
##D # Visually assess linearity assumption using quantile levels
##D # Notes: Can assess across all MCMC draws or for single iteration (as below)
##D 
##D # Loess gives visual of center of quantile levels across covariate;
##D # trend line should be near 0.5
##D library(ggplot2)
##D use <- sample(1:ncol(summ$ql),1)
##D plasma$qlsamp <- summ$ql[,use]
##D ggplot(data=plasma, aes(x=Age, y=qlsamp)) + geom_point() + geom_smooth(se=F,
##D method="loess")
##D 
##D # Violin plot allows for assessment of entire distribution across covariate;
##D # densities within decile bins should be blocky-uniform 
##D cut_dec <- function(x) factor(cut(x, quantile(x,0:10/10),inc=TRUE),labels=1:10)
##D ggplot(data=plasma, aes(x=cut_dec(Age), y=qlsamp)) + geom_violin() +
##D xlab("Age Decile Bins")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.qrjoint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("waic")
### * waic

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: waic
### Title: Watanabe Information Criterion
### Aliases: waic
### Keywords: programming

### ** Examples
 
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

## Not run: 
##D summary(fit.qrj, more = TRUE)
##D ll <- sm$ll
##D waic(ll[,21:40])
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("waic", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
