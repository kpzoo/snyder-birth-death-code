# Script tests several time-varying BDP functions in non-batch form

# Assumptions and modifications
# - adds a scaled parameter version of logistic
# - includes tests for convergence and runs MCMC twice on the same tree
# - made to accommodate several functions of choice
# - times are referenced to tips but this required for MCMC
# - uses same priors as Snyder code
# - assumes 1 tree for now and writes it to csv
# - uses time varying BDP functions from Hohna document
# - true values x0 from Hohna 2014

# Clean the workspace and console
closeAllConnections()
rm(list=ls())
cat("\014")  
graphics.off()

# Packages of interest, TESS is the most important for BDP simulation
library(laser) # Rabosky code
library(ape) # Paradis etc phylo code
library(geiger) # other packages depend on this and coda
library(TreeSim) # Stadler code
library(TESS) # Hohna code
library(ggplot2) # polished plotting

# Clock the code
startTiming <- proc.time() # this starts the clock

#---------------------------Initialisation-------------------------------------------------

# Decide function to simulate
rateType = 'logistic'
print(paste(c('Simulating and estimating:', rateType)))
secRun <- TRUE # run two independent MCMCs

# Set working directory
pathf <- paste(c('/Users/kris/Desktop/tess code/', rateType, '/'), collapse = '')
setwd(pathf)

# Boolean for if start at MRCA - must be consistent, and nIters for MCMC
useMRCA <- TRUE
nIters <- 10^4
if(useMRCA){print('Using TMRCA = TRUE')}
nosave <- FALSE # control overwriting of data for testing

# Simulation parameters
nBatch <- 1
T <- 0 # pretty much a dummy variable, except for rateshift

# Function to write simple csv files to correct path
tableWrite <- function(val, name, pathname) {
  # Add path to name
  str0 <- paste(c(pathname, name), collapse = "")
  # Write table
  write.table(val, str0, row.names=FALSE, col.names=FALSE, sep=",")
}

# Control speciation and extinction functions based on rateType
switch(rateType,
       # Hohna 2014 BDP
       hohna = {
         # True value of BDP
         x0 <- c(1, 4, 1)
         # Exponential decreasing speciation
         spec <- function(t){
           y <- x0[1] + x0[2]*exp(-x0[3]*t)
           return(y)
         }
         # Constant extinction rate, +0*t keeps dimensions of t
         exti <- function(t){ 
           y <- x0[1] + 0*t
           return(y)
         }
       },
       # Logistic BDP from Paradis 2010
       logistic = {
         # True value to be estimated
         #x0 <- c(0.03, 0.2, 0.02, 0.25)
         x0 <- c(0.2, 1, 0.3, 2)
         # Logistically falling speciations
         spec <- function(t){
           y <- 1/(1 + exp(-x0[1]*t + x0[2]))
           #y <- 1/(1 + x0[2]*(x0[1]^t))
           return(y)
         }
         # Logistically falling extinctions
         exti <- function(t){ 
           y <- 1/(1 + exp(-x0[3]*t + x0[4]))
           #y <- 1/(1 + x0[4]*(x0[3]^t))
           return(y)
         }
       }
)

# No. parameters to estimate and no. taxa
numRV <- length(x0)
nTax <- 100

# Expected no. taxa at end time T
expected <- function(t){ tess.nTaxa.expected(begin = 0, t = t, end = t, lambda = spec, mu = exti,
                      MRCA = useMRCA, reconstructed = TRUE)
}

#---------------------------Simulation of tree and branching------------------------------------------------- 

# Simulate the tree of interest 
ntree <- tess.sim.taxa(n=nBatch, nTaxa=nTax, max=nTax, lambda = spec, mu = exti, MRCA=useMRCA)  
times <- branching.times(ntree[[1]])
#if(useMRCA){times <- as.numeric(sort(times))}else{times <- as.numeric(sort(c(0, times)))}

# Write branching times to csv files and tree
if(!nosave){
  tableWrite(times, 'branch.csv', pathf)
  write.nexus(ntree[[1]], file = 'tesstree.txt')
}

# Get expected no. taxa to T and use as nTax
nExp <- floor(expected(max(times)))
print(paste(c('Expected no. taxa is ', nExp), collapse='')) 


#---------------------------MCMC and likelihood-------------------------------------------------

# Bounds on parameter space to match Matlab
xmin <- 0.1*x0
#xmax <- xmin + 2*x0
xmax <- 5*x0

# Write initial conditions
xinit <- matrix(0, 5, numRV)
xinit[1, 1:numRV] <- x0
xinit[2, 1:numRV] <- xmin
xinit[3, 1:numRV] <- xmax
xinit[4, 1:3] <- c(T, nBatch, nTax) 
if(useMRCA){useMRCAbool <- 1} else {useMRCAbool <- 0}
xinit[5, 1:2] <- c(useMRCAbool, nIters) # assign certain entries

# Write branching times to csv files and tree
if(!nosave){tableWrite(xinit, 'init.csv', pathf)}

# Set uniform priors and entries into initialisation data structure
switch(rateType, 
       hohna = {
         prior_1 <- function(x) { dunif(x, min=xmin[1], max=xmax[1], log=TRUE) }
         prior_2 <- function(x) { dunif(x, min=xmin[2], max=xmax[2], log=TRUE) }
         prior_3 <- function(x) { dunif(x, min=xmin[3], max=xmax[3], log=TRUE) }
         priorsBD <- c("x1"=prior_1, "x2"=prior_2, "x3"=prior_3)
       },
       logistic = {
         prior_1 <- function(x) { dunif(x, min=xmin[1], max=xmax[1], log=TRUE) }
         prior_2 <- function(x) { dunif(x, min=xmin[2], max=xmax[2], log=TRUE) }
         prior_3 <- function(x) { dunif(x, min=xmin[3], max=xmax[3], log=TRUE) }
         prior_4 <- function(x) { dunif(x, min=xmin[4], max=xmax[4], log=TRUE) }
         priorsBD <- c("x1"=prior_1, "x2"=prior_2, "x3"=prior_3, "x4"=prior_4)
       }
)

# We now specify the speciation and extinction rates as functions and pass them
# into the likelihood, which again must be provided as a function.
switch(rateType,
       # Hohna 2014 BDP
       hohna = {
         likelihoodBD <- function(params) {
           speciation <- function(t) params[1] + params[2]*exp(-params[3]*t)
           extinction <- function(t) params[1]
           lnl <- tess.likelihood(times, lambda = speciation, mu = extinction, samplingProbability = 1.0, 
                                  MRCA = useMRCA, CONDITION = 'taxa', log = TRUE)
           return (lnl)
         }
         # Parameter for log transforms
         logT <- c(TRUE,TRUE,FALSE)
       },
       # Logistic BDP from Paradis 2010
       logistic = {
         likelihoodBD <- function(params) {
           speciation <- function(t) 1/(1 + exp(-params[1]*t + params[2]))
           extinction <- function(t) 1/(1 + exp(-params[3]*t + params[4]))
           #speciation <- function(t) 1/(1 + params[2]*(params[1]^t))
           #extinction <- function(t) 1/(1 + params[4]*(params[3]^t))
           lnl <- tess.likelihood(times, lambda = speciation, mu = extinction,
                                  samplingProbability = 1.0, MRCA = useMRCA,
                                  CONDITION = 'survival', log = TRUE)
           return (lnl)
         }
       # Parameter for log transforms
       logT <- c(FALSE,TRUE,FALSE,TRUE)       
       }
)


# Main MCMC code for estimation of parameters
samples <- tess.mcmc(likelihoodFunction = likelihoodBD, priors = priorsBD,
                           parameters = runif(numRV,0,1), logTransforms = logT,
                           delta = rep(1, numRV), iterations = nIters, burnin = nIters/4,
                           thinning = 10, adaptive = TRUE, verbose = TRUE)

# Second run for testing convergence
if(secRun){
  samples2 <- tess.mcmc(likelihoodFunction = likelihoodBD, priors = priorsBD,
                       parameters = runif(numRV,0,1), logTransforms = logT,
                       delta = rep(1, numRV), iterations = nIters, burnin = nIters/4,
                       thinning = 10, adaptive = TRUE, verbose = TRUE)
}


#---------------------------Results and plots-------------------------------------------------

# Sample statistics, ess > 100 for enough sample points for posterior
ess <- effectiveSize(samples)
xmcmc <- colMeans(samples)
# 2.5% and 97.5% percentiles for HPD
xlow <- quantile(samples, 0.025)
xhigh <- quantile(samples, 0.975)
if(secRun){
  ess2 <- effectiveSize(samples2)
  xmcmc2 <- colMeans(samples2)
  # 2.5% and 97.5% percentiles for HPD
  xlow2 <- quantile(samples2, 0.025)
  xhigh2 <- quantile(samples2, 0.975)
}

# Summarise results
su <- summary(samples)
postscript('mcmc.ps')
plot(samples)
dev.off()
print(paste(c('True values:', x0)))
print(paste(c('Estimates:', signif(xmcmc, 4))))
if(secRun){print(paste(c('Estimates 2:', signif(xmcmc2, 4))))}

# Write all the ess values, means and other statistics to csv files
if(!nosave){
  tableWrite(samples, 'samples.csv', pathf)
  if(secRun){tableWrite(samples2, 'samples2.csv', pathf)}
  tableWrite(ess, 'ess.csv', pathf)
  tableWrite(xmcmc, 'meanMCMC.csv', pathf)
}

# Plot tree, ltt, autocorrelation and mcmc densities in a single pane
#quartz() # opens a new window on mac os
if(!secRun){
  if(!nosave){postscript('rundata.ps')}
  par(mfrow=c(3, numRV)) # plotting environment of dim = mfrow
  # First put densities - numRV plots used
  densplot(samples)
  # My LTT and tree
  plot(ntree[[1]], use.edge.length = TRUE, show.tip.label = FALSE)
  plot(times, 2:nTax, 's', xlab = 'time from crown', ylab = 'no. taxa') 
  # Autocorrelation across lag
  lag <- 1:25
  corrVal <- autocorr.diag(samples, lag)
  for(i in 1:numRV){plot(lag, corrVal[1:nrow(corrVal), i], 'l', ylab = 'acf')}
  dev.off()
} else {
  quartz()
  par(mfrow=c(5, numRV)) # plotting environment of dim = mfrow
  # First put densities - numRV plots used
  densplot(samples)
  densplot(samples2)
  
  # Autocorrelation across lag
  lag <- 1:25
  corrVal <- autocorr.diag(samples, lag)
  corrVal2 <- autocorr.diag(samples2, lag)
  for(i in 1:numRV){plot(lag, corrVal[1:nrow(corrVal), i], 'l', ylab = 'acf')}
  for(i in 1:numRV){plot(lag, corrVal2[1:nrow(corrVal), i], 'l', ylab = 'acf')}
  
  # My LTT and tree - same for both samples
  plot(ntree[[1]], use.edge.length = TRUE, show.tip.label = FALSE)
  plot(times, 2:nTax, 's', xlab = 'time from crown', ylab = 'no. taxa') 
  if(!nosave){quartz.save('rundata2.pdf', type = 'pdf')}
}

# Tests for 5% significance of convergence (null hypothesis) with Gewecke
# Geweke tests if samples from early and later fractions of run taken from
# same distribution as a test of convergence
ub <- qnorm(1 - 0.05/2)
lb <- qnorm(0.05/2)
# Check for failure to convergence using the Geweke diagnostic
zscore <- geweke.diag(samples)
# Any zscore falling outside this range then reject convergence
if(any(zscore[[1]] < lb) | any(zscore[[1]] > ub)){
  # The [[1]] is as zscore is a list
  print('MCMC has not converged at 5% according to Gewecke')} else 
  {print('Cannot reject convergence hypothesis')}
if(secRun){
  zscore2 <- geweke.diag(samples2)
  if(any(zscore2[[1]] < lb) | any(zscore2[[1]] > ub)){
    print('MCMC 2 has not converged at 5% according to Gewecke')} else 
    {print('Cannot reject convergence hypothesis on 2nd run')}
}

# Gelman-Rubin test - compares the variance of sampled parameter values
# within each simulation to that between two simulations. It effectively
# assesses whether we can reject the null hypothesis that samples from
# the 2 independent MCMC simulations are drawn from the same distribution
if(secRun){
  gel <- gelman.diag(x=list(run1 = samples,run2 = samples2),
              confidence = 0.95, transform = FALSE,
              autoburnin = FALSE, multivariate=TRUE)$mpsrf
  # Want this value close to 1 for convergence
  print(paste(c('The Gelman value is: ', gel)))
}

# Display time, rounded to 4 sig figs
endTiming <- proc.time()
runTime <- (endTiming - startTiming)/60
print(paste(c('Code time in mins:', signif(runTime[3], 4))))
