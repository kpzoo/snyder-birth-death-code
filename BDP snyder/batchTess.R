# Batch version of tessSingle

# Hohna 2014 exponentially decreasing diversification rate function

# Assumptions and modifications
# - while loop ensures only runs that pass the Gewecke test are used
# - replaced condition on taxa to CONDITION = 'survival'
# - uses Gewecke statistic to test MCMC convergence
# - made to accommodate several functions of choice
# - writes multiple trees into a single csv
# - times are referenced to tips but this required for MCMC
# - uses priors that are exported to Snyder code
# - uses time varying BDP functions from Hohna document

# Clean the workspace and console
closeAllConnections()
rm(list=ls())
cat("\014")  
graphics.off()

# Packages of interest, TESS is the most important for BDP simulation
library(laser)
library(ape)
library(geiger)
library(TreeSim)
library(TESS)
library(ggplot2) # polished plotting

# Clock the code
startTiming <- proc.time() # this starts the clock

#---------------------------Initialisation-------------------------------------------------

# Decide function to simulate and likelihood condition
rateType = 'hohna'
print(paste(c('Simulating and estimating:', rateType)))
cond <- 'survival'
print(paste(c('Conditioning likelihood on:', cond)))

# Set working directory
pathf <- paste(c('/Users/kris/Desktop/tess code/', rateType, '/batch/'), collapse = '')
setwd(pathf)

# Boolean for if start at MRCA - must be consistent, and nIters for MCMC
useMRCA <- TRUE
nIters <- 10^4
if(useMRCA){print('Using TMRCA = TRUE')}

# True parameters
nBatch <- 100
T <- 5 # pretty much a dummy variable, except for rateshift

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
           return(y)
         }
         # Logistically falling extinctions
         exti <- function(t){ 
           y <- 1/(1 + exp(-x0[3]*t + x0[4]))
           return(y)
         }
       }
)

# No. parameters to estimate and no. taxa
numRV <- length(x0)
nTax <- 100

# Expected no. taxa at end time T
expected <- function(t){
  tess.nTaxa.expected(begin = 0, t = t, end = t, lambda = spec, mu = exti,
                      MRCA = useMRCA, reconstructed = TRUE)
}

#---------------------------Tree and branching times-------------------------------------------------

# Simulate the trees of interest
ntree <- tess.sim.taxa(n=nBatch, nTaxa=nTax, max=nTax, lambda = spec, mu = exti, MRCA=useMRCA)  
# Convert batch of trees into branching times, no sorting or 0 adding <--------
brn <- matrix(nrow = nBatch, ncol = nTax-1)
for (i in 1:nBatch) {
  # Sims holding n constant, add 0 as a branching time and sort
  brn[i,] <- branching.times(ntree[[i]])
}

# Write branching times to single csv file
tableWrite(brn, 'branchBatch.csv', pathf)
# Write trees to single txt file
write.nexus(ntree, file = 'batchtree.txt')


#---------------------------MCMC and likelihood-------------------------------------------------

# Uniform priors on parameters to match Matlab
xmin <- 0.1*x0
#xmax <- 3*x0
xmax <- 5*x0

# Write initial conditions
xinit <- matrix(0, 5, numRV)
xinit[1, 1:numRV] <- x0
xinit[2, 1:numRV] <- xmin
xinit[3, 1:numRV] <- xmax
xinit[4, 1:3] <- c(T, nBatch, nTax) 
if(useMRCA){useMRCAbool <- 1} else {useMRCAbool <- 0}
xinit[5, 1:2] <- c(useMRCAbool, nIters) # assign certain entries

# Write initial conditions to csv
tableWrite(xinit, 'init.csv', pathf)

# Set uniform priors and entries into initialisation data structure
switch(rateType, 
       hohna = {
         prior_1 <- function(x) { dunif(x, min=xmin[1], max=xmax[1], log=TRUE) }
         prior_2 <- function(x) { dunif(x, min=xmin[2], max=xmax[2], log=TRUE) }
         prior_3 <- function(x) { dunif(x, min=xmin[3], max=xmax[3], log=TRUE) }
         priorsBD <- c("x1"=prior_1, "x2"=prior_2, "x3"=prior_3)
       },
       rateshift = {
         prior_1 <- function(x) { dunif(x, min=xmin[1], max=xmax[1], log=TRUE) }
         prior_2 <- function(x) { dunif(x, min=xmin[2], max=xmax[2], log=TRUE) }
         prior_3 <- function(x) { dunif(x, min=xmin[3], max=xmax[3], log=TRUE) }
         prior_4 <- function(x) { dunif(x, min=xmin[4], max=xmax[4], log=TRUE) }
         priorsBD <- c("x1"=prior_1, "x2"=prior_2, "x3"=prior_3, "x4"=prior_4)
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
                                  MRCA = useMRCA, CONDITION = cond, log = TRUE)
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
           lnl <- tess.likelihood(times, lambda = speciation, mu = extinction,
                                  samplingProbability = 1.0, MRCA = useMRCA,
                                  CONDITION = cond, log = TRUE)
           return (lnl)
         }
         # Parameter for log transforms
         logT <- c(FALSE,TRUE,FALSE,TRUE)        
       }
)

# Storage variables for runs
ess <- matrix(nrow = nBatch, ncol = numRV)
xmcmc <- matrix(nrow = nBatch, ncol = numRV)
xlow <- matrix(nrow = nBatch, ncol = numRV)
xhigh <- matrix(nrow = nBatch, ncol = numRV)

# Bounds for Geewecke
ub <- qnorm(1 - 0.05/2)
lb <- qnorm(0.05/2)
iFail <- 0*(1:nBatch)

# Main MCMC code for estimation of parameters done for each branching time set
i = 1
while(i <= nBatch){
  # Current branching time set
  times <- brn[i,]
  samples <- tess.mcmc(likelihoodFunction = likelihoodBD, priors = priorsBD,
                       parameters = runif(numRV,0,1), logTransforms = logT,
                       delta = rep(1, numRV), iterations = nIters, burnin = nIters/4,
                       thinning = 10, adaptive = TRUE, verbose = TRUE)
  # Get ess 
  essTemp <- effectiveSize(samples)
  
  # Get Gewecke statistic and print if satisfied
  zscore <- geweke.diag(samples)
  # Condition in case run fails
  if(!any(essTemp == 0)){
    if(any(zscore[[1]] < lb) | any(zscore[[1]] > ub)){
      # The [[1]] is as zscore is a list
      print('MCMC has not converged at 5% according to Gewecke')} 
    else{
      print('Cannot reject convergence hypothesis')
      # Get ess for each tree and mean values
      ess[i,] <- essTemp
      xmcmc[i,] <- colMeans(samples)
      # 2.5% and 97.5% percentiles for HPD
      for(j in 1:numRV){
        xlow[i, j] <- quantile(samples[,j], 0.025)
        xhigh[i, j] <- quantile(samples[,j], 0.975)
      }  
      # Write samples and ess to csv files
      nameSam <- paste(c('samples_', i, '.csv'), collapse = '')
      tableWrite(samples, nameSam, pathf)
      # Print the run that was completed
      print(paste(c('Completed run: ', i, ' of ', nBatch), collapse = ''))
      # Update index
      i = i + 1
      iFail[i] = 0
    }
  } 
  else 
  {
    # Failed run should eventually be corrected
    iFail[i] = 1
  }
}

# Write all the ess values, means and other statistics to csv files
tableWrite(ess, 'ess.csv', pathf)
tableWrite(xmcmc, 'meanMCMC.csv', pathf)
tableWrite(xlow, 'lowMCMC.csv', pathf)
tableWrite(xhigh, 'highMCMC.csv', pathf)
