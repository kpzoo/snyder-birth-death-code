# Examine the Agamid dataset from laser

# Clean the workspace and console
closeAllConnections()
rm(list=ls())
cat("\014")  
graphics.off()

# Load package with agamid data
library(laser)
data("agamids")

# Set working directory
pathf <- paste(c('/Users/kris/Desktop/Zoology/Zoology2017/Lent term/Code/agamid code/'), collapse = '')
setwd(pathf)

# Function to write simple csv files to correct path
tableWrite <- function(val, name, pathname) {
  # Add path to name
  str0 <- paste(c(pathname, name), collapse = "")
  # Write table
  write.table(val, str0, row.names=FALSE, col.names=FALSE, sep=",")
}

# Write the tree and get branching times (distance of each node from the tips)
tableWrite(agamids, 'agamidsTree.txt', pathf)
agbtimes <- getBtimes(string = agamids)
tableWrite(agbtimes, 'branch.csv', pathf)

# Single run to be checked later for initial condition dependence
mod1 <- bd(agbtimes, ai=runif(2, min=0, max=1))
mod2 <- fitSPVAR(agbtimes, init=runif(3, min=0, max=100))
mod3 <- fitEXVAR(agbtimes, init=runif(3, min=0, max=100))
mod4 <- fitBOTHVAR(agbtimes, init=runif(4, min=0, max=100))


# Initialise matrix of 100 rows (for each initial condition) for each model
n = 100
const = matrix(nrow = n, ncol = 2)
spvar = matrix(nrow = n, ncol = 3)
exvar = matrix(nrow = n, ncol = 3)
bothvar = matrix(nrow = n, ncol = 4)

# Store AIC values and log max likelihoods
aic = matrix(nrow = n, ncol = 4)
lik = matrix(nrow = n, ncol = 4)

# Loop and perform optimisations at each initial condition
for(i in 1:n){
  # Constant rate BDP
  # (runif produces random numbers between limits in a list)
  mod <- bd(agbtimes, ai=runif(3, min=0, max=1))
  const[i, 1:2] = c(mod$r, mod$a)
  lik[i, 1] = mod$LH
  aic[i, 1] = mod$aic
  # Time varying birth and constant death
  mod <- fitSPVAR(agbtimes, init=runif(3, min=0, max=100))
  spvar[i, 1:3] = c(mod$lam0, mod$k, mod$mu0)
  lik[i, 2] = mod$LH
  aic[i, 2] = mod$aic
  # Constant birth, time varying death
  mod <- fitEXVAR(agbtimes, init=runif(3, min=0, max=100))
  exvar[i, 1:3] = c(mod$lam0, mod$z, mod$mu0)
  lik[i, 3] = mod$LH
  aic[i, 3] = mod$aic
  # Time varying birth and death rate
  mod <- fitBOTHVAR(agbtimes, init=runif(4, min=0, max=100))
  bothvar[i, 1:4] = c(mod$lam0, mod$k, mod$z, mod$mu0)
  lik[i, 4] = mod$LH
  aic[i, 4] = mod$aic
}

# Write model selection data
tableWrite(lik, 'lik.csv', pathf)
tableWrite(aic, 'aic.csv', pathf)

# Write complete data from repeat runs
tableWrite(const, 'const.csv', pathf)
tableWrite(spvar, 'spvar.csv', pathf)
tableWrite(exvar, 'exvar.csv', pathf)
tableWrite(bothvar, 'bothvar.csv', pathf)

# Get means of parameters for each model
constMean <- colMeans(const)
spvarMean <- colMeans(spvar)
exvarMean <- colMeans(exvar)
bothvarMean <- colMeans(bothvar)