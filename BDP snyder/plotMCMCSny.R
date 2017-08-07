# Script to plot Snyder and MCMC samples comparatively

# Assumptions and modifications
# - assumes that datapath has 2 csv files in it

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
library(gridExtra) # grid options on plots
library(cowplot) # add on to ggplot for publication ready plots

# Set working directory
pathf <- '/Users/kris/Desktop/tess code/rateshift/'
setwd(pathf)
print(paste(c('Working directory:', pathf), collapse = ''))

# Ensure both files available then extract data and plot
if(file.exists('sampSny.csv') & file.exists('samples.csv')){
  # Snyder and MCMC samples extracted
  sampSny <- read.csv('sampSny.csv', sep = ',')
  samples <- read.csv('samples.csv', sep = ',')
  # Check data dimensions the same
  nc1 <- ncol(samples)
  nr1 <- nrow(samples)
  nc2 <- ncol(sampSny)
  nr2 <- nrow(sampSny)
  if(nc1 != nc2 | nr1 != nr2){
    print('Error in dimensions')}
  else{
    # No. parameters in the BDP and no. samples
    numRV <- nc1
    nSamp <- nr1
    
    # Function to set data frame for ggplot on a given variable
    setdataframe <- function(mcmcData, snyData, nS, id){
      # Input data must be a single stream for a variable
      if(length(mcmcData) != length(snyData) & length(snyData) != nS){
        print('Dimensions of inputs not sensible')
        plt <- 0
      } else{
        # Create a labelled matrix by columns
        m <- matrix(nrow = nS, ncol = 2)
        m[,1] = mcmcData
        m[,2] = snyData # assign by columns
        colnames(m) = c('mcmc', 'snyder') # column names
        # Convert to data frame
        df <- as.data.frame(m)
        # Want data in 1 column, name in 2nd so stack
        dfs <- stack(df)
        # Values labelled values and named by ind in dfs from stack
        # Now group by ind and draw smoothed densities with geom_density
        plt <- ggplot(dfs, aes(x=values)) + 
          geom_density(aes(group=ind,  fill=ind), alpha=0.3)+scale_fill_discrete(name="")
        # Change labels on individual plot
        xi <- bquote(x[.(id)])
        plt$labels$x <- xi # bquote uses value of id for string
        plt$labels$y <- 'pdf'
      }
      return(plt)
    }
    
    # Matrix of lists, lapply gets data frame for each parameter in a list entry
    dfset <- lapply(1:numRV, function(x) setdataframe(samples[,x], sampSny[,x], nSamp, x))
    
    # Plot all the combined densities on a single pane and save
    quartz()
    grid.arrange(grobs=dfset, nrow = ceiling(numRV/2), ncol = ceiling(numRV/2))
    # marrange allows plots to be in form for saving with a title
    ml <- marrangeGrob(dfset, nrow=2, ncol=2, top='Smoothed marginal posteriors from MCMC and Snyder') 
    ggsave('sampDensity.pdf', ml, path = pathf)
  }
}

# Other way to create a data frame
# Create a data frame with labelling, each controls repetition len
#dat <- data.frame(dens = c(s1, s2), lines = rep(c("mcmc", "snyder"), 
#                                                each = length(s1)))
# Create a matrix with s1 and s2