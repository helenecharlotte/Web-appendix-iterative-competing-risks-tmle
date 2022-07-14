### follic.simulation.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Jul 14 2022 (11:53) 
## Version: 
## Last-Updated: Jul 14 2022 (11:58) 
##           By: Helene
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# NB: Note that ordering of code below is quite important

######################################################################

setwd("~/research/TMLE-from-2020june/survival-baseline/Web-appendix-iterative-competing-risks-tmle")

######################################################################

library(data.table)
library(zoo)
library(glmnet)
library(survival)
library(stringr)
library(ggplot2) 
library(ltmle)
library(nleqslv)
library(parallel)
library(foreach)
library(doParallel)
library(prodlim)
library(gridExtra)
library(survival)
library(riskRegression)
library(Matrix)
library(hdnom)
library(MASS)
library(xtable)
library(timereg)
library(cmprsk)
library(randomForestSRC)

######################################################################

source("./R/sim.data.continuous.R") 
source("./R/contmle.R") 
source("./R/cox.loss.fun.R") 
source("./R/lebesgue.loss.fun.R")
source("./R/cv.fun.R")     
source("./R/basis.fun.R")
source("./R/hal.screening.R")
source("./R/fit.hal.R")   
source("./R/cox.sl.R")  
source("./R/fit.categorical.R")
source("./R/predict.catfit.R")

######################################################################

data(follic, package="randomForestSRC")
follic <- data.table(follic) 

#-- convert variables:
follic[, stage:=as.numeric(clinstg==2)] 
follic[, chemo:=as.numeric(ch=="Y")]
follic <- follic[, -c("clinstg", "ch"), with=FALSE]

######################################################################

source("./simulation-functions/estimate.weibulls.R")
source("./simulation-functions/follic.simulation.functions.R") 

######################################################################

if (FALSE) { 

    #--- test covariate dependent censoring:
    
    seed <- 291 
    sim.follic.3 <- simulate.follic.3(observed.covars = TRUE,
                                      sim.sample = nrow(follic),
                                      counterfactual = NULL,
                                      seed = seed, 
                                      informative.censoring = TRUE,
                                      keep.times = FALSE)

    sim.follic.3[, table(status)]
    follic[, table(status)]

    sim.follic.3[, summary(time)]  
    follic[, summary(time)]

    #--- test with independent censoring:

    seed <- 291 
    sim.follic.3 <- simulate.follic.3(observed.covars = TRUE,
                                      sim.sample = nrow(follic),
                                      counterfactual = NULL,
                                      seed = seed, 
                                      informative.censoring = FALSE,
                                      keep.times = FALSE)

    sim.follic.3[, table(status)]
    follic[, table(status)]

    sim.follic.3[, summary(time)]  
    follic[, summary(time)]

}

######################################################################
### follic.simulation.R ends here
