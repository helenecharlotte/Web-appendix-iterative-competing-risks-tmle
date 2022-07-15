### follic.simulation.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Jul 14 2022 (11:53) 
## Version: 
## Last-Updated: Jul 15 2022 (11:32) 
##           By: Helene
##     Update #: 75
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

#-- define new variables:
follic[, agec := age-mean(age)]
follic[, hgbc := hgb-mean(hgb)]

######################################################################

source("./simulation/estimate.weibulls.R")
source("./simulation/follic.simulation.functions.R")
source("./simulation/follic.run.fun.R")
source("./simulation/follic.output.fun.R") 

######################################################################
#
#--- get true values of parameters: 

for (tau in c(0.25, 0.5, 0.75, 1:10)) run.follic(get.truth = TRUE, tau = tau)
#for (tau in c(0.25, 0.5, 0.75, 1:10)) run.follic(get.truth = TRUE, parameter = "1", tau = tau)

run.follic(get.truth = TRUE)
run.follic(get.truth = TRUE, parameter = "1")
run.follic(get.truth = TRUE, parameter = "0")

run.follic(get.truth = TRUE, observed.covars = FALSE)
run.follic(get.truth = TRUE, observed.covars = FALSE, parameter = "1")
run.follic(get.truth = TRUE, observed.covars = FALSE, parameter = "0")

######################################################################

no.cores <- detectCores() - 1

run.follic(M = 499, verbose = TRUE, fit.initial = "cox", no.cores = no.cores)
run.follic(M = 499, verbose = TRUE, fit.initial = "cox", no.cores = no.cores,
           observed.covars = FALSE, informative.censoring = FALSE)


run.follic(M = 500, verbose = TRUE, fit.initial = "cox", no.cores = no.cores)
run.follic(M = 500, verbose = TRUE, fit.initial = "cox", no.cores = no.cores,
           observed.covars = FALSE, informative.censoring = FALSE)

run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
           observed.covars = FALSE, informative.censoring = TRUE, sim.sample = 1000)
run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 1,
           observed.covars = FALSE, informative.censoring = TRUE, sim.sample = 1000)
run.follic(M = 500, verbose = TRUE, fit.initial = "cox", no.cores = no.cores,
           observed.covars = FALSE, informative.censoring = FALSE, sim.sample = 1000)
run.follic(M = 500, verbose = TRUE, fit.initial = "cox", no.cores = no.cores,
           observed.covars = FALSE, informative.censoring = TRUE, sim.sample = 1000)

######################################################################

cox.output <- follic.output.fun(M = 499,
                                fit.initial = "cox",
                                informative.censoring = TRUE,
                                observed.covars = TRUE,
                                sim.sample = nrow(follic))

cox.output <- follic.output.fun(M = 500,
                                fit.initial = "cox",
                                informative.censoring = FALSE,
                                observed.covars = FALSE,
                                sim.sample = nrow(follic))

rf.output <- follic.output.fun(M = 500,
                               fit.initial = "rf",
                               informative.censoring = TRUE,
                               observed.covars = FALSE,
                               sim.sample = 1000)

hal.output <- follic.output.fun(M = 500,
                               fit.initial = "hal",
                               informative.censoring = TRUE,
                               observed.covars = FALSE,
                               sim.sample = 1000)

cox.output <- follic.output.fun(M = 500,
                               fit.initial = "cox",
                               informative.censoring = TRUE,
                               observed.covars = FALSE,
                               sim.sample = 1000)

######################################################################

if (FALSE) { 

    #--- test covariate dependent censoring:
    
    seed <- sample(10000, 1) 
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

    sim.follic.3[status == 1, summary(time)]  
    follic[status == 1, summary(time)]

    #--- test with independent censoring:

    #seed <- 291 
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

    sim.follic.3[status == 1, summary(time)]  
    follic[status == 1, summary(time)]

}

######################################################################
### follic.simulation.R ends here
