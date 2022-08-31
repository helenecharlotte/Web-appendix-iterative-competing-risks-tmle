### follic.simulation.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Jul 14 2022 (11:53) 
## Version: 
## Last-Updated: Aug 24 2022 (19:24) 
##           By: Helene
##     Update #: 683
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
library(survtmle)
library(SuperLearner)

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

source("./simulation/estimate.weibulls.R")
source("./simulation/follic.simulation.functions.R")
source("./simulation/follic.run.fun.R")
source("./simulation/follic.output.fun.R")

######################################################################
#
#--- get true values of parameters: 

#for (tau in c(0.25, 0.5, 0.75, 1:10)) run.follic(get.truth = TRUE, tau = tau)
#for (tau in c(0.25, 0.5, 0.75, 1:10)) run.follic(get.truth = TRUE, parameter = "1", tau = tau)
#for (tau in c(0.25, 0.5, 0.75, 1:10)) run.follic(get.truth = TRUE, parameter = "0", tau = tau)

run.follic(get.truth = TRUE)
run.follic(get.truth = TRUE, tau = 5)
#run.follic(get.truth = TRUE, parameter = "1")
#run.follic(get.truth = TRUE, parameter = "0")


######################################################################

no.cores <- detectCores() - 1

run.follic(M = 500, verbose = TRUE, fit.initial = "cox", no.cores = no.cores)
run.follic(M = 500, verbose = TRUE, fit.initial = "cox", no.cores = no.cores,
           observed.covars = FALSE, informative.censoring = FALSE, sim.sample = 1000)

run.follic(M = 500, verbose = TRUE, fit.initial = "cox", no.cores = no.cores,
           informative.censoring = FALSE)

run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
           observed.covars = FALSE, informative.censoring = TRUE, sim.sample = 1000)
run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 1,
           observed.covars = FALSE, informative.censoring = TRUE, sim.sample = 1000)

#-- here
run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 6,
           observed.covars = FALSE, informative.censoring = FALSE, sim.sample = 1000)
run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 6,
           observed.covars = FALSE, informative.censoring = TRUE, sim.sample = 1000)
run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 6,
           observed.covars = TRUE, informative.censoring = TRUE)
run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 6,
           observed.covars = TRUE, informative.censoring = FALSE)

run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
           observed.covars = FALSE, informative.censoring = TRUE, sim.sample = 1000)
run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
           observed.covars = FALSE, informative.censoring = FALSE, sim.sample = 1000)

#-- here2
run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 6,
           observed.covars = TRUE, informative.censoring = TRUE, sim.sample = 1000)
run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 6,
           observed.covars = TRUE, informative.censoring = FALSE, sim.sample = 1000)
run.follic(M = 500, verbose = FALSE, fit.initial = "cox", no.cores = 6,
           observed.covars = TRUE, informative.censoring = FALSE, sim.sample = 1000)
#-- here2.5
run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
           observed.covars = TRUE, informative.censoring = TRUE, sim.sample = 1000)
run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
           observed.covars = TRUE, informative.censoring = FALSE, sim.sample = 1000)


run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
           observed.covars = TRUE, informative.censoring = FALSE, sim.sample = 1000)
run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
           observed.covars = FALSE, informative.censoring = FALSE, sim.sample = 1000)


run.follic(M = 1, verbose = FALSE, fit.initial = "rf", no.cores = 1,
           observed.covars = FALSE, informative.censoring = FALSE, sim.sample = 1000)

run.follic(M = 500, verbose = TRUE, fit.initial = "cox", no.cores = no.cores,
           observed.covars = FALSE, informative.censoring = TRUE, sim.sample = 1000)

#-- here3
run.follic(M = 6, verbose = TRUE, fit.initial = "hal", no.cores = 6,
           hal.sl = TRUE,
           observed.covars = FALSE, informative.censoring = TRUE, sim.sample = 1000)

######################################################################
#-- HERE (July, 18)

run.follic(M = 500, verbose = TRUE, fit.initial = "cox", no.cores = 6,
           observed.covars = TRUE,
           randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           sim.sample = 1000)

run.follic(M = 1, verbose = TRUE, fit.initial = "hal", no.cores = 6,
           observed.covars = TRUE,
           cut.two.way = 0,
           #browse = TRUE,
           randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           sim.sample = 1000)

run.follic(M = 1, verbose = TRUE, fit.initial = "hal", no.cores = 6,
           observed.covars = TRUE,
           hal.sl = TRUE,
           randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           sim.sample = 1000)

# new version where simulate *randomized* treatment (*independent* censoring):
if (FALSE) run.follic(M = 1, verbose = TRUE, fit.initial = "cox", no.cores = 6,
                      observed.covars = TRUE,
                      randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                      sim.sample = 1000)
if (FALSE) run.follic(M = 1, verbose = TRUE, fit.initial = "hal", no.cores = 6,
                      observed.covars = TRUE,
                      randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                      sim.sample = 1000)
if (FALSE) run.follic(M = 1, verbose = TRUE, fit.initial = "hal", no.cores = 6,
                      observed.covars = TRUE,
                      cut.time = 12, cut.time.A = 10, cut.time.covar = 10,
                      randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                      sim.sample = 1000)
if (FALSE) run.follic(M = 1, verbose = TRUE, fit.initial = "hal", no.cores = 6,
                      observed.covars = TRUE,
                      hal.sl = TRUE,
                      randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                      sim.sample = 1000)
if (FALSE) run.follic(M = 1, verbose = TRUE, fit.initial = "hal", no.cores = 6,
                      observed.covars = TRUE,
                      hal.sl = TRUE,
                      randomized.treatment = TRUE, observed.treatment = FALSE,
                      informative.censoring = TRUE,
                      sim.sample = 1000)
if (FALSE) run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
                      observed.covars = TRUE,
                      randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                      sim.sample = 1000)
if (FALSE) run.follic(M = 3, verbose = TRUE, fit.initial = "hal", no.cores = 1,
                      observed.covars = TRUE,
                      sim.sample = 1000)

run.follic(M = 1, verbose = FALSE, fit.initial = "hal", no.cores = 1,
           observed.covars = TRUE,
           #hal.sl = TRUE,
           browse = TRUE,
           #randomized.treatment = TRUE, observed.treatment = FALSE,
           informative.censoring = TRUE, #FALSE,
           sim.sample = nrow(follic))

# informative censoring: 
if (FALSE) run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 6,
                      observed.covars = TRUE, informative.censoring = TRUE, sim.sample = 1000)
if (FALSE) run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
           observed.covars = TRUE, informative.censoring = TRUE, sim.sample = 1000)
if (FALSE) run.follic(M = 500, verbose = FALSE, fit.initial = "cox", no.cores = 6,
                      observed.covars = TRUE, informative.censoring = TRUE, sim.sample = 1000)

run.follic(M = 500, verbose = TRUE, fit.initial = "hal", no.cores = 6,
           hal.sl = TRUE,
           observed.covars = TRUE, informative.censoring = TRUE, sim.sample = 1000)

# perhaps (simulated non-randomized treatment):
run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 6,
           observed.covars = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           sim.sample = 1000)
run.follic(M = 500, verbose = FALSE, fit.initial = "cox", no.cores = 6,
           observed.covars = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           sim.sample = 1000)

if (FALSE) {
    # this was never run:
    run.follic(M = 500, verbose = FALSE, fit.initial = "cox", no.cores = 6,
               observed.covars = TRUE, informative.censoring = TRUE, sim.sample = 1000)

    # need to re-run independent censoring:
    run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 6,
               observed.covars = TRUE, informative.censoring = FALSE, sim.sample = 1000)
    run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
               observed.covars = TRUE, informative.censoring = FALSE, sim.sample = 1000)
    run.follic(M = 500, verbose = FALSE, fit.initial = "cox", no.cores = 6,
               observed.covars = TRUE, informative.censoring = FALSE, sim.sample = 1000)

    # fun to have too? 
    run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 6,
               observed.covars = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
               sim.sample = 1000)

    # these should be done:
    run.follic(M = 500, verbose = FALSE, fit.initial = "rf", no.cores = 1,
               observed.covars = TRUE, informative.censoring = TRUE, sim.sample = 1000)
    run.follic(M = 500, verbose = FALSE, fit.initial = "hal", no.cores = 6,
               observed.covars = TRUE, informative.censoring = TRUE, sim.sample = 1000)
}


######################################################################
#-- trying survtmle (July, 19)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           informative.censoring = TRUE,
           fit.survtmle = TRUE,
           grid.survtmle = (0:10)*4,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           informative.censoring = TRUE,
           fit.survtmle = TRUE,
           grid.survtmle = (0:20)*2,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           informative.censoring = TRUE,
           fit.survtmle = TRUE,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           informative.censoring = TRUE,
           fit.survtmle = TRUE,
           grid.survtmle = (0:80)/2,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           informative.censoring = TRUE,
           fit.survtmle = TRUE,
           sl.survtmle = TRUE,
           grid.survtmle = (0:10)*4,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           informative.censoring = TRUE,
           fit.survtmle = TRUE,
           sl.survtmle = TRUE,
           grid.survtmle = (0:20)*2,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           informative.censoring = TRUE,
           fit.survtmle = TRUE,
           sl.survtmle = TRUE,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           informative.censoring = TRUE,
           fit.survtmle = TRUE,
           sl.survtmle = TRUE,
           grid.survtmle = (0:80)/2,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           fit.survtmle = TRUE,
           grid.survtmle = (0:10)*4,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           fit.survtmle = TRUE,
           grid.survtmle = (0:20)*2,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           fit.survtmle = TRUE,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           fit.survtmle = TRUE,
           grid.survtmle = (0:80)/2,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           informative.censoring = TRUE,
           fit.survtmle = TRUE,
           grid.survtmle = (0:80)/2,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           fit.survtmle = TRUE,
           sl.survtmle = TRUE,
           grid.survtmle = (0:10)*4,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           fit.survtmle = TRUE,
           sl.survtmle = TRUE,
           grid.survtmle = (0:20)*2,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           fit.survtmle = TRUE,
           sl.survtmle = TRUE,
           sim.sample = 1000)

run.follic(M = 500, verbose = FALSE, no.cores = 6,
           observed.covars = TRUE,
           randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
           fit.survtmle = TRUE,
           sl.survtmle = TRUE,
           grid.survtmle = (0:80)/2,
           sim.sample = 1000)


######################################################################

cox.output <- follic.output.fun(M = 500,
                                fit.initial = "cox",
                                informative.censoring = TRUE,
                                observed.covars = TRUE,
                                sim.sample = nrow(follic))

cox.output <- follic.output.fun(M = 500,
                                fit.initial = "cox",
                                informative.censoring = FALSE,
                                observed.covars = FALSE,
                                sim.sample = 1000)

hal.output <- follic.output.fun(M = 500,
                               fit.initial = "hal",
                               informative.censoring = FALSE,
                               observed.covars = FALSE,
                               sim.sample = 1000)

rf.output <- follic.output.fun(M = 500,
                               fit.initial = "rf",
                               informative.censoring = FALSE,
                               observed.covars = FALSE,
                               sim.sample = 1000)

hal.output <- follic.output.fun(M = 500,
                               fit.initial = "hal",
                               informative.censoring = TRUE,
                               observed.covars = FALSE,
                               sim.sample = 1000)

rf.output <- follic.output.fun(M = 500,
                               fit.initial = "rf",
                               informative.censoring = TRUE,
                               observed.covars = FALSE,
                               sim.sample = 1000)

hal.output <- follic.output.fun(M = 500,
                               fit.initial = "hal",
                               informative.censoring = TRUE,
                               observed.covars = TRUE,
                               sim.sample = nrow(follic))

hal.output <- follic.output.fun(M = 500,
                               fit.initial = "hal",
                               informative.censoring = FALSE,
                               observed.covars = TRUE,
                               sim.sample = nrow(follic))

hal.output <- follic.output.fun(M = 500,
                               fit.initial = "hal",
                               informative.censoring = TRUE,
                               observed.covars = TRUE,
                               sim.sample = 1000)

hal.output <- follic.output.fun(M = 500,
                               fit.initial = "hal",
                               informative.censoring = FALSE,
                               observed.covars = TRUE,
                               sim.sample = 1000)

cox.output <- follic.output.fun(M = 500,
                               fit.initial = "cox",
                               informative.censoring = FALSE,
                               observed.covars = TRUE,
                               sim.sample = 1000)

cox.output <- follic.output.fun(M = 500,
                               fit.initial = "cox",
                               informative.censoring = TRUE,
                               observed.covars = FALSE,
                               sim.sample = 1000)

#-- (SIMULATED covariates)

rf.output <- follic.output.fun(M = 500,
                               fit.initial = "rf",
                               informative.censoring = TRUE,
                               observed.covars = FALSE,
                               sim.sample = 1000)

rf.output <- follic.output.fun(M = 500,
                               fit.initial = "rf",
                               informative.censoring = FALSE,
                               observed.covars = FALSE,
                               sim.sample = 1000)

hal.output <- follic.output.fun(M = 500,
                                fit.initial = "hal",
                                informative.censoring = TRUE,
                                observed.covars = FALSE,
                                sim.sample = 1000)

hal.output <- follic.output.fun(M = 500,
                                fit.initial = "hal",
                                informative.censoring = FALSE,
                                observed.covars = FALSE,
                                sim.sample = 1000)


#####################################################################

#-- HERE? (OBSERVED covariates)

rf.inf.output <- follic.output.fun(M = 500,
                                   fit.initial = "rf",
                                   informative.censoring = TRUE,
                                   observed.covars = TRUE,
                                   output.directory = "backup-july-19/simulation",
                                   sim.sample = 1000)

rf.inf.output <- follic.output.fun(M = 500,
                                   fit.initial = "rf",
                                   informative.censoring = TRUE,
                                   observed.covars = TRUE,
                                   sim.sample = 1000)

rf.ind.output <- follic.output.fun(M = 500,
                                   fit.initial = "rf",
                                   informative.censoring = FALSE,
                                   observed.covars = TRUE,
                                   sim.sample = 1000)

rf.rand.output <- follic.output.fun(M = 500,
                                    fit.initial = "rf",
                                    informative.censoring = FALSE,
                                    observed.covars = TRUE,
                                    randomized.treatment = TRUE,
                                    observed.treatment = FALSE,
                                    output.directory = "backup-july-19/simulation",
                                    sim.sample = 1000)

rf.rand.output <- follic.output.fun(M = 500,
                                    fit.initial = "rf",
                                    informative.censoring = FALSE,
                                    observed.covars = TRUE,
                                    randomized.treatment = TRUE,
                                    observed.treatment = FALSE,
                                    sim.sample = 1000)

hal.inf.output <- follic.output.fun(M = 500,
                                    fit.initial = "hal",
                                    informative.censoring = TRUE,
                                    observed.covars = TRUE,
                                    output.directory = "backup-july-18/simulation",
                                    sim.sample = 1000)

hal.inf.output <- follic.output.fun(M = 500,
                                    fit.initial = "hal",
                                    informative.censoring = TRUE,
                                    observed.covars = TRUE,
                                    sim.sample = 1000)

hal.sl.inf.output <- follic.output.fun(M = 500,
                                       fit.initial = "hal",
                                       hal.sl = TRUE,
                                       informative.censoring = TRUE,
                                       observed.covars = TRUE,
                                       output.directory = "backup-july-27/simulation",
                                       sim.sample = 1000)

hal.sl.inf.output <- follic.output.fun(M = 500,
                                       fit.initial = "hal",
                                       hal.sl = TRUE,
                                       informative.censoring = TRUE,
                                       observed.covars = TRUE,
                                       sim.sample = 1000)

hal.ind.output <- follic.output.fun(M = 1000,#500,
                                    fit.initial = "hal",
                                    informative.censoring = FALSE,
                                    observed.covars = TRUE,
                                    sim.sample = 1000)

rf.inf.output.n2000 <- follic.output.fun(M = 500,#1000,
                                         fit.initial = "rf",
                                         onlyM = 500,
                                         informative.censoring = TRUE,
                                         observed.covars = TRUE,
                                         sim.sample = 2000)

hal.inf.output.n2000 <- follic.output.fun(M = 500,#1000,
                                          fit.initial = "hal",
                                          # onlyM = 500,
                                          informative.censoring = TRUE,
                                          observed.covars = TRUE,
                                          sim.sample = 2000)

hal.sl.inf.output.n2000 <- follic.output.fun(M = 500,
                                             fit.initial = "hal",
                                             hal.sl = TRUE,
                                             #browse = TRUE,
                                             se.bound = 0.037,
                                             informative.censoring = TRUE,
                                             observed.covars = TRUE,
                                             sim.sample = 2000)

cox.inf.output.n2000 <- follic.output.fun(M = 500,#1000,
                                          onlyM = 500,
                                          fit.initial = "cox",
                                          informative.censoring = TRUE,
                                          observed.covars = TRUE,
                                          sim.sample = 2000)

rf.rand.output.n2000 <- follic.output.fun(M = 500,
                                          fit.initial = "rf",
                                          # onlyM = 500,
                                          informative.censoring = FALSE,
                                          observed.covars = TRUE,
                                          randomized.treatment = TRUE,
                                          observed.treatment = FALSE,
                                          #output.directory = "backup-aug-1/simulation",
                                          sim.sample = 2000)

hal.rand.output.n2000 <- follic.output.fun(M = 1000,
                                           fit.initial = "hal",
                                           onlyM = 500,
                                           informative.censoring = FALSE,
                                           observed.covars = TRUE,
                                           randomized.treatment = TRUE,
                                           observed.treatment = FALSE,
                                           output.directory = "backup-aug-1/simulation",
                                           sim.sample = 2000)

hal.rand.output.n2000 <- follic.output.fun(M = 500,#1000,
                                           fit.initial = "hal",
                                           onlyM = 500,
                                           informative.censoring = FALSE,
                                           observed.covars = TRUE,
                                           randomized.treatment = TRUE,
                                           observed.treatment = FALSE,
                                           #output.directory = "backup-aug-1/simulation",
                                           sim.sample = 2000)

hal.sl.rand.output.n2000 <- follic.output.fun(M = 500,
                                              fit.initial = "hal",
                                              hal.sl = TRUE,
                                              informative.censoring = FALSE,
                                              observed.covars = TRUE,
                                              #browse = TRUE,
                                              randomized.treatment = TRUE,
                                              observed.treatment = FALSE,
                                              #output.directory = "backup-aug-1/simulation",
                                              sim.sample = 2000)

cox.rand.output.n2000 <- follic.output.fun(M = 1000,
                                           onlyM = 500,
                                           fit.initial = "cox",
                                           informative.censoring = FALSE,
                                           observed.covars = TRUE,
                                           randomized.treatment = TRUE,
                                           observed.treatment = FALSE,
                                           output.directory = "backup-aug-1/simulation",
                                           sim.sample = 2000)

cox.rand.output.n2000 <- follic.output.fun(M = 500,#1000,
                                           onlyM = 500,
                                           fit.initial = "cox",
                                           informative.censoring = FALSE,
                                           observed.covars = TRUE,
                                           randomized.treatment = TRUE,
                                           observed.treatment = FALSE,
                                           #output.directory = "backup-aug-1/simulation",
                                           sim.sample = 2000)

hal.inf.output.nfollic <- follic.output.fun(M = 500,
                                            fit.initial = "hal",
                                            informative.censoring = TRUE,
                                            observed.covars = TRUE,
                                            #output.directory = "backup-aug-1/simulation",
                                            sim.sample = nrow(follic))

hal.sl.inf.output.nfollic <- follic.output.fun(M = 500,
                                             fit.initial = "hal",
                                             hal.sl = TRUE,
                                             informative.censoring = TRUE,
                                             observed.covars = TRUE,
                                             sim.sample = nrow(follic))

cox.inf.output.nfollic <- follic.output.fun(M = 500,
                                            fit.initial = "cox",
                                            informative.censoring = TRUE,
                                            observed.covars = TRUE,
                                            #output.directory = "backup-aug-1/simulation",
                                            sim.sample = nrow(follic))


hal.rand.output.nfollic <- follic.output.fun(M = 500,
                                             fit.initial = "hal",
                                             informative.censoring = FALSE,
                                             observed.covars = TRUE,
                                             randomized.treatment = TRUE,
                                             observed.treatment = FALSE,
                                            # output.directory = "backup-aug-1/simulation",
                                             sim.sample = nrow(follic))

hal.sl.rand.output.nfollic <- follic.output.fun(M = 500,
                                                fit.initial = "hal",
                                                hal.sl = TRUE,
                                                informative.censoring = FALSE,
                                                observed.covars = TRUE,
                                                #browse = TRUE,
                                                randomized.treatment = TRUE,
                                                observed.treatment = FALSE,
                                                output.directory = "backup-aug-1/simulation",
                                                sim.sample = nrow(follic))

cox.rand.output.nfollic <- follic.output.fun(M = 500,
                                             fit.initial = "cox",
                                             informative.censoring = FALSE,
                                             observed.covars = TRUE,
                                             randomized.treatment = TRUE,
                                             observed.treatment = FALSE,
                                             #output.directory = "backup-aug-1/simulation",
                                             sim.sample = nrow(follic))

hal.rand.output <- follic.output.fun(M = 1000,
                                     fit.initial = "hal",
                                     informative.censoring = FALSE,
                                     observed.covars = TRUE,
                                     randomized.treatment = TRUE,
                                     observed.treatment = FALSE,
                                     #output.directory = "backup-july-19/simulation",
                                     sim.sample = 1000)

hal.rand.output <- follic.output.fun(M = 500,
                                     fit.initial = "hal",
                                     informative.censoring = FALSE,
                                     observed.covars = TRUE,
                                     randomized.treatment = TRUE,
                                     observed.treatment = FALSE,
                                     sim.sample = 1000)

hal.sl.rand.output <- follic.output.fun(M = 500,
                                        fit.initial = "hal",
                                        hal.sl = TRUE,
                                        informative.censoring = FALSE,
                                        observed.covars = TRUE,
                                        #browse = TRUE,
                                        randomized.treatment = TRUE,
                                        observed.treatment = FALSE,
                                        #output.directory = "backup-july-27/simulation",
                                        sim.sample = 1000)

hal.sim.output <- follic.output.fun(M = 500,
                                    fit.initial = "hal",
                                    informative.censoring = FALSE,
                                    observed.covars = TRUE,
                                    observed.treatment = FALSE,
                                    output.directory = "backup-july-19/simulation",
                                    sim.sample = 1000)

cox.inf.output <- follic.output.fun(M = 500,
                                    fit.initial = "cox",
                                    informative.censoring = TRUE,
                                    observed.covars = TRUE,
                                    output.directory = "backup-july-19/simulation",
                                    sim.sample = 1000)

cox.inf.output <- follic.output.fun(M = 500,
                                    fit.initial = "cox",
                                    informative.censoring = TRUE,
                                    observed.covars = TRUE,
                                    sim.sample = 1000)

cox.ind.output <- follic.output.fun(M = 500,
                                    fit.initial = "cox",
                                    informative.censoring = FALSE,
                                    observed.covars = TRUE,
                                    sim.sample = 1000)


cox.rand.output <- follic.output.fun(M = 500,
                                     fit.initial = "cox",
                                     informative.censoring = FALSE,
                                     observed.covars = TRUE,
                                     randomized.treatment = TRUE,
                                     observed.treatment = FALSE,
                                     #output.directory = "backup-july-19/simulation",
                                     sim.sample = 1000)

cox.sim.output <- follic.output.fun(M = 500,
                                    fit.initial = "cox",
                                    informative.censoring = FALSE,
                                    observed.covars = TRUE,
                                    observed.treatment = FALSE,
                                    output.directory = "backup-july-19/simulation",
                                    sim.sample = 1000)


cox.rand.output.5 <- follic.output.fun(M = 500,
                                       tau = 5,
                                       fit.initial = "cox",
                                       informative.censoring = FALSE,
                                       observed.covars = TRUE,
                                       randomized.treatment = TRUE,
                                       observed.treatment = FALSE,
                                       sim.sample = 1000)

hal.rand.output.5 <- follic.output.fun(M = 1000,#500,
                                       tau = 5,
                                       fit.initial = "hal",
                                       informative.censoring = FALSE,
                                       observed.covars = TRUE,
                                       randomized.treatment = TRUE,
                                       observed.treatment = FALSE,
                                       sim.sample = 1000)

cox.inf.output.5 <- follic.output.fun(M = 500,
                                      tau = 5,
                                      fit.initial = "cox",
                                      informative.censoring = TRUE,
                                      observed.covars = TRUE,
                                      sim.sample = 1000)

hal.inf.output.5 <- follic.output.fun(M = 1000,#500,
                                      tau = 5,
                                      fit.initial = "hal",
                                      informative.censoring = TRUE,
                                      observed.covars = TRUE,
                                      sim.sample = 1000)




#--- informative censoring:
follic.compare.results(hal.inf.output, rf.inf.output, cox.inf.output)

#--- randomized treatment + independent censoring:
follic.compare.results(hal.rand.output, rf.rand.output, cox.rand.output)

#--- independent censoring:
follic.compare.results(hal.ind.output, rf.ind.output, cox.ind.output)

#--- simulated treatment + independent censoring (OBS! Not run for RF):
follic.compare.results(hal.sim.output, cox.sim.output, cox.sim.output)

#--- compare to survtmle!

survtmle.rand.grid0 <- follic.output.survtmle(M = 500,
                                              randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                              observed.covars = TRUE,
                                              sl.survtmle = TRUE,
                                              grid.survtmle = (0:10)*4, 
                                              sim.sample = 1000)

survtmle.rand.grid1 <- follic.output.survtmle(M = 500,
                                              randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                              observed.covars = TRUE,
                                              grid.survtmle = (0:20)*2, 
                                              sim.sample = 1000)

survtmle.rand.grid2 <- follic.output.survtmle(M = 500,
                                              randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                              observed.covars = TRUE,
                                              grid.survtmle = (0:40), 
                                              sim.sample = 1000)

survtmle.rand.grid3 <- follic.output.survtmle(M = 500,
                                              randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                              observed.covars = TRUE,
                                              grid.survtmle = (0:80)/2, 
                                              sim.sample = 1000)

survtmle.inf.grid0 <- follic.output.survtmle(M = 500,
                                             informative.censoring = TRUE,
                                             observed.covars = TRUE,
                                             grid.survtmle = (0:10)*4, 
                                             sim.sample = 1000)

survtmle.inf.grid1 <- follic.output.survtmle(M = 500,
                                             informative.censoring = TRUE,
                                             observed.covars = TRUE,
                                             grid.survtmle = (0:20)*2, 
                                             sim.sample = 1000)

survtmle.inf.grid2 <- follic.output.survtmle(M = 500,
                                             informative.censoring = TRUE,
                                             observed.covars = TRUE,
                                             grid.survtmle = (0:40), 
                                             sim.sample = 1000)

survtmle.inf.grid3 <- follic.output.survtmle(M = 500,
                                             informative.censoring = TRUE,
                                             observed.covars = TRUE,
                                             grid.survtmle = (0:80)/2, 
                                             sim.sample = 1000)

# informative censoring: 
cbind(unlist(survtmle.inf.grid0),unlist(survtmle.inf.grid1), unlist(survtmle.inf.grid2), unlist(survtmle.inf.grid3),
      follic.compare.results(hal.inf.output,hal.inf.output,hal.inf.output)[,1])

# randomized treatment + independent censoring: 
cbind(unlist(survtmle.rand.grid0),unlist(survtmle.rand.grid1), unlist(survtmle.rand.grid2), unlist(survtmle.rand.grid3),
      follic.compare.results(hal.rand.output,hal.rand.output,hal.rand.output)[,1])

######################################################################
# plotting results

contmle.results <- do.call("rbind", lapply(list("Randomized treatment + independent censoring", "Informative censoring"), function(outer.which) {
    out.inner <- data.table(do.call("rbind", lapply(list("hal", "cox", "rf"), function(inner.initial) {
        try(unlist(follic.output.fun(M = 500,
                                     randomized.treatment = outer.which == "Randomized treatment + independent censoring",
                                     observed.treatment = outer.which != "Randomized treatment + independent censoring",
                                     informative.censoring = outer.which != "Randomized treatment + independent censoring",
                                     observed.covars = TRUE,
                                     fit.initial = inner.initial,
                                     sim.sample = 1000)))
    })))
    out.inner[, initial := c("HAL+TMLE", "Cox+TMLE", "RF+TMLE")][, setting := outer.which]
    return(out.inner[substr(bias.tmle, 1, 5) != "Error"])
}))

contmle.results <- rbind(contmle.results, do.call("rbind", lapply(unique(contmle.results[["setting"]]), function(which) {
    contmle.results[setting == which & tolower(substr(initial,1,3)) == "hal"][, initial := "KM"][, bias.tmle := bias.km][, se.tmle := se.km][, sd.tmle := sd.km][, mse.tmle := mse.km]
})))

grid.arrange(ggplot(contmle.results) + theme_bw() +
             geom_point(aes(x = initial, y = as.numeric(bias.tmle)), shape = 4) +
             facet_grid(. ~ setting) +
             geom_hline(yintercept = 0, linetype = "dashed", col = "red", alpha = 0.5) +
             ylab("bias") + xlab(""),
             ggplot(contmle.results) + theme_bw() +
             geom_point(aes(x = initial, y = as.numeric(cov.tmle)), shape = 4) +
             facet_grid(. ~ setting) +
             geom_hline(yintercept = 0.95, linetype = "dashed", col = "red", alpha = 0.5) +
             ylab("coverage") + xlab(""),
             ggplot(contmle.results) + theme_bw() +
             geom_point(aes(x = initial, y = as.numeric(sd.tmle)/as.numeric(se.tmle)), shape = 4) +
             facet_grid(. ~ setting) +
             geom_hline(yintercept = 1, linetype = "dashed", col = "red", alpha = 0.5) +
             ylab("SD/SE") + xlab(""),
             ggplot(contmle.results) + theme_bw() +
             geom_point(aes(x = initial, y = as.numeric(mse.tmle)/as.numeric(mse.km)), shape = 4) +
             facet_grid(. ~ setting) +
             geom_hline(yintercept = 1, linetype = "dashed", col = "red", alpha = 0.5) +
             ylab("MSE/MSE(KM)") + xlab(""))


survtmle.results <- do.call("rbind", lapply(list("Randomized treatment + independent censoring", "Informative censoring"), function(outer.which) {
    do.call("rbind", lapply(list((0:10)*4, (0:20)*2, (0:40), (0:80)/2), function(inner.grid) {
        out.inner <- data.table(do.call("rbind", lapply(list("glm", "sl"), function(inner.sl) {
            try(unlist(follic.output.survtmle(M = 500,
                                       randomized.treatment = outer.which == "Randomized treatment + independent censoring",
                                       observed.treatment = outer.which != "Randomized treatment + independent censoring",
                                       informative.censoring = outer.which != "Randomized treatment + independent censoring",
                                       observed.covars = TRUE,
                                       sl.survtmle = inner.sl == "sl",
                                       grid.survtmle = inner.grid, 
                                       sim.sample = 1000)))
        })))
        out.inner[, initial := c("GLM", "SL")][, grid.length.used := length(inner.grid)][, setting := outer.which]
        return(out.inner[substr(bias.survtmle, 1, 5) != "Error"])
    }))
}))

survtmle.results.fixed.grid <- survtmle.results[grid.length.used == 41]
survtmle.results.fixed.grid[, initial := paste0("survtmle (", initial, ")")]
names(survtmle.results.fixed.grid) <- gsub("surv", "", names(survtmle.results.fixed.grid))
tmle.results <- rbind(contmle.results, survtmle.results.fixed.grid, fill = TRUE)
tmle.results[, mse.km := na.locf(mse.km), by = setting]

grid.arrange(ggplot(tmle.results) + theme_bw() +
             geom_point(aes(x = initial, y = as.numeric(bias.tmle)), shape = 4) +
             facet_grid(. ~ setting) +
             geom_hline(yintercept = 0, linetype = "dashed", col = "red", alpha = 0.5) +
             ylab("bias") + xlab("") +
             theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                   strip.text = element_text(size = 14),
                   strip.background = element_blank()),
             ggplot(tmle.results) + theme_bw() +
             geom_point(aes(x = initial, y = as.numeric(cov.tmle)), shape = 4) +
             facet_grid(. ~ setting) +
             geom_hline(yintercept = 0.95, linetype = "dashed", col = "red", alpha = 0.5) +
             ylab("coverage") + xlab("") +
             theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                   strip.text = element_text(size = 14),
                   strip.background = element_blank()),
             ggplot(tmle.results) + theme_bw() +
             geom_point(aes(x = initial, y = as.numeric(sd.tmle)/as.numeric(se.tmle)), shape = 4) +
             facet_grid(. ~ setting) +
             geom_hline(yintercept = 1, linetype = "dashed", col = "red", alpha = 0.5) +
             ylab("SD/SE") + xlab("") +
             theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                   strip.text = element_text(size = 14),
                   strip.background = element_blank()),
             ggplot(tmle.results) + theme_bw() +
             geom_point(aes(x = initial, y = as.numeric(mse.tmle)/as.numeric(mse.km)), shape = 4) +
             facet_grid(. ~ setting) +
             geom_hline(yintercept = 1, linetype = "dashed", col = "red", alpha = 0.5) +
             ylab("MSE/MSE(KM)") + xlab("") +
             theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                   strip.text = element_text(size = 14),
                   strip.background = element_blank()))


grid.arrange(ggplot(survtmle.results) + theme_bw() +
             geom_point(aes(x = grid.length.used, y = as.numeric(bias.survtmle)), shape = 4) +
             facet_grid(initial ~ setting) +
             geom_hline(yintercept = 0, linetype = "dashed", col = "red", alpha = 0.5) +
             ylab("bias") + xlab("grid length used") +
             theme(strip.text = element_text(size = 14),
                   strip.background = element_blank()),
             ggplot(survtmle.results) + theme_bw() +
             geom_point(aes(x = grid.length.used, y = as.numeric(cov.survtmle)), shape = 4) +
             facet_grid(initial ~ setting) +
             geom_hline(yintercept = 0.95, linetype = "dashed", col = "red", alpha = 0.5) +
             ylab("coverage") + xlab("grid length used") +
             theme(strip.text = element_text(size = 14),
                   strip.background = element_blank()),
             ggplot(survtmle.results) + theme_bw() +
             geom_point(aes(x = grid.length.used, y = as.numeric(sd.survtmle)/as.numeric(se.survtmle)), shape = 4) +
             facet_grid(initial ~ setting) +
             geom_hline(yintercept = 1, linetype = "dashed", col = "red", alpha = 0.5) +
             ylab("SD/SE") + xlab("grid length used") +
             theme(strip.text = element_text(size = 14),
                   strip.background = element_blank()))


######################################################################
# collecting again ? (aug 2, 2022)

hal.inf.output.n1000 <- follic.output.fun(M = 500,#1000,#500,
                                           fit.initial = "hal",
                                           #onlyM = 500,
                                           informative.censoring = TRUE,
                                           observed.covars = TRUE,
                                           sim.sample = 1000)

hal2.inf.output.n1000 <- follic.output.fun(M = 500,#1000,#500,
                                           fit.initial = "hal",
                                           #onlyM = 500,
                                           informative.censoring = TRUE,
                                           observed.covars = TRUE,
                                           cut.time = 12, cut.time.A = 10, cut.time.covar = 10,
                                           sim.sample = 1000)

hal2.rand.output.n1000 <- follic.output.fun(M = 500,#1000,
                                            fit.initial = "hal",
                                            #onlyM = 500,
                                            informative.censoring = FALSE,
                                            observed.covars = TRUE,
                                            randomized.treatment = TRUE,
                                            observed.treatment = FALSE,
                                            cut.time = 12, cut.time.A = 10, cut.time.covar = 10,
                                            sim.sample = 1000)

hal3.rand.output.n1000 <- follic.output.fun(M = 500,#1000,
                                            fit.initial = "hal",
                                            #onlyM = 500,
                                            informative.censoring = FALSE,
                                            observed.covars = TRUE,
                                            randomized.treatment = TRUE,
                                            observed.treatment = FALSE,
                                            cut.time = 25,
                                            sim.sample = 1000)

hal4.rand.output.n1000 <- follic.output.fun(M = 500,#1000,
                                            fit.initial = "hal",
                                            #onlyM = 500,
                                            informative.censoring = FALSE,
                                            observed.covars = TRUE,
                                            randomized.treatment = TRUE,
                                            observed.treatment = FALSE,
                                            cut.two.way = 15,
                                            sim.sample = 1000)

hal.rand.output.n1000 <- follic.output.fun(M = 500,#1000,
                                            fit.initial = "hal",
                                            #onlyM = 500,
                                            informative.censoring = FALSE,
                                            observed.covars = TRUE,
                                            randomized.treatment = TRUE,
                                            observed.treatment = FALSE,
                                            sim.sample = 1000)



hal.rand.output.n2000 <- follic.output.fun(M = 500,#1000,
                                           fit.initial = "hal",
                                           #onlyM = 500,
                                           informative.censoring = FALSE,
                                           observed.covars = TRUE,
                                           randomized.treatment = TRUE,
                                           observed.treatment = FALSE,
                                           #output.directory = "backup-aug-1/simulation",
                                           sim.sample = 2000)

rf.rand.output.n2000 <- follic.output.fun(M = 500,
                                          fit.initial = "rf",
                                          # onlyM = 500,
                                          informative.censoring = FALSE,
                                          observed.covars = TRUE,
                                          randomized.treatment = TRUE,
                                          observed.treatment = FALSE,
                                          #output.directory = "backup-aug-1/simulation",
                                          sim.sample = 2000)

hal.sl.rand.output.n2000 <- follic.output.fun(M = 500,
                                              fit.initial = "hal",
                                              hal.sl = TRUE,
                                              informative.censoring = FALSE,
                                              observed.covars = TRUE,
                                              #browse = TRUE,
                                              randomized.treatment = TRUE,
                                              observed.treatment = FALSE,
                                              #output.directory = "backup-aug-1/simulation",
                                              sim.sample = 2000)

hal.rand.output.nfollic <- follic.output.fun(M = 500,
                                             fit.initial = "hal",
                                             informative.censoring = FALSE,
                                             observed.covars = TRUE,
                                             randomized.treatment = TRUE,
                                             observed.treatment = FALSE,
                                             #output.directory = "backup-aug-1/simulation",
                                             sim.sample = nrow(follic))

hal.sl.rand.output.nfollic <- follic.output.fun(M = 500,
                                                fit.initial = "hal",
                                                hal.sl = TRUE,
                                                informative.censoring = FALSE,
                                                observed.covars = TRUE,
                                                #browse = TRUE,
                                                randomized.treatment = TRUE,
                                                observed.treatment = FALSE,
                                                #output.directory = "backup-aug-1/simulation",
                                                sim.sample = nrow(follic))

hal.inf.output.n2000 <- follic.output.fun(M = 500,#1000,
                                          fit.initial = "hal",
                                          onlyM = 500,
                                          informative.censoring = TRUE,
                                          observed.covars = TRUE,
                                          #output.directory = "backup-aug-1/simulation",
                                          sim.sample = 2000)

hal.sl.inf.output.n2000 <- follic.output.fun(M = 500,
                                             fit.initial = "hal",
                                             hal.sl = TRUE,
                                             informative.censoring = TRUE,
                                             observed.covars = TRUE,
                                             #output.directory = "backup-aug-1/simulation",
                                             sim.sample = 2000)

rf.inf.output.n2000 <- follic.output.fun(M = 500,
                                         fit.initial = "rf",
                                         onlyM = 500,
                                         informative.censoring = TRUE,
                                         observed.covars = TRUE,
                                         #output.directory = "backup-aug-1/simulation",
                                         sim.sample = 2000)

hal.inf.output.n1000 <- follic.output.fun(M = 500,#1000,#500,
                                          fit.initial = "hal",
                                          #onlyM = 500,
                                          informative.censoring = TRUE,
                                          observed.covars = TRUE,
                                          #output.directory = "backup-aug-1/simulation",
                                          sim.sample = 1000)

cox.inf.output.n1000 <- follic.output.fun(M = 500,
                                          fit.initial = "cox",
                                          #onlyM = 500,
                                          informative.censoring = TRUE,
                                          observed.covars = TRUE,
                                          #output.directory = "backup-aug-1/simulation",
                                          sim.sample = 1000)

hal.sl.inf.output.n1000 <- follic.output.fun(M = 500,
                                             fit.initial = "hal",
                                             hal.sl = TRUE,
                                             informative.censoring = TRUE,
                                             observed.covars = TRUE,
                                             #output.directory = "backup-aug-1/simulation",
                                             sim.sample = 1000)

rf.inf.output.n1000 <- follic.output.fun(M = 500,
                                         fit.initial = "rf",
                                         onlyM = 500,
                                         informative.censoring = TRUE,
                                         observed.covars = TRUE,
                                         #output.directory = "backup-aug-1/simulation",
                                         sim.sample = 1000)

rf.rand.output.n1000 <- follic.output.fun(M = 500,#1000,
                                           fit.initial = "rf",
                                           #onlyM = 500,
                                           informative.censoring = FALSE,
                                           observed.covars = TRUE,
                                           randomized.treatment = TRUE,
                                           observed.treatment = FALSE,
                                           #output.directory = "backup-aug-1/simulation",
                                          sim.sample = 1000)

hal.sl.rand.output.n1000 <- follic.output.fun(M = 500,#1000,
                                              fit.initial = "hal",
                                              hal.sl = TRUE,
                                              #onlyM = 500,
                                              informative.censoring = FALSE,
                                              observed.covars = TRUE,
                                              randomized.treatment = TRUE,
                                              observed.treatment = FALSE,
                                              #output.directory = "backup-aug-1/simulation",
                                              sim.sample = 1000)


hal.rand.output.n1000 <- follic.output.fun(M = 500,#1000,
                                           fit.initial = "hal",
                                           #onlyM = 500,
                                           informative.censoring = FALSE,
                                           observed.covars = TRUE,
                                           randomized.treatment = TRUE,
                                           observed.treatment = FALSE,
                                           #output.directory = "backup-aug-1/simulation",
                                           sim.sample = 1000)

cox.rand.output.n1000 <- follic.output.fun(M = 500,#1000,
                                           fit.initial = "cox",
                                           #onlyM = 500,
                                           informative.censoring = FALSE,
                                           observed.covars = TRUE,
                                           randomized.treatment = TRUE,
                                           observed.treatment = FALSE,
                                           #output.directory = "backup-aug-1/simulation",
                                           sim.sample = 1000)

hal.inf.output.nfollic <- follic.output.fun(M = 500,
                                          fit.initial = "hal",
                                          onlyM = 500,
                                          informative.censoring = TRUE,
                                          observed.covars = TRUE,
                                          #output.directory = "backup-aug-1/simulation",
                                          sim.sample = nrow(follic))

hal.sl.inf.output.nfollic <- follic.output.fun(M = 500,
                                             fit.initial = "hal",
                                             hal.sl = TRUE,
                                             informative.censoring = TRUE,
                                             observed.covars = TRUE,
                                             output.directory = "backup-aug-1/simulation",
                                             sim.sample = nrow(follic))

rf.inf.output.nfollic <- follic.output.fun(M = 500,
                                           fit.initial = "rf",
                                           onlyM = 500,
                                           informative.censoring = TRUE,
                                           observed.covars = TRUE,
                                           #output.directory = "backup-aug-1/simulation",
                                           sim.sample = nrow(follic))

cox.inf.output.n2000 <- follic.output.fun(M = 1000,
                                          onlyM = 500,
                                          fit.initial = "cox",
                                          informative.censoring = TRUE,
                                          observed.covars = TRUE,
                                          output.directory = "backup-aug-1/simulation",
                                          sim.sample = 2000)

cox.inf.output.n2000 <- follic.output.fun(M = 500,#1000,
                                          onlyM = 500,
                                          fit.initial = "cox",
                                          informative.censoring = TRUE,
                                          observed.covars = TRUE,
                                          #output.directory = "backup-aug-1/simulation",
                                          sim.sample = 2000)

survtmle.inf.grid0.n2000 <- follic.output.survtmle(M = 500,
                                                   informative.censoring = TRUE,
                                                   observed.covars = TRUE,
                                                   grid.survtmle = (0:10)*4, 
                                                   sim.sample = 2000)

survtmle.rand.grid0.n2000 <- follic.output.survtmle(M = 500,
                                                    randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                                    observed.covars = TRUE,
                                                    sl.survtmle = FALSE,
                                                    grid.survtmle = (0:10)*4, 
                                                    sim.sample = 2000)

survtmle.inf.grid1.n2000 <- follic.output.survtmle(M = 500,
                                                   informative.censoring = TRUE,
                                                   observed.covars = TRUE,
                                                   grid.survtmle = (0:20)*2, 
                                                   sim.sample = 2000)

survtmle.rand.grid1.n2000 <- follic.output.survtmle(M = 500,
                                                    randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                                    observed.covars = TRUE,
                                                    sl.survtmle = FALSE,
                                                    grid.survtmle = (0:20)*2, 
                                                    sim.sample = 2000)

survtmle.inf.grid2.n2000 <- follic.output.survtmle(M = 500,
                                                   informative.censoring = TRUE,
                                                   observed.covars = TRUE,
                                                   grid.survtmle = (0:40), 
                                                   sim.sample = 2000)

survtmle.rand.grid2.n2000 <- follic.output.survtmle(M = 500,
                                                    randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                                    observed.covars = TRUE,
                                                    sl.survtmle = FALSE,
                                                    grid.survtmle = (0:40), 
                                                    sim.sample = 2000)

survtmle.inf.grid3.n2000 <- follic.output.survtmle(M = 500,
                                                   informative.censoring = TRUE,
                                                   observed.covars = TRUE,
                                                   grid.survtmle = (0:80)/2, 
                                                   sim.sample = 2000)

survtmle.rand.grid3.n2000 <- follic.output.survtmle(M = 500,
                                                    randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                                    observed.covars = TRUE,
                                                    sl.survtmle = FALSE,
                                                    grid.survtmle = (0:80)/2, 
                                                    sim.sample = 2000)

survtmle.inf.grid0.n1000 <- follic.output.survtmle(M = 500,
                                                   informative.censoring = TRUE,
                                                   observed.covars = TRUE,
                                                   grid.survtmle = (0:10)*4, 
                                                   sim.sample = 1000)

survtmle.rand.grid0.n1000 <- follic.output.survtmle(M = 500,
                                                    randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                                    observed.covars = TRUE,
                                                    sl.survtmle = FALSE,
                                                    grid.survtmle = (0:10)*4, 
                                                    sim.sample = 1000)

survtmle.inf.grid1.n1000 <- follic.output.survtmle(M = 500,
                                                   informative.censoring = TRUE,
                                                   observed.covars = TRUE,
                                                   grid.survtmle = (0:20)*2, 
                                                   sim.sample = 1000)

survtmle.rand.grid1.n1000 <- follic.output.survtmle(M = 500,
                                                    randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                                    observed.covars = TRUE,
                                                    sl.survtmle = FALSE,
                                                    grid.survtmle = (0:20)*2, 
                                                    sim.sample = 1000)

survtmle.inf.grid2.n1000 <- follic.output.survtmle(M = 500,
                                                   informative.censoring = TRUE,
                                                   observed.covars = TRUE,
                                                   grid.survtmle = (0:40), 
                                                   sim.sample = 1000)

survtmle.rand.grid2.n1000 <- follic.output.survtmle(M = 500,
                                                    randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                                    observed.covars = TRUE,
                                                    sl.survtmle = FALSE,
                                                    grid.survtmle = (0:40), 
                                                    sim.sample = 1000)

survtmle.inf.grid3.n1000 <- follic.output.survtmle(M = 500,
                                                   informative.censoring = TRUE,
                                                   observed.covars = TRUE,
                                                   grid.survtmle = (0:80)/2, 
                                                   sim.sample = 1000)

survtmle.rand.grid3.n1000 <- follic.output.survtmle(M = 500,
                                                    randomized.treatment = TRUE, observed.treatment = FALSE, informative.censoring = FALSE,
                                                    observed.covars = TRUE,
                                                    sl.survtmle = FALSE,
                                                    grid.survtmle = (0:80)/2, 
                                                    sim.sample = 1000)


cox.inf.output.n1000 <- follic.output.fun(M = 500,#1000,
                                          onlyM = 500,
                                          fit.initial = "cox",
                                          informative.censoring = TRUE,
                                          observed.covars = TRUE,
                                          output.directory = "backup-aug-1/simulation",
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

    print("status = 1:") 
    sim.follic.3[status==1, summary(time)] 
    follic[status==1, summary(time)]

    cbind(sim.data = sim.follic.3[status == 1, summary(time), by = "chemo"], 
          obs.data = follic[status == 1, summary(time), by = "chemo"])

    print("status = 2:")
    sim.follic.3[status==2, summary(time)] 
    follic[status==2, summary(time)]

    cbind(sim.data = sim.follic.3[status == 2, summary(time), by = "chemo"], 
          obs.data = follic[status == 2, summary(time), by = "chemo"])
    
    print("status = 0:")
    sim.follic.3[status==0, summary(time)] 
    follic[status==0, summary(time)]

    cbind(sim.data = sim.follic.3[status == 0, summary(time), by = "chemo"], 
          obs.data = follic[status == 0, summary(time), by = "chemo"])

    sim.follic.3[time <= 10, table(status)]
    follic[time <= 10, table(status)]

    #--- test with independent censoring:

    seed <- sample(10000, 1)
    #seed <- 4660
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

    print("status = 1:") 
    sim.follic.3[status==1, summary(time)] 
    follic[status==1, summary(time)]

    print("status = 2:")
    sim.follic.3[status==2, summary(time)] 
    follic[status==2, summary(time)]

    print("status = 0:")
    sim.follic.3[status==0, summary(time)] 
    follic[status==0, summary(time)]
    
    cbind(sim.data = sim.follic.3[status == 0, summary(time), by = "chemo"], 
          obs.data = follic[status == 0, summary(time), by = "chemo"])

    sim.follic.3[time <= 10, table(status)]
    follic[time <= 10, table(status)]

    #--- test with independent censoring and randomized treatment:

    seed <- sample(10000, 1)
    #seed <- 4660
    sim.follic.3 <- simulate.follic.3(observed.covars = TRUE,
                                      sim.sample = nrow(follic),
                                      counterfactual = NULL,
                                      seed = seed, 
                                      informative.censoring = FALSE,
                                      observed.treatment = FALSE, 
                                      randomized.treatment = TRUE,
                                      keep.times = FALSE)

    sim.follic.3[, table(chemo, status)]
    follic[, table(chemo, status)]

    sim.follic.3[, table(status)]
    follic[, table(status)]

    cbind(sim.data = sim.follic.3[, summary(time), by = "chemo"], 
          obs.data = follic[, summary(time), by = "chemo"])

    print("status = 1:") 
    sim.follic.3[status==1, summary(time)] 
    follic[status==1, summary(time)]

    cbind(sim.data = sim.follic.3[status == 1, summary(time), by = "chemo"], 
          obs.data = follic[status == 1, summary(time), by = "chemo"])

    print("status = 2:")
    sim.follic.3[status==2, summary(time)] 
    follic[status==2, summary(time)]

    cbind(sim.data = sim.follic.3[status == 2, summary(time), by = "chemo"], 
          obs.data = follic[status == 2, summary(time), by = "chemo"])

    print("status = 0:")
    sim.follic.3[status==0, summary(time)] 
    follic[status==0, summary(time)]

    cbind(sim.data = sim.follic.3[status == 0, summary(time), by = "chemo"], 
          obs.data = follic[status == 0, summary(time), by = "chemo"])
    
    sim.follic.3[time <= 10, table(status)]
    follic[time <= 10, table(status)]

    
    #--- test covariate dependent censoring, simulated covariates:
    
    seed <- sample(10000, 1) 
    sim.follic.3 <- simulate.follic.3(observed.covars = FALSE,
                                      sim.sample = nrow(follic),
                                      counterfactual = NULL,
                                      seed = seed,
                                      informative.censoring = TRUE,
                                      keep.times = FALSE)

    sim.follic.3[, table(status)]
    follic[, table(status)]

    sim.follic.3[, summary(time)]  
    follic[, summary(time)]

    print("status = 1:") 
    sim.follic.3[status==1, summary(time)] 
    follic[status==1, summary(time)]

    print("status = 2:")
    sim.follic.3[status==2, summary(time)] 
    follic[status==2, summary(time)]

    print("status = 0:")
    sim.follic.3[status==0, summary(time)] 
    follic[status==0, summary(time)]

    #--- test with independent censoring, simulated covariates:

    seed <- sample(10000, 1) 
    sim.follic.3 <- simulate.follic.3(observed.covars = FALSE,
                                      sim.sample = nrow(follic),
                                      counterfactual = NULL,
                                      seed = seed,
                                      informative.censoring = FALSE,
                                      keep.times = FALSE)

    sim.follic.3[, table(status)]
    follic[, table(status)]

    sim.follic.3[, summary(time)]  
    follic[, summary(time)]

    print("status = 1:") 
    sim.follic.3[status==1, summary(time)] 
    follic[status==1, summary(time)]

    print("status = 2:")
    sim.follic.3[status==2, summary(time)] 
    follic[status==2, summary(time)]

    print("status = 0:")
    sim.follic.3[status==0, summary(time)] 
    follic[status==0, summary(time)]



    
}

######################################################################
### follic.simulation.R ends here
