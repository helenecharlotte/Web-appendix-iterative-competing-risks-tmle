### follic.run.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Jul 14 2022 (12:12) 
## Version: 
## Last-Updated: Jul 18 2022 (09:35) 
##           By: Helene
##     Update #: 13
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

run.follic <- function(M = 1, no_cores = 1, print.m = TRUE, seed.init = 100, no.cores = 1,
                       verbose = FALSE, rf.seed = 1, save.output = TRUE,
                       get.truth = FALSE, parameter = "ate",
                       fit.initial = "cox", tau = 10,
                       observed.covars = TRUE,
                       observed.treatment = TRUE,
                       hal.sl = FALSE, browse.hal = FALSE,
                       sim.sample = nrow(follic),
                       informative.censoring = TRUE,
                       browse = FALSE) {

    if (get.truth) {
        
        set.seed(seed.init)

        if (parameter == "ate" | parameter == "1") {

            sim.follic.1 <- simulate.follic.3(counterfactual = 1, sim.sample = 1e6,
                                              observed.covars = observed.covars)
            true.psi.1 <- mean(sim.follic.1[["time"]] <= tau & sim.follic.1[["status"]] == 1)
            
        }

        if (parameter == "ate" | parameter == "0") {

            sim.follic.0 <- simulate.follic.3(counterfactual = 0, sim.sample = 1e6,
                                              observed.covars = observed.covars)
            true.psi.0 <- mean(sim.follic.0[["time"]] <= tau & sim.follic.0[["status"]] == 1)

        }

        if (parameter == "ate") {
            
            print(paste0("true ate = ",
                         true.ate <- true.psi.1 - true.psi.0))

            saveRDS(true.ate,
                    file=paste0("./simulation/output/",
                                "outlist-follic-true-ate",
                                paste0("-seed.init", seed.init),
                                ifelse(observed.covars, "", "-simulatedcovars"),
                                paste0("-tau", tau),
                                ".rds"))

            return(true.ate)
            
        } else if (parameter == "1") {

            print(paste0("true psi.1 = ", true.psi.1))

            saveRDS(true.psi.1,
                    file=paste0("./simulation/output/",
                                "outlist-follic-true-psi1",
                                paste0("-seed.init", seed.init),
                                ifelse(observed.covars, "", "-simulatedcovars"),
                                paste0("-tau", tau),
                                ".rds"))

            return(true.psi.1)            

        } else {

            print(paste0("true psi.0 = ", true.psi.0))

            saveRDS(true.psi.0,
                    file=paste0("./simulation/output/",
                                "outlist-follic-true-psi0",
                                paste0("-seed.init", seed.init),
                                ifelse(observed.covars, "", "-simulatedcovars"),
                                paste0("-tau", tau),
                                ".rds"))

            return(true.psi.0)            

        }
    }

    if (browse) browser()
    
    registerDoParallel(no.cores)
    
    out <- foreach(m=1:M, .errorhandling="pass"
                   ) %dopar% {

                       if (print.m) print(paste0("m=", m))

                       sim.follic <- simulate.follic.3(seed = m+seed.init,
                                                       sim.sample = sim.sample,
                                                       observed.covars = observed.covars,
                                                       observed.treatment = observed.treatment,
                                                       informative.censoring = informative.censoring)
                     
                       out <- list("est"=contmle(sim.follic, estimation=list("outcome"=list(fit=fit.initial,
                                                                                        model=Surv(time, status==1)~chemo+stage+hgb+age,
                                                                                        lambda.cvs=seq(0.008, 0.02, length=10)),
                                                                         "cens"=list(fit=fit.initial,
                                                                                     model=Surv(time, status==0)~chemo+stage+hgb+age),
                                                                         "cr2"=list(fit=fit.initial,
                                                                                    model=Surv(time, status==2)~chemo+stage+hgb+age)
                                                                         ),
                                                 treat.model=chemo~stage+hgb+age,
                                                 treat.effect=parameter,
                                                 no.small.steps=500,
                                                 sl.models=list(mod1=list(Surv(time, status==1)~chemo+stage+hgb+age, t0 = (1:50)/2000)), 
                                                 output.km=TRUE,
                                                 rf.seed=rf.seed,
                                                 hal.screening=TRUE,
                                                 hal.sl=hal.sl, browse.hal=browse.hal,
                                                 cut.time.grid=7:9, 
                                                 cut.time.covar=5,
                                                 V=3, lambda.cvs=seq(0.1, 0.03, length=10), maxit=1e5, penalize.time=FALSE,
                                                 verbose=verbose,
                                                 iterative=TRUE,
                                                 tau=tau, target=1))

                       names(out) <- paste0("m=", m)
                       print(out)
                       return(out)
                   }

    
    stopImplicitCluster()

    if (save.output) {
        saveRDS(out,
                file=paste0("./simulation/output/",
                            "outlist-follic-contmle",
                            paste0("-", parameter),
                            paste0("-seed.init", seed.init),
                            paste0("-fit.initial", fit.initial),
                            ifelse(hal.sl, "-hal.sl", ""),
                            paste0("-tau", tau),
                            ifelse(informative.censoring, "", "-independentcens"),
                            ifelse(observed.covars, "", "-simulatedcovars"),
                            ifelse(observed.treatment, "", "-simulatedtreatment"),
                            ifelse(sim.sample == nrow(follic), "", paste0("-n", sim.sample)),
                            "-M", M, ".rds"))
    } else {
        return(out)
    }


}

######################################################################
### follic.run.fun.R ends here
