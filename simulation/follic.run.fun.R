### follic.run.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Jul 14 2022 (12:12) 
## Version: 
## Last-Updated: Sep 15 2022 (18:33) 
##           By: Helene
##     Update #: 111
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
                       randomized.treatment = FALSE,
                       hal.sl = FALSE, browse.hal = FALSE,
                       sim.sample = nrow(follic),
                       informative.censoring = TRUE,
                       browse = FALSE,
                       cut.time=18,
                       cut.time.covar=15,
                       cut.time.A=15,
                       cut.covars=15,
                       cut.two.way=10,
                       fit.survtmle = FALSE,
                       grid.survtmle = 0:40,
                       sl.survtmle = FALSE) {

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
                                                       randomized.treatment = randomized.treatment,
                                                       informative.censoring = informative.censoring)

                       if (fit.survtmle) {

                           #-- discretize time variable: 
                           sim.follic[, dtime := as.numeric(cut(time, breaks = grid.survtmle))]

                           glm.time <- paste0(paste0(sapply(1:floor(tau*((length(grid.survtmle)-1)/40)), function(t) {
                               paste0("I(t==", t, ")")
                           }), collapse = "+"), "+",
                           paste0(sapply(1:floor(tau*((length(grid.survtmle)-1)/40)), function(t) {
                               paste0("I(trt*t==", t, ")")
                           }), collapse = "+"), "+ trt + hgb + stage + age")

                           #-- apply survtmle:
                           if (sl.survtmle) {
                               fit.survtmle <- survtmle(ftime = sim.follic[["dtime"]],
                                                        ftype = sim.follic[["status"]],
                                                        trt = sim.follic[["chemo"]],
                                                        adjustVars = data.frame(hgb = sim.follic[["hgb"]],
                                                                                age = sim.follic[["age"]],
                                                                                stage = sim.follic[["stage"]]),
                                                        t0 = floor(tau*((length(grid.survtmle)-1)/40)),
                                                        SL.trt = c("SL.glm","SL.mean","SL.step","SL.ranger"),
                                                        SL.ftime = c("SL.glm","SL.mean","SL.step","SL.ranger"),
                                                        SL.ctime = c("SL.glm","SL.mean","SL.step","SL.ranger"),
                                                        ftypeOfInterest = 1)
                               if (verbose) {
                                   print(fit.survtmle$ftimeMod)
                                   print(fit.survtmle$ctimeMod)
                                   print(fit.survtmle$trtMod)
                               }
                           } else if (informative.censoring) {
                               fit.survtmle <- survtmle(ftime = sim.follic[["dtime"]],
                                                        ftype = sim.follic[["status"]],
                                                        trt = sim.follic[["chemo"]],
                                                        adjustVars = data.frame(hgb = sim.follic[["hgb"]],
                                                                                age = sim.follic[["age"]],
                                                                                stage = sim.follic[["stage"]]),
                                                        t0 = floor(tau*((length(grid.survtmle)-1)/40)),
                                                        glm.trt = "hgb + age + stage",
                                                        glm.ftime = glm.time,
                                                        glm.ctime = glm.time,                         
                                                        ftypeOfInterest = 1)
                           } else {
                               fit.survtmle <- survtmle(ftime = sim.follic[["dtime"]],
                                                        ftype = sim.follic[["status"]],
                                                        trt = sim.follic[["chemo"]],
                                                        adjustVars = data.frame(hgb = sim.follic[["hgb"]],
                                                                                age = sim.follic[["age"]],
                                                                                stage = sim.follic[["stage"]]),
                                                        t0 = floor(tau*((length(grid.survtmle)-1)/40)),
                                                        glm.trt = "hgb + age + stage",
                                                        glm.ftime = glm.time,
                                                        glm.ctime = glm.time,                         
                                                        ftypeOfInterest = 1)
                           }

                           if (parameter == "1") {
                               est <- fit.survtmle$est[[2]]
                               se <- mean(unlist(fit.survtmle$ic[,"D.j1.z1"])^2)/sqrt(nrow(sim.follic))
                           } else if (parameter == "0") {
                               est <- fit.survtmle$est[[1]]
                               se <- mean(unlist(fit.survtmle$ic[,"D.j1.z0"])^2)/sqrt(nrow(sim.follic))
                           } else {
                               est <- fit.survtmle$est[[2]] - fit.survtmle$est[[1]]
                               se <- mean((unlist(fit.survtmle$ic[,"D.j1.z1"])-unlist(fit.survtmle$ic[,"D.j1.z0"]))^2)/sqrt(nrow(sim.follic))
                           }

                           out <- list("est" = c(survtmle.est = est, survtmle.se = se))
                           names(out) <- paste0("m=", m)
                           
                           print(out)
                           return(out)
                           
                       } else {                       
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
                                                     #cut.time.grid=7:9,
                                                     cut.time=cut.time,
                                                     cut.time.covar=cut.time.covar,
                                                     cut.time.A=cut.time.A,
                                                     cut.covars=cut.covars,
                                                     cut.two.way=cut.two.way,
                                                     V=3, lambda.cvs=seq(0.1, 0.03, length=10), maxit=1e5, penalize.time=FALSE,
                                                     verbose=verbose,
                                                     iterative=TRUE,
                                                     tau=tau, target=1))
                       }

                       names(out) <- paste0("m=", m)
                       print(out)
                       return(out)
                   }

    
    stopImplicitCluster()

    if (save.output) {
        saveRDS(out,
                file=paste0("./simulation/output/",
                            "outlist-follic-contmle",
                            ifelse(fit.survtmle, "-survtmle", ""),
                            ifelse(fit.survtmle, paste0("-gridlength", length(grid.survtmle)), ""),
                            ifelse(fit.survtmle & sl.survtmle, "-sl", ""),                            
                            paste0("-", parameter),
                            paste0("-seed.init", seed.init),
                            paste0("-fit.initial", fit.initial),
                            ifelse(hal.sl, "-hal.sl", ""),
                            ifelse(cut.time == 18, "", paste0("-cut.time", cut.time)),
                            ifelse(cut.time.covar == 15, "", paste0("-cut.time.covar", cut.time.covar)),
                            ifelse(cut.time.A == 15, "", paste0("-cut.time.A", cut.time.A)),
                            ifelse(cut.covars == 15, "", paste0("-cut.covars", cut.covars)),
                            ifelse(cut.two.way == 10, "", paste0("-cut.two.way", cut.two.way)),
                            paste0("-tau", tau),
                            ifelse(informative.censoring, "", "-independentcens"),
                            ifelse(observed.covars, "", "-simulatedcovars"),
                            ifelse(observed.treatment, "",
                            ifelse(randomized.treatment, "-simulatedrandomizedtreatment",
                                   "-simulatedtreatment")),
                            ifelse(sim.sample == nrow(follic), "", paste0("-n", sim.sample)),
                            "-M", M, ".rds"))
    } else {
        return(out)
    }


}

######################################################################
### follic.run.fun.R ends here
