### follic.output.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Jul 14 2022 (12:51) 
## Version: 
## Last-Updated: Sep 16 2022 (08:04) 
##           By: Helene
##     Update #: 165
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#######################################################################################

mse <- function(x) mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE)

#######################################################################################

follic.output.fun <- function(M = 500,
                              parameter = "ate",
                              seed.init = 100,
                              fit.initial = "rf",
                              hal.sl = FALSE,
                              cut.time=18,
                              cut.time.covar=15,
                              cut.time.A=15,
                              cut.covars=15,
                              cut.two.way=10,
                              informative.censoring = FALSE,
                              tau = 10,
                              observed.covars = FALSE,
                              observed.treatment = TRUE,
                              randomized.treatment = FALSE,
                              sim.sample = 1000,
                              output.directory = "simulation",
                              browse = FALSE,
                              onlyM = M,
                              se.bound = 0.06
                              ) {

    print(paste0("./", output.directory, "/output/",
                 "outlist-follic-contmle",
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

    print(paste0("modified at: ", file.info(paste0("./", output.directory, "/output/",
                                                   "outlist-follic-contmle",
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
                                                   "-M", M, ".rds"))$mtime))


    print("#----------------------------------------------------------------------")
    print(paste0("# results for initial fit = * ", ifelse(hal.sl, "hal-SL", fit.initial), " * "))
    print(paste0("# ", ifelse(informative.censoring, "informative censoring", "independent censoring")))
    print(paste0("# ", ifelse(observed.covars, "observed covariates", "simulated covariates")))
    print(paste0("# ", ifelse(observed.treatment, "observed treatment",
                       ifelse(randomized.treatment, "simulated randomized treatment",
                              "simulated treatment"))))
    print("#----------------------------------------------------------------------")

    out <- readRDS(file=paste0("./", output.directory, "/output/",
                               "outlist-follic-contmle",
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

    if (onlyM<M) {
        out <- out[1:onlyM]
        print(paste0("(only look at first ", length(out), " runs)"))
        print("#----------------------------------------------------------------------")
    }
    
    print(paste0("true.psi = ", true.psi <- readRDS(file=paste0("./", output.directory, "/output/",
                                                                "outlist-follic-true-",
                                                                ifelse(parameter == "ate", parameter, paste0("psi", parameter)),
                                                                paste0("-seed.init", seed.init),
                                                                ifelse(observed.covars, "", "-simulatedcovars"),
                                                                paste0("-tau", tau),
                                                                ".rds"))))

    print(paste0("bias tmle = ", mean(tmle.est <- unlist(lapply(out, function(out1) {
        if (is.list(out1[[1]])) {
            return(out1[[1]]$tmle$F1["tmle.est",][[1]])
        } else return(NA)
    })), na.rm = TRUE)-true.psi))
    print(paste0("bias init = ", mean(init.est <- unlist(lapply(out, function(out1) {
        if (is.list(out1[[1]])) {
            return(out1[[1]]$init$F1["init.est",][[1]])
        } else return(NA)
    })), na.rm = TRUE)-true.psi))
    print(paste0("bias km   = ", mean(km.est <- unlist(lapply(out, function(out1) {
        if (is.list(out1[[1]])) {
            return(out1[[1]]$km$F1["km.est",][[1]])
        } else return(NA)
    })), na.rm = TRUE)-true.psi))

    print(paste0("fraction failed to produce output = ", mean(is.na(tmle.est))))
    if (mean(is.na(tmle.est))>0) print(paste0("(for m = ", paste0((1:length(out))[is.na(tmle.est)], collapse = ","), ")"))
    if (browse) browser()

    par(mfrow = c(1,4))
    hist(tmle.est)
    abline(v = true.psi, col = "red")
    abline(v = mean(tmle.est), col = "blue")
    hist(init.est)
    abline(v = true.psi, col = "red")
    abline(v = mean(init.est), col = "blue")
    hist(km.est)
    abline(v = true.psi, col = "red")
    abline(v = mean(km.est), col = "blue")

    #-- se 
    print(paste0("mean tmle.se = ", mean(tmle.se <- unlist(lapply(out, function(out1) {
        if (is.list(out1[[1]])) {
            return(out1[[1]]$tmle$F1["tmle.se",][[1]])
        } else return(NA)
    })), na.rm = TRUE)))
    hist(tmle.se)
    abline(v = sd(tmle.est, na.rm = TRUE), col = "red")
    abline(v = mean(tmle.se, na.rm = TRUE), col = "blue")
    print(paste0("tmle mc sd = ", sd(tmle.est, na.rm = TRUE)))
    print(paste0("init mc sd = ", sd(init.est, na.rm = TRUE)))
    print(paste0("km mc sd   = ", sd(km.est, na.rm = TRUE)))

    km.se <- unlist(lapply(out, function(out1) {
        if (is.list(out1[[1]])) {
            return(out1[[1]]$km$F1["km.se",][[1]])
        } else return(NA)
    }))

    print(paste0("mse (tmle) = ", mse(tmle.est)))
    print(paste0("mse (tmle/km) = ", mse(tmle.est)/mse(km.est)))

    #-- km coverage
    print(paste0("km coverage = ", km.coverage <- mean(km.est - 1.96*km.se <= true.psi &
                                                       true.psi <= km.est + 1.96*km.se, na.rm = TRUE)))

    #-- oracle coverage
    print(paste0("oracle coverage = ", mean(tmle.est - 1.96*tmle.se <=  mean(tmle.est, na.rm = TRUE) &
         mean(tmle.est, na.rm = TRUE) <= tmle.est + 1.96*tmle.se, na.rm = TRUE)))

    #-- coverage
    print(paste0("coverage = ", coverage <- mean(tmle.est - 1.96*tmle.se <= true.psi &
                                                 true.psi <= tmle.est + 1.96*tmle.se, na.rm = TRUE)))

    if (any(tmle.se[!is.na(tmle.se)]>se.bound) & sim.sample>nrow(follic)) {
        print("---------------------")
        print(paste0("no of tmle.se>", se.bound, " = ", sum(tmle.se>se.bound)))
        print(paste0("(for m = ", paste0((1:length(out))[tmle.se>se.bound], collapse = ","), ")"))
        print(paste0("bias tmle.est = ", mean(tmle.est[tmle.se<se.bound]-true.psi, na.rm = TRUE)))
        print(paste0("mean tmle.se = ", mean(tmle.se[tmle.se<se.bound], na.rm = TRUE)))
        print(paste0("init mc sd = ", sd(init.est[tmle.se<se.bound], na.rm = TRUE)))
        print(paste0("tmle mc sd = ", sd(tmle.est[tmle.se<se.bound], na.rm = TRUE)))
        print(paste0("mse (tmle/km) = ", mse(tmle.est[tmle.se<se.bound])/mse(km.est[tmle.se<se.bound])))
        print(paste0("coverage = ", mean((tmle.est - 1.96*tmle.se <= true.psi &
                                          true.psi <= tmle.est + 1.96*tmle.se)[tmle.se<se.bound], na.rm = TRUE)))
        print("---------------------")
    }

    return(list(bias = list(tmle = mean(tmle.est-true.psi, na.rm = TRUE),
                            init = mean(init.est-true.psi, na.rm = TRUE),
                            km = mean(km.est-true.psi, na.rm = TRUE)),
                se = list(tmle = mean(tmle.se, na.rm = TRUE),
                          km = mean(km.se, na.rm = TRUE)),
                sd = list(tmle = sd(tmle.est, na.rm = TRUE),
                          init = sd(init.est, na.rm = TRUE),
                          km = sd(km.est, na.rm = TRUE)),
                mse = list(tmle = mse(tmle.est)*100,
                           km = mse(km.est)*100),
                cov = list(tmle = coverage, km = km.coverage)))

}

follic.compare.results <- function(hal.output, rf.output, cox.output) {
    if (rf.output$mse$km / hal.output$mse$km != 1) warning("check if hal and rf were run on same data")
    return(cbind(hal = c(bias = hal.output$bias$tmle,
                         cov = hal.output$cov$tmle,
                         sd = hal.output$sd$tmle,
                         se = hal.output$se$tmle, 
                         mse = hal.output$mse$tmle,
                         mse.km = hal.output$mse$tmle/hal.output$mse$km),
                 rf = c(bias = rf.output$bias$tmle,
                        cov = rf.output$cov$tmle,
                        sd = rf.output$sd$tmle,
                        se = rf.output$se$tmle, 
                        mse = rf.output$mse$tmle,
                        mse.km = rf.output$mse$tmle/rf.output$mse$km),
                 km = c(bias = hal.output$bias$km,
                        cov = hal.output$cov$km,
                        sd = hal.output$sd$km,
                        se = hal.output$se$km, 
                        mse = hal.output$mse$km,
                        mse.km = rf.output$mse$km / hal.output$mse$km),
                 cox = c(bias = cox.output$bias$tmle,
                         cov = cox.output$cov$tmle,
                         sd = cox.output$sd$tmle,
                         se = cox.output$se$tmle, 
                         mse = cox.output$mse$tmle,
                         mse.km = cox.output$mse$tmle/cox.output$mse$km)))
}

follic.output.survtmle <- function(M = 500,
                                   parameter = "ate",
                                   seed.init = 100,
                                   fit.initial = "cox",
                                   informative.censoring = FALSE,
                                   tau = 10,
                                   browse = FALSE,
                                   observed.covars = FALSE,
                                   observed.treatment = TRUE,
                                   randomized.treatment = FALSE,
                                   fit.survtmle = TRUE,
                                   grid.survtmle = (0:40),
                                   survtmle.glm.misspecify = FALSE,
                                   sl.survtmle = FALSE,
                                   sim.sample = 1000,
                                   output.directory = "simulation"
                                   ) {

    print(paste0("./", output.directory, "/output/",
                 "outlist-follic-contmle-survtmle",
                 ifelse(fit.survtmle & survtmle.glm.misspecify & !sl.survtmle, "-misspecifyGLM", ""),
                 ifelse(fit.survtmle, paste0("-gridlength", length(grid.survtmle)), ""),
                 ifelse(fit.survtmle & sl.survtmle, "-sl", ""),                            
                 paste0("-", parameter),
                 paste0("-seed.init", seed.init),
                 paste0("-fit.initial", fit.initial),
                 #ifelse(hal.sl, "-hal.sl", ""),
                 paste0("-tau", tau),
                 ifelse(informative.censoring, "", "-independentcens"),
                 ifelse(observed.covars, "", "-simulatedcovars"),
                 ifelse(observed.treatment, "",
                 ifelse(randomized.treatment, "-simulatedrandomizedtreatment",
                        "-simulatedtreatment")),
                 ifelse(sim.sample == nrow(follic), "", paste0("-n", sim.sample)),
                 "-M", M, ".rds"))

    print(paste0("modified at: ", file.info(paste0("./", output.directory, "/output/",
                                                   "outlist-follic-contmle-survtmle",
                                                   ifelse(fit.survtmle & survtmle.glm.misspecify & !sl.survtmle, "-misspecifyGLM", ""),
                                                   ifelse(fit.survtmle, paste0("-gridlength", length(grid.survtmle)), ""),
                                                   ifelse(fit.survtmle & sl.survtmle, "-sl", ""),                            
                                                   paste0("-", parameter),
                                                   paste0("-seed.init", seed.init),
                                                   paste0("-fit.initial", fit.initial),
                                                   #ifelse(hal.sl, "-hal.sl", ""),
                                                   paste0("-tau", tau),
                                                   ifelse(informative.censoring, "", "-independentcens"),
                                                   ifelse(observed.covars, "", "-simulatedcovars"),
                                                   ifelse(observed.treatment, "",
                                                   ifelse(randomized.treatment, "-simulatedrandomizedtreatment",
                                                          "-simulatedtreatment")),
                                                   ifelse(sim.sample == nrow(follic), "", paste0("-n", sim.sample)),
                                                   "-M", M, ".rds"))$mtime))

    out <- readRDS(file=paste0("./", output.directory, "/output/",
                               "outlist-follic-contmle-survtmle",
                               ifelse(fit.survtmle & survtmle.glm.misspecify & !sl.survtmle, "-misspecifyGLM", ""),
                               ifelse(fit.survtmle, paste0("-gridlength", length(grid.survtmle)), ""),
                               ifelse(fit.survtmle & sl.survtmle, "-sl", ""),                            
                               paste0("-", parameter),
                               paste0("-seed.init", seed.init),
                               paste0("-fit.initial", fit.initial),
                               #ifelse(hal.sl, "-hal.sl", ""),
                               paste0("-tau", tau),
                               ifelse(informative.censoring, "", "-independentcens"),
                               ifelse(observed.covars, "", "-simulatedcovars"),
                               ifelse(observed.treatment, "",
                               ifelse(randomized.treatment, "-simulatedrandomizedtreatment",
                                      "-simulatedtreatment")),
                               ifelse(sim.sample == nrow(follic), "", paste0("-n", sim.sample)),
                               "-M", M, ".rds"))

    if (browse) browser()

    print("#----------------------------------------------------------------------")
    print(paste0("# results for survtmle"))
    print(paste0("# ", ifelse(informative.censoring, "informative censoring", "independent censoring")))
    print(paste0("# ", ifelse(observed.covars, "observed covariates", "simulated covariates")))
    print(paste0("# ", ifelse(observed.treatment, "observed treatment",
                       ifelse(randomized.treatment, "simulated randomized treatment",
                              "simulated treatment"))))
    print("#----------------------------------------------------------------------")

    print(paste0("true.psi = ", true.psi <- readRDS(file=paste0("./", output.directory, "/output/",
                                                                "outlist-follic-true-",
                                                                ifelse(parameter == "ate", parameter, paste0("psi", parameter)),
                                                                paste0("-seed.init", seed.init),
                                                                ifelse(observed.covars, "", "-simulatedcovars"),
                                                                paste0("-tau", tau),
                                                                ".rds"))))

    survtmle.nas <- is.na(unlist(sapply(out, function(out1) out1[[1]]["survtmle.est"])))
    print(paste0("number of na's = ", sum(survtmle.nas)))
    print(paste0("(for m = ", paste0((1:length(out))[survtmle.nas], collapse = ","), ")"))
    if (sum(survtmle.nas)>0) out <- out[-survtmle.nas]
    
    print(paste0("bias survtmle = ", mean(survtmle.est <- as.numeric(unlist(sapply(out, function(out1) out1[[1]]["survtmle.est"]))), na.rm = TRUE)-true.psi))
    print(paste0("mean survtmle.se = ", mean(survtmle.se <- as.numeric(unlist(sapply(out, function(out) out[[1]]["survtmle.se"]))), na.rm = TRUE)))
    print(paste0("sd survtmle = ", sd(survtmle.est, na.rm = TRUE)))
    print(paste0("mse survtmle = ", mse(survtmle.est)))
    print(paste0("coverage = ", coverage <- mean(survtmle.est - 1.96*survtmle.se <= true.psi &
                                                 true.psi <= survtmle.est + 1.96*survtmle.se, na.rm = TRUE)))

    return(list(bias = list(survtmle = mean(survtmle.est-true.psi, na.rm = TRUE)),
                cov = list(survtmle = coverage),
                se = list(survtmle = mean(survtmle.se, na.rm = TRUE)),
                sd = list(survtmle = sd(survtmle.est, na.rm = TRUE)),
                mse = list(survtmle = mse(survtmle.est)*100),
                mse.km = NA))
    
}

######################################################################
### follic.output.fun.R ends here
