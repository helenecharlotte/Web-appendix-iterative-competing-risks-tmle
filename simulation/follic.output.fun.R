### follic.output.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Jul 14 2022 (12:51) 
## Version: 
## Last-Updated: Jul 19 2022 (08:30) 
##           By: Helene
##     Update #: 46
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
                              informative.censoring = FALSE,
                              tau = 10,
                              observed.covars = FALSE,
                              observed.treatment = TRUE,
                              randomized.treatment = FALSE,
                              sim.sample = 1000
                              ) {

    print(paste0("./simulation/output/",
                 "outlist-follic-contmle",
                 paste0("-", parameter),
                 paste0("-seed.init", seed.init),
                 paste0("-fit.initial", fit.initial),
                 paste0("-tau", tau),
                 ifelse(informative.censoring, "", "-independentcens"),
                 ifelse(observed.covars, "", "-simulatedcovars"),
                 ifelse(observed.treatment, "",
                 ifelse(randomized.treatment, "-simulatedrandomizedtreatment",
                        "-simulatedtreatment")),
                 ifelse(sim.sample == nrow(follic), "", paste0("-n", sim.sample)),
                 "-M", M, ".rds"))

    print(paste0("modified at: ", file.info(paste0("./simulation/output/",
                                                   "outlist-follic-contmle",
                                                   paste0("-", parameter),
                                                   paste0("-seed.init", seed.init),
                                                   paste0("-fit.initial", fit.initial),
                                                   paste0("-tau", tau),
                                                   ifelse(informative.censoring, "", "-independentcens"),
                                                   ifelse(observed.covars, "", "-simulatedcovars"),
                                                   ifelse(observed.treatment, "",
                                                   ifelse(randomized.treatment, "-simulatedrandomizedtreatment",
                                                          "-simulatedtreatment")),
                                                   ifelse(sim.sample == nrow(follic), "", paste0("-n", sim.sample)),
                                                   "-M", M, ".rds"))$mtime))


    print("#----------------------------------------------------------------------")
    print(paste0("# results for initial fit = * ", fit.initial, " * "))
    print(paste0("# ", ifelse(informative.censoring, "informative censoring", "independent censoring")))
    print(paste0("# ", ifelse(observed.covars, "observed covariates", "simulated covariates")))
    print(paste0("# ", ifelse(observed.treatment, "observed treatment",
                       ifelse(randomized.treatment, "simulated randomized treatment",
                              "simulated treatment"))))
    print("#----------------------------------------------------------------------")

    print(paste0("true.psi = ", true.psi <- readRDS(file=paste0("./simulation/output/",
                                                                "outlist-follic-true-",
                                                                ifelse(parameter == "ate", parameter, paste0("psi", parameter)),
                                                                paste0("-seed.init", seed.init),
                                                                ifelse(observed.covars, "", "-simulatedcovars"),
                                                                paste0("-tau", tau),
                                                                ".rds"))))

    out <- readRDS(file=paste0("./simulation/output/",
                               "outlist-follic-contmle",
                               paste0("-", parameter),
                               paste0("-seed.init", seed.init),
                               paste0("-fit.initial", fit.initial),
                               paste0("-tau", tau),
                               ifelse(informative.censoring, "", "-independentcens"),
                               ifelse(observed.covars, "", "-simulatedcovars"),
                               ifelse(observed.treatment, "",
                               ifelse(randomized.treatment, "-simulatedrandomizedtreatment",
                                      "-simulatedtreatment")),
                               ifelse(sim.sample == nrow(follic), "", paste0("-n", sim.sample)),
                               "-M", M, ".rds"))

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

    return(list(bias = list(tmle = mean(tmle.est-true.psi),
                            init = mean(init.est-true.psi),
                            km = mean(km.est-true.psi)),
                se = list(tmle = mean(tmle.se),
                          km = mean(km.se)),
                sd = list(tmle = sd(tmle.est),
                          init = sd(init.est),
                          km = sd(km.est)),
                mse = list(tmle = mse(tmle.est),
                           km = mse(km.est)),
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

######################################################################
### follic.output.fun.R ends here
