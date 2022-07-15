### follic.simulation.functions.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Jul 14 2022 (11:52) 
## Version: 
## Last-Updated: Jul 15 2022 (13:25) 
##           By: Helene
##     Update #: 24
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#######################################################################################

#--- function to simulate from Weibull (change-points!) 
#

simulate.change.weibull.3 <- function(seed = 100,
                                      follic.sim = follic, #sim.sample = nrow(follic),
                                      change.weibull.parameters = list(
                                          nu.A1.t1 = gamma.status1.t1.1,
                                          nu.A1.t2 = gamma.status1.t2.1,
                                          nu.A1.t3 = gamma.status1.t3.1,
                                          nu.A1.t4 = gamma.status1.t4.1,
                                          nu.A0.t1 = gamma.status1.t1.0,
                                          nu.A0.t2 = gamma.status1.t2.0,
                                          nu.A0.t3 = gamma.status1.t3.0,
                                          nu.A0.t4 = gamma.status1.t4.0,
                                          eta.A1.t1 = lambda.status1.t1.1^gamma.status1.t1.1,
                                          eta.A1.t2 = lambda.status1.t2.1^gamma.status1.t2.1,
                                          eta.A1.t3 = lambda.status1.t3.1^gamma.status1.t3.1,
                                          eta.A1.t4 = lambda.status1.t4.1^gamma.status1.t4.1,
                                          eta.A0.t1 = lambda.status1.t1.0^gamma.status1.t1.0,
                                          eta.A0.t2 = lambda.status1.t2.0^gamma.status1.t2.0,
                                          eta.A0.t3 = lambda.status1.t3.0^gamma.status1.t3.0,
                                          eta.A0.t4 = lambda.status1.t4.0^gamma.status1.t4.0                                            
                                      ),
                                      t0 = log(bhazs[chaz1>0 & chemo==0][["time"]])[kmax.1.t1.0],
                                      t1 = log(bhazs[chaz1>0 & chemo==0][["time"]])[kmax.1.t2.0],
                                      t2 = log(bhazs[chaz1>0 & chemo==0][["time"]])[kmax.1.t3.0],
                                      counterfactual = NULL,
                                      observed.covars = TRUE,
                                      fit.cox = bhaz.cox[["1"]]) {

    set.seed(seed)

    if (length(names(coef(fit.cox)))>0) {
        exp.covar <- follic.sim[, exp(rowSums(sapply(names(coef(fit.cox))[names(coef(fit.cox)) != "chemo" & sapply(names(coef(fit.cox)), function(name) length(grep("period", name)) == 0)], function(covar) {
            coef(fit.cox)[covar]*follic.sim[[covar]]
        })))]
    } else {
        exp.covar <- exp(0)
    }

    if (any(names(change.weibull.parameters)%in%"nu.A1.t4")) {

        Lambda.inv <- function(u, t, nu, eta, nu2=nu, eta2=eta, nu3=nu, eta3=eta, nu4=nu, eta4=eta) {
            return( rowSums(cbind((u <= (eta*exp.covar)*t0^{nu}) *
                                  (( (u + eta*exp.covar*t^{nu}) /
                                     (eta*exp.covar) )^{1/nu} - t),
            (u > (eta*exp.covar)*t0^{nu} & u <= (eta2*exp.covar)*t1^{nu2}) *
            (( (u - (eta2*exp.covar)*t0^{nu2} +
                eta2*exp.covar*t0^{nu2}) /
               (eta2*exp.covar) )^{1/nu2} - t),
            (u > (eta2*exp.covar)*t1^{nu2} & u <= (eta3*exp.covar)*t2^{nu3}) *
            (( (u - (eta3*exp.covar)*t1^{nu3} +
                eta3*exp.covar*t1^{nu3}) /
               (eta3*exp.covar) )^{1/nu3} - t),
            (u > (eta3*exp.covar)*t2^{nu3}) *
            (( (u - (eta4*exp.covar)*t2^{nu4} +
                eta4*exp.covar*t1^{nu4}) /
               (eta4*exp.covar) )^{1/nu4} - t)), na.rm=TRUE) )
        }

        U <- -log(runif(nrow(follic.sim)))

        Tout.A1 <- Lambda.inv(U, 0,
                              nu = change.weibull.parameters[["nu.A1.t1"]],
                              eta = change.weibull.parameters[["eta.A1.t1"]],
                              nu2 = change.weibull.parameters[["nu.A1.t2"]],
                              eta2 = change.weibull.parameters[["eta.A1.t2"]],
                              nu3 = change.weibull.parameters[["nu.A1.t3"]],
                              eta3 = change.weibull.parameters[["eta.A1.t3"]],
                              nu4 = change.weibull.parameters[["nu.A1.t4"]],
                              eta4 = change.weibull.parameters[["eta.A1.t4"]])

        Tout.A0 <- Lambda.inv(U, 0,
                              nu = change.weibull.parameters[["nu.A0.t1"]],
                              eta = change.weibull.parameters[["eta.A0.t1"]],
                              nu2 = change.weibull.parameters[["nu.A0.t2"]],
                              eta2 = change.weibull.parameters[["eta.A0.t2"]],
                              nu3 = change.weibull.parameters[["nu.A0.t3"]],
                              eta3 = change.weibull.parameters[["eta.A0.t3"]],
                              nu4 = change.weibull.parameters[["nu.A0.t4"]],
                              eta4 = change.weibull.parameters[["eta.A0.t4"]])
        
    } else {

        Lambda.inv <- function(u, t, nu, eta, nu2=nu, eta2=eta, nu3=nu, eta3=eta) {
            return( rowSums(cbind((u <= (eta*exp.covar)*t0^{nu}) *
                                  (( (u + eta*exp.covar*t^{nu}) /
                                     (eta*exp.covar) )^{1/nu} - t),
            (u > (eta*exp.covar)*t0^{nu} & u <= (eta2*exp.covar)*t1^{nu2}) *
            (( (u - (eta2*exp.covar)*t0^{nu2} +
                eta2*exp.covar*t0^{nu2}) /
               (eta2*exp.covar) )^{1/nu2} - t),
            (u > (eta2*exp.covar)*t1^{nu2}) *
            (( (u - (eta3*exp.covar)*t1^{nu3} +
                eta3*exp.covar*t1^{nu3}) /
               (eta3*exp.covar) )^{1/nu3} - t)), na.rm=TRUE) )
        }

        U <- -log(runif(nrow(follic.sim)))

        Tout.A1 <- Lambda.inv(U, 0,
                              nu = change.weibull.parameters[["nu.A1.t1"]],
                              eta = change.weibull.parameters[["eta.A1.t1"]],
                              nu2 = change.weibull.parameters[["nu.A1.t2"]],
                              eta2 = change.weibull.parameters[["eta.A1.t2"]],
                              nu3 = change.weibull.parameters[["nu.A1.t3"]],
                              eta3 = change.weibull.parameters[["eta.A1.t3"]])

        Tout.A0 <- Lambda.inv(U, 0,
                              nu = change.weibull.parameters[["nu.A0.t1"]],
                              eta = change.weibull.parameters[["eta.A0.t1"]],
                              nu2 = change.weibull.parameters[["nu.A0.t2"]],
                              eta2 = change.weibull.parameters[["eta.A0.t2"]],
                              nu3 = change.weibull.parameters[["nu.A0.t3"]],
                              eta3 = change.weibull.parameters[["eta.A0.t3"]])
        
    }

    
    if (length(counterfactual) == 0) {
        return(Tout <- follic.sim[["chemo"]]*Tout.A1 + (1-follic.sim[["chemo"]])*Tout.A0)
    } else if (counterfactual == 1) {
        return(Tout.A1) 
    } else {
        return(Tout.A0) 
    }
}

#######################################################################################

#--- function to simulate follic data 
#

simulate.follic.3 <- function(observed.covars = TRUE,
                              sim.sample = nrow(follic),
                              counterfactual = NULL,
                              seed = 100,
                              informative.censoring = TRUE,
                              keep.times = FALSE) {

    set.seed(seed)

    if (observed.covars) {
        if (length(sim.sample)>0) {
            follic.sim <- follic[sample(1:nrow(follic), sim.sample, replace=TRUE)]
        } else {
            follic.sim <- follic
        }
    } else {
        follic.sim <- data.table(stage = rbinom(sim.sample, 1, p.stage),
                                 chemo = rbinom(sim.sample, 1, p.chemo),
                                 age = rnorm(sim.sample, p.age[[1]][["mean"]], p.age[[1]][["sd"]]),
                                 hgb = rnorm(sim.sample, p.hgb[[1]][["mean"]], p.hgb[[1]][["sd"]]))
        follic.sim[, age.squared := age^2]
        follic.sim[, age1 := 1*(age >= 45)]
        follic.sim[, age2 := 1*(age >= 58)]
        follic.sim[, age3 := 1*(age >= 65)]
    }

    if (informative.censoring) {
        T0 <- simulate.change.weibull.3(follic.sim = follic.sim,
                                        seed = seed*1010,
                                        change.weibull.parameters = list(
                                            nu.A1.t1 = gamma.status0.t1.1,
                                            nu.A1.t2 = gamma.status0.t2.1,
                                            nu.A1.t3 = gamma.status0.t3.1,
                                            nu.A0.t1 = gamma.status0.t1.0,
                                            nu.A0.t2 = gamma.status0.t2.0,
                                            nu.A0.t3 = gamma.status0.t3.0,
                                            eta.A1.t1 = lambda.status0.t1.1^gamma.status0.t1.1,
                                            eta.A1.t2 = lambda.status0.t2.1^gamma.status0.t2.1,
                                            eta.A1.t3 = lambda.status0.t3.1^gamma.status0.t3.1,
                                            eta.A0.t1 = lambda.status0.t1.0^gamma.status0.t1.0,
                                            eta.A0.t2 = lambda.status0.t2.0^gamma.status0.t2.0,
                                            eta.A0.t3 = lambda.status0.t3.0^gamma.status0.t3.0                                            
                                        ),
                                        counterfactual = counterfactual,
                                        t0 = (bhazs[chaz0>0 & chemo==0][["time"]])[kmax.0.t1.0],
                                        t1 = (bhazs[chaz0>0 & chemo==0][["time"]])[kmax.0.t2.0], 
                                        fit.cox = bhaz.cox[["0"]])
    } else {
        T0 <- simulate.change.weibull.3(follic.sim = follic.sim,
                                        seed = seed*1010,
                                        change.weibull.parameters = list(
                                            nu.A1.t1 = gamma.status0.independent.t1.1,
                                            nu.A1.t2 = gamma.status0.independent.t2.1,
                                            nu.A1.t3 = gamma.status0.independent.t3.1,
                                            nu.A0.t1 = gamma.status0.independent.t1.0,
                                            nu.A0.t2 = gamma.status0.independent.t2.0,
                                            nu.A0.t3 = gamma.status0.independent.t3.0,
                                            eta.A1.t1 = lambda.status0.independent.t1.1^gamma.status0.independent.t1.1,
                                            eta.A1.t2 = lambda.status0.independent.t2.1^gamma.status0.independent.t2.1,
                                            eta.A1.t3 = lambda.status0.independent.t3.1^gamma.status0.independent.t3.1,
                                            eta.A0.t1 = lambda.status0.independent.t1.0^gamma.status0.independent.t1.0,
                                            eta.A0.t2 = lambda.status0.independent.t2.0^gamma.status0.independent.t2.0,
                                            eta.A0.t3 = lambda.status0.independent.t3.0^gamma.status0.independent.t3.0                                            
                                        ),
                                        counterfactual = counterfactual,
                                        t0 = (bhazs.uninformative.cens[chaz0>0 & chemo==0][["time"]])[kmax.0.t1.0],
                                        t1 = (bhazs.uninformative.cens[chaz0>0 & chemo==0][["time"]])[kmax.0.t2.0], 
                                        fit.cox = bhaz.uninformative.cens[["0"]])
    }

    T1 <- simulate.change.weibull.3(follic.sim = follic.sim,
                                    seed = seed*9010,
                                    change.weibull.parameters = list(
                                        nu.A1.t1 = gamma.status1.t1.1,
                                        nu.A1.t2 = gamma.status1.t2.1,
                                        nu.A1.t3 = gamma.status1.t3.1,
                                        nu.A1.t4 = gamma.status1.t4.1,
                                        nu.A0.t1 = gamma.status1.t1.0,
                                        nu.A0.t2 = gamma.status1.t2.0,
                                        nu.A0.t3 = gamma.status1.t3.0,
                                        nu.A0.t4 = gamma.status1.t4.0,
                                        eta.A1.t1 = lambda.status1.t1.1^gamma.status1.t1.1,
                                        eta.A1.t2 = lambda.status1.t2.1^gamma.status1.t2.1,
                                        eta.A1.t3 = lambda.status1.t3.1^gamma.status1.t3.1,
                                        eta.A1.t4 = lambda.status1.t4.1^gamma.status1.t4.1,
                                        eta.A0.t1 = lambda.status1.t1.0^gamma.status1.t1.0,
                                        eta.A0.t2 = lambda.status1.t2.0^gamma.status1.t2.0,
                                        eta.A0.t3 = lambda.status1.t3.0^gamma.status1.t3.0,
                                        eta.A0.t4 = lambda.status1.t4.0^gamma.status1.t4.0                                            
                                    ),
                                    counterfactual = counterfactual,
                                    t0 = (bhazs[chaz1>0 & chemo==0][["time"]])[kmax.1.t1.0],
                                    t1 = (bhazs[chaz1>0 & chemo==0][["time"]])[kmax.1.t2.0],
                                    t2 = log(bhazs[chaz1>0 & chemo==0][["time"]])[kmax.1.t3.0],
                                    fit.cox = bhaz.cox[["1"]])

    T2 <- simulate.change.weibull.3(follic.sim = follic.sim,
                                    seed = seed*5010,
                                    change.weibull.parameters = list(
                                        nu.A1.t1 = gamma.status2.t1.1,
                                        nu.A1.t2 = gamma.status2.t2.1,
                                        nu.A1.t3 = gamma.status2.t3.1,
                                        nu.A0.t1 = gamma.status2.t1.0,
                                        nu.A0.t2 = gamma.status2.t2.0,
                                        nu.A0.t3 = gamma.status2.t3.0,
                                        eta.A1.t1 = lambda.status2.t1.1^gamma.status2.t1.1,
                                        eta.A1.t2 = lambda.status2.t2.1^gamma.status2.t2.1,
                                        eta.A1.t3 = lambda.status2.t3.1^gamma.status2.t3.1,
                                        eta.A0.t1 = lambda.status2.t1.0^gamma.status2.t1.0,
                                        eta.A0.t2 = lambda.status2.t2.0^gamma.status2.t2.0,
                                        eta.A0.t3 = lambda.status2.t3.0^gamma.status2.t3.0
                                    ),
                                    counterfactual = counterfactual,
                                    t0 = (bhazs[chaz2>0 & chemo==0][["time"]])[kmax.2.t1.0],
                                    t1 = (bhazs[chaz2>0 & chemo==0][["time"]])[kmax.2.t2.0], 
                                    fit.cox = bhaz.cox[["2"]])

    if (length(counterfactual)>0) {
        follic.sim[, time := apply(cbind(T1, T2), 1, min)]
    } else {
        follic.sim[, time := apply(cbind(T1, T2, T0), 1, min)]
    }

    follic.sim[, status := 1*(time == T1) + 2*(time == T2)]

    return(follic.sim)
}
 


######################################################################
### follic.simulation.functions.R ends here
