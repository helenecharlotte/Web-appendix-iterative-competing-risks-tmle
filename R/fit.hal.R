##' bla
##'
##' blub
##' @title atitle
##' @param covars covariates to include in HAL model. 
##' @param dt dataset. 
##' @param treatment name of treatment variable. 
##' @param V number of folds in cross-validation. 
##' @param cut.one.way cut-offs for main effect indicators. 
##' @param method.risk option to pick different cross-validation schemes for cox models; basically it picks the
##' risk set for the partial likelihood. Should be chosen as 'test'.
##' @param cv.glmnet if TRUE, default of glmnet is used to find penalization (do not want this). 
##' @param seed random seed :). 
##' @param grouped option to set cross-validation method in glmnet. Keep grouped=TRUE. 
##' @param use.min combind with cv.glmnet, uses the minimal optimal value of penalization.
##' @param penalize.treatment option to penalize/not penalize treatment indicators. 
##' @param penalize.time option to penalize/not penalize time indicators. 
##' @param treatment.prediction name of treatment variable.
##' @param intervention interventional values for treatment variable. 
##' @param delta.var name of event type variable. 
##' @param delta.value type of event of interest here-
##' @param time.var name of time varaible. 
##' @param cut.time cut-off for time variable. 
##' @param cut.time.treatment cut-off for time/treatment interaction.
##' @param browse browser() for helene :). 
##' @param predict time-horizon at which to predict (if want prediction). 
##' @param mat dataset with information for HAL.
##' @param verbose produce comments for code throughout. 
##' @param lambda.cvs grid over which to choose penalization. 
##' @param two.way two-way interactions of variables to include. 
##' @param cut.two.way cut-offs for two-way interactions.
##' @return a
##' @seealso b
##' @examples c
##' @export 
##' @author Helene C. W. Rytgaard <hely@@biostat.ku.dk>
fit.hal <- function(covars, dt, treatment=NULL, V=5, cut.one.way=8, method.risk="test", cv.glmnet=FALSE,
                    seed=13349, grouped=TRUE, use.min=TRUE, penalize.treatment=FALSE,
                    penalize.time=TRUE, treatment.prediction=treatment,
                    intervention=c(0,1), 
                    delta.var="delta", delta.value=1,
                    time.var=NULL, cut.time=5, cut.time.treatment=NULL,
                    cut.time.covar=NULL, cut.one.way.time=cut.one.way,
                    browse=FALSE, browse1=FALSE, print.extra=FALSE,
                    predict=FALSE, mat=NULL, verbose=FALSE,
                    lambda.cvs=c(sapply(1:5, function(jjj) (9:1)/(10^jjj))),
                    two.way=cbind(var1="", var2=""), cut.two.way=5,
                    return.cve=FALSE) {

    set.seed(seed)

    if (browse1) browser()
    
    X.hal <- basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=dt,
                       time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                       cut.time.covar=cut.time.covar, cut.one.way.time=cut.one.way.time,
                       delta.var=delta.var, delta.value=delta.value, pseudo.dt=mat,
                       two.way=two.way, cut.two.way=cut.two.way)
    
    if (is.list(X.hal) & length(X.hal)==2) {
        pseudo.dt <- X.hal$pseudo.dt
        X.hal <- X.hal$X
        pseudo.dt[, D:=sum(time.obs==grid.time & get(delta.var)==delta.value), by="x"]
        pseudo.dt[, RT:=sum(risk.time), by="x"]
        if (print.extra) print(paste0("nrow pseudo.dt = ", nrow(pseudo.dt)))
    } else {
        Y.hal <- dt[, Surv(time, get(delta.var)==delta.value)]
    }

    penalty.factor <- rep(1, ncol(X.hal))
   
    if (!penalize.treatment & length(grep(treatment, colnames(X.hal)))>0) {
        penalty.factor[colnames(X.hal) %in% grep(treatment, colnames(X.hal), value=TRUE)] <- 0
        penalty.factor <- 1-(cumprod(1-penalty.factor))
    }

    if (!penalize.time & length(grep("grid.period", colnames(X.hal)))>0) {
        not.penalize <- grep("grid.period", colnames(X.hal), value=TRUE)
        if (length(grep(":", not.penalize))>0) not.penalize <- not.penalize[-grep(":", not.penalize)]
        penalty.factor[colnames(X.hal) %in% not.penalize] <- 0
    }    

    if (length(time.var)>0) {
        tmp.dt <- unique(pseudo.dt[, c("x", "D", "RT")])
        Y2.hal <- tmp.dt[RT>0, D]
        offset2 <- tmp.dt[RT>0, log(RT)]
        X2.hal <- unique.matrix(X.hal)[tmp.dt$RT>0,]
        if (cv.glmnet) {
            hal.fit <- cv.glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                                 offset=offset2,
                                 family="poisson",
                                 penalty.factor=penalty.factor,
                                 maxit=10000)
            if (use.min) {
                lambda.cv <- hal.fit$lambda.min
                hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                                  offset=offset2,
                                  family="poisson",
                                  lambda=lambda.cv,
                                  penalty.factor=penalty.factor,
                                  maxit=10000)
                cve.hal <- cv.fun(loss.fun=lebesgue.loss.fun, dt=pseudo.dt, X=X.hal, Y=NULL, time.var=time.var,
                                  penalty.factor=penalty.factor,
                                  offset=TRUE, V=V,lambda.cvs=lambda.cv, delta.var=delta.var, delta.value=delta.value)
            } else {
                cve.hal <- cv.fun(loss.fun=lebesgue.loss.fun, dt=pseudo.dt, X=X.hal, Y=NULL, time.var=time.var,
                                  penalty.factor=penalty.factor,
                                  offset=TRUE, V=V,lambda.cvs=hal.fit$lambda.1se, delta.var=delta.var, delta.value=delta.value)
            }
        } else {
            cve.hal <- cv.fun(loss.fun=lebesgue.loss.fun, dt=pseudo.dt, X=X.hal, Y=NULL, time.var=time.var,
                              penalty.factor=penalty.factor, lambda.cvs=lambda.cvs,
                              offset=TRUE, V=V, delta.var=delta.var, delta.value=delta.value)
            hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                              offset=offset2,,
                              lambda=cve.hal$min$lambda.cv, 
                              family="poisson",
                              penalty.factor=penalty.factor,
                              maxit=10000)
            repeat {
                if (sum(coef(hal.fit))==0) { # if just empty model is fit
                    cve.hal <- cv.fun(loss.fun=lebesgue.loss.fun, dt=pseudo.dt, X=X.hal, Y=NULL, time.var=time.var,
                                      penalty.factor=penalty.factor, lambda.cvs=lambda.cvs[lambda.cvs>cve.hal$min$lambda.cv],
                                      offset=TRUE, V=V, delta.var=delta.var, delta.value=delta.value)
                    hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                                      offset=offset2,,
                                      lambda=cve.hal$min$lambda.cv, 
                                      family="poisson",
                                      penalty.factor=penalty.factor,
                                      maxit=10000)
                } else {
                    break
                }
            }
            lambda.cv <- cve.hal$min$lambda.cv
            if (browse) {
                lambda.cvs <- c(sapply(1:5, function(jjj) (9:1)/(10^jjj)))
                hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                                  offset=offset2,,
                                  lambda=lambda.cvs,
                                  #cve.hal$min$lambda.cv, 
                                  family="poisson",
                                  penalty.factor=penalty.factor,
                                  maxit=10000)
                lapply(lambda.cvs, function(lambda.cv) covars[sapply(covars, function(covar) length(grep(covar, coef(hal.fit, s=lambda.cv)@Dimnames[[1]][coef(hal.fit, s=lambda.cv)@i+1]))>0)])
            }
        }
    } else if (cv.glmnet) {
        hal.fit <- cv.glmnet(x=as.matrix(X.hal), y=Y.hal,
                             family="cox",
                             grouped=grouped,
                             penalty.factor=penalty.factor,
                             maxit=10000)
        if (use.min) {
            lambda.cv <- hal.fit$lambda.min
            hal.fit <- glmnet(x=as.matrix(X.hal), y=Y.hal,
                              family="cox",
                              lambda=lambda.cv,
                              grouped=grouped,
                              penalty.factor=penalty.factor,
                              maxit=10000)
            cve.hal <- cv.fun(loss.fun=cox.loss.fun, dt=dt, X=X.hal, Y=Y.hal, method.risk=method.risk, V=V,
                              penalty.factor=penalty.factor,
                              lambda.cvs=lambda.cv, delta.var=delta.var, delta.value=delta.value)
        } else {
            cve.hal <- cv.fun(loss.fun=cox.loss.fun, dt=dt, X=X.hal, Y=Y.hal, method.risk=method.risk, V=V,
                              penalty.factor=penalty.factor,
                              lambda.cvs=hal.fit$lambda.1se, delta.var=delta.var, delta.value=delta.value)
        }
    } else {
        cve.hal <- cv.fun(loss.fun=cox.loss.fun, dt=dt, X=X.hal, Y=Y.hal, lambda.cvs=lambda.cvs,
                          method.risk=method.risk, V=V, delta.var=delta.var, delta.value=delta.value)
        hal.fit <- glmnet(x=as.matrix(X.hal), y=Y.hal,
                          lambda=cve.hal$min$lambda.cv, 
                          family="cox",
                          penalty.factor=penalty.factor,
                          maxit=10000)
        lambda.cv <- cve.hal$min$lambda.cv
        if (browse) {
            browser()
            lambda.cvs <- c(sapply(1:5, function(jjj) (9:1)/(10^jjj)))
            hal.fit <- glmnet(x=as.matrix(X.hal), y=Y.hal,
                              lambda=cve.hal$min$lambda.cvs, 
                              family="cox",
                              penalty.factor=penalty.factor,
                              maxit=10000)
            lapply(lambda.cvs, function(lambda.cv) covars[sapply(covars, function(covar) length(grep(covar, coef(hal.fit, s=lambda.cv)@Dimnames[[1]][coef(hal.fit, s=lambda.cv)@i+1]))>0)])
        }
    }

    if (verbose) {
        print(coef(hal.fit, s=lambda.cv))
    }

    if (predict) {

        if (length(time.var)>0) {

            if (length(mat)>0) {
                X.hal.a <- lapply(intervention, function(a) {
                    return(basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=copy(dt)[,(treatment):=a],
                                     time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                                     cut.time.covar=cut.time.covar, cut.one.way.time=cut.one.way.time,
                                     predict=predict,
                                     pseudo.dt=copy(mat)[,(treatment):=a],
                                     delta.var=delta.var, delta.value=delta.value,
                                     two.way=two.way, cut.two.way=cut.two.way))
                })
                ## X.hal1 <- basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=copy(dt)[,(treatment):=1],
                ##                     time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                ##                     predict=predict,
                ##                     pseudo.dt=copy(mat)[,(treatment):=1],
                ##                     delta.var=delta.var, delta.value=delta.value,
                ##                     two.way=two.way, cut.two.way=cut.two.way)
                ## X.hal0 <- basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=copy(dt)[,(treatment):=0],
                ##                     time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                ##                     predict=predict,
                ##                     pseudo.dt=copy(mat)[,(treatment):=0], 
                ##                     delta.var=delta.var, delta.value=delta.value,
                ##                     two.way=two.way, cut.two.way=cut.two.way)$X
            } else {
                X.hal.a <- lapply(intervention, function(a) {
                    return(basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=copy(dt)[,(treatment):=a],
                                     time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                                     cut.time.covar=cut.time.covar, cut.one.way.time=cut.one.way.time,
                                     predict=predict,
                                     delta.var=delta.var, delta.value=delta.value,
                                     two.way=two.way, cut.two.way=cut.two.way))
                })
                ## X.hal1 <- basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=copy(dt)[,(treatment):=1],
                ##                     time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                ##                     predict=predict,
                ##                     delta.var=delta.var, delta.value=delta.value,
                ##                     two.way=two.way, cut.two.way=cut.two.way)
                ## X.hal0 <- basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=copy(dt)[,(treatment):=0],
                ##                     time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                ##                     predict=predict,
                ##                     delta.var=delta.var, delta.value=delta.value,
                ##                     two.way=two.way, cut.two.way=cut.two.way)$X
            }

            pseudo.dt <- X.hal.a[[1]]$pseudo.dt
            X.hal.a <- lapply(1:length(intervention), function(aa) X.hal.a[[aa]]$X)            
            ## pseudo.dt <- X.hal1$pseudo.dt
            ## X.hal1 <- X.hal1$X

            for (aa in 1:length(intervention)) {
                pseudo.dt[, (paste0("fit.lambda", intervention[aa])):=exp(predict(hal.fit, X.hal.a[[aa]],
                                                                                  newoffset=0, s=lambda.cv))]
            }
            ## pseudo.dt[, fit.lambda1:=exp(predict(hal.fit, X.hal1,
            ##                                      newoffset=0, s=lambda.cv))]
            ## pseudo.dt[, fit.lambda0:=exp(predict(hal.fit, X.hal0,
            ##                                      newoffset=0, s=lambda.cv))]
           
            if (length(mat)>0) {

                mat <- do.call("rbind", lapply(intervention, function(a) {
                    merge(mat[get(treatment.prediction)==a],
                          pseudo.dt[, c("id", covars, "grid.period", paste0("fit.lambda", a)), with=FALSE],
                          by=c("id", covars, "grid.period"))[, fit.lambda:=get(paste0("fit.lambda", a))][, (paste0("fit.lambda", a)):=NULL]
                }))
                ## mat <- rbind(merge(mat[get(treatment.prediction)==1],
                ##                    pseudo.dt[, c("id", covars, "grid.period", "fit.lambda1"), with=FALSE],
                ##                    by=c("id", covars, "grid.period")),
                ##              merge(mat[get(treatment.prediction)==0],
                ##                    pseudo.dt[, c("id", covars, "grid.period", "fit.lambda0"), with=FALSE],
                ##                    by=c("id", covars, "grid.period")), fill=TRUE)
                ## mat[get(treatment.prediction)==1, fit.lambda:=fit.lambda1]
                ## mat[get(treatment.prediction)==0, fit.lambda:=fit.lambda0]

                ## mat[, fit.lambda1:=NULL]
                ## mat[, fit.lambda0:=NULL]
                
                mat[, risk.time:=c(0,diff(time)), by=c("id", treatment.prediction)]
                mat[, fit.dLambda:=fit.lambda*risk.time]
                mat[, fit.Lambda:=cumsum(fit.dLambda), by=c("id", treatment.prediction)]

                if (delta.value>0) {
                    mat[, (paste0("dhaz", delta.value)):=fit.dLambda]
                    mat[, (paste0("fit.cox", delta.value)):=1]
                } else {
                    mat[, surv.C.pois:=exp(-cumsum(fit.dLambda)), by=c("id", treatment.prediction)]
                    mat[, surv.C1.pois:=c(1, surv.C.pois[-.N]), by=c("id", treatment.prediction)]
                    mat[, Ht:=Ht*surv.C1/surv.C1.pois]
                    mat[, surv.C1:=surv.C1.pois]
                }
                
                if (return.cve) {
                    return(list(mat=mat, cve=cve.hal))
                } else {
                    return(mat)
                }
                
            } else {

                for (a in intervention) {
                    pseudo.dt[, (paste0("fit.dLambda", a)):=get(paste0("fit.lambda", a))*risk.time]
                    pseudo.dt[, (paste0("fit.Lambda", a)):=cumsum(get((paste0("fit.dLambda", a)))), by="id"]
                    pseudo.dt[, (paste0("surv", a)):=exp(-get(paste0("fit.Lambda", a))), by="id"]
                }
                
                ## pseudo.dt[, fit.dLambda1:=fit.lambda1*risk.time]
                ## pseudo.dt[, fit.Lambda1:=cumsum(fit.dLambda1), by="id"]

                ## pseudo.dt[, fit.dLambda0:=fit.lambda0*risk.time]
                ## pseudo.dt[, fit.Lambda0:=cumsum(fit.dLambda0), by="id"]

                ## pseudo.dt[, surv1:=exp(-fit.Lambda1), by="id"]
                ## pseudo.dt[, surv0:=exp(-fit.Lambda0), by="id"]

                if (length(intervention)>1) {
                    return(list(fit=hal.fit,
                                cve=cve.hal,
                                target.fit=pseudo.dt[time==predict, -mean(surv1 - surv0)]))
                } else {
                    return(list(fit=hal.fit,
                                cve=cve.hal,
                                target.fit=pseudo.dt[time==predict, -1-mean(get(paste0("surv", intervention)))]))
                }
            }

        } else {

            X.hal.a <- lapply(intervention, function(a) {
                return(basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=copy(dt)[,(treatment):=a],
                                 time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                                 cut.time.covar=cut.time.covar, cut.one.way.time=cut.one.way.time,
                                 predict=predict,
                                 delta.var=delta.var, delta.value=delta.value,
                                 two.way=two.way, cut.two.way=cut.two.way))
            })
            ## X.hal1 <- basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=copy(dt)[,(treatment):=1],
            ##                     time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
            ##                     predict=predict,
            ##                     delta.var=delta.var, delta.value=delta.value,
            ##                     two.way=two.way, cut.two.way=cut.two.way)
            ## X.hal0 <- basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=copy(dt)[,(treatment):=0],
            ##                     time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
            ##                     predict=predict,
            ##                     delta.var=delta.var, delta.value=delta.value,
            ##                     two.way=two.way, cut.two.way=cut.two.way)


            basehaz.a <- lapply(1:length(intervention), function(aa) {
                return(glmnet_basesurv(dt[, time], dt[, get(delta.var)==delta.value],
                                       predict(hal.fit, X.hal.a[[aa]], type="link"), centered=TRUE))
            })
            ## basehaz1 <- glmnet_basesurv(dt[, time], dt[, get(delta.var)==delta.value],
            ##                             predict(hal.fit, X.hal1, type="link"), centered=TRUE)
            ## basehaz0 <- glmnet_basesurv(dt[, time], dt[, get(delta.var)==delta.value],
            ##                             predict(hal.fit, X.hal0, type="link"), centered=TRUE)

            basehaz <-  data.table(time=c(0, basehaz.a[[1]]$time))

            for (aa in 1:length(intervention)) {
                basehaz[, (paste0("hazard", intervention[aa])):=c(0, basehaz.a[[aa]]$cumulative_base_hazard)]
            }

            basehaz <- basehaz[time<=predict]
            
            ## basehaz <- data.table(time=c(0, basehaz1$time),
            ##                       hazard1=c(0, basehaz1$cumulative_base_hazard),
            ##                       hazard0=c(0, basehaz0$cumulative_base_hazard))[time<=predict]
            
            if (length(mat)>0) {

                mat <- do.call("rbind", lapply(intervention, function(a) {
                    return(merge(mat[get(treatment)==a], basehaz,
                                 by="time", all.x=TRUE)[, (paste0("chaz", delta.value)):=get(paste0("hazard", a))][, (paste0("hazard", a)):=NULL])
                }))
                ## mat <- rbind(merge(mat[get(treatment)==1], basehaz,
                ##                    by="time", all.x=TRUE),
                ##              merge(mat[get(treatment)==0], basehaz,
                ##                    by="time", all.x=TRUE))

                ## mat[get(treatment)==1, (paste0("chaz", delta.value)):=hazard1]
                ## mat[get(treatment)==0, (paste0("chaz", delta.value)):=hazard0]

                ## mat[, hazard1:=NULL]
                ## mat[, hazard0:=NULL]

                mat[, (paste0("chaz", delta.value)):=na.locf(get(paste0("chaz", delta.value))), by=c("id", treatment)]

                mat[, (paste0("dhaz", delta.value)):=
                          c(0, diff(get(paste0("chaz", delta.value)))), by=c("id", treatment)]

                for (aa in 1:length(intervention)) {
                    dt[, (paste0("fit.hal", intervention[aa],".", delta.value)):=exp(predict(hal.fit, newx=X.hal.a[[aa]], type="link", s=lambda.cv))]
                }
                ## dt[, (paste0("fit.hal1.", delta.value)):=exp(predict(hal.fit, newx=X.hal1, type="link", s=lambda.cv))]
                ## dt[, (paste0("fit.hal0.", delta.value)):=exp(predict(hal.fit, newx=X.hal0, type="link", s=lambda.cv))]

                mat <- do.call("rbind", lapply(intervention, function(a) {
                    return(merge(mat[get(treatment)==a, -paste0("fit.cox", delta.value), with=FALSE],
                                 dt[, (paste0("fit.cox", delta.value)):=
                                          get(paste0("fit.hal", a,".", delta.value))][, c("id", paste0("fit.cox", delta.value)),
                                                                                      with=FALSE],
                                 by="id", all.x=TRUE))
                }))
                ## mat <- rbind(merge(mat[get(treatment)==1, -paste0("fit.cox", delta.value), with=FALSE],
                ##                    dt[, (paste0("fit.cox", delta.value)):=
                ##                             get(paste0("fit.hal1.", delta.value))][, c("id", paste0("fit.cox", delta.value)),
                ##                                                                    with=FALSE],
                ##                    by="id", all.x=TRUE),
                ##              merge(mat[get(treatment)==0, -paste0("fit.cox", delta.value), with=FALSE],
                ##                    dt[, (paste0("fit.cox", delta.value)):=
                ##                             get(paste0("fit.hal0.", delta.value))][, c("id", paste0("fit.cox", delta.value)),
                ##                                                                    with=FALSE],
                ##                    by="id", all.x=TRUE))

                if (delta.value==0) {
                    mat[, "surv.C.hal":=
                              exp(-get(paste0("chaz", delta.value))*
                                  get(paste0("fit.cox", delta.value)))]
                    mat[, "surv.C1.hal":=c(1, surv.C.hal[-.N]), by=c("id", treatment)]
                    mat[, Ht:=Ht*surv.C1/surv.C1.hal]
                    mat[, surv.C1:=surv.C1.hal]
                }

                if (return.cve) {
                    return(list(mat=mat, cve=cve.hal))
                } else {
                    return(mat)
                }
                
            }

            if (length(intervention)>1) {
                return(list(fit=hal.fit,
                            cve=cve.hal,
                            target.fit=-mean(exp(-exp(predict(hal.fit, type="link", newx=X.hal1, s=lambda.cv))*basehaz[nrow(basehaz), hazard1])-
                                             exp(-exp(predict(hal.fit, type="link", newx=X.hal0, s=lambda.cv))*basehaz[nrow(basehaz), hazard0]))))
            } else {
                return(list(fit=hal.fit,
                            cve=cve.hal,
                            target.fit=1-mean(exp(-exp(predict(hal.fit, type="link", newx=X.hal.a[[1]], s=lambda.cv))*basehaz[nrow(basehaz), get(paste0("hazard", intervention))]))))
            }
            
        }

    }
    
    return(list(fit=hal.fit,
                cve=cve.hal))
}
