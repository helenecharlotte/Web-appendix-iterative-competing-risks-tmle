### estimate.weibulls.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Jul 14 2022 (11:51) 
## Version: 
## Last-Updated: Jul 14 2022 (13:50) 
##           By: Helene
##     Update #: 9
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

######################################################################

#--- try see what super learner would pick
bhaz.sl <- contmle(follic, estimation=list("outcome"=list(fit="sl",
                                                          model=Surv(time, status==1)~chemo+stage+hgb+age,
                                                          lambda.cvs=seq(0.008, 0.02, length=10)),
                                           "cens"=list(fit="sl",
                                                       model=Surv(time, status==0)~chemo+stage+hgb+age),
                                           "cr2"=list(fit="sl",
                                                      model=Surv(time, status==2)~chemo+stage+hgb+age)
                                           ),
                   treat.model=chemo~stage,
                   treat.effect="ate",
                   no.small.steps=500,
                   sl.models=list(mod1=list(Surv(time, status==1)~chemo+stage+hgb+age, t0 = (1:50)/2000)), 
                   output.km=TRUE,
                   output.bhaz=TRUE, 
                   V=3, lambda.cvs=seq(0.1, 0.03, length=10), maxit=1e5, penalize.time=FALSE,
                   verbose=TRUE,
                   iterative=TRUE,
                   tau=20, target=1)

#--- uninformative censoring
bhaz.uninformative.cens <-
    contmle(follic, estimation=list("outcome"=list(fit="cox",
                                                   model=Surv(time, status==1)~chemo+stage+hgb+age, lambda.cvs=seq(0.008,
                                                                                                                   0.02, length=10)), "cens"=list(fit="cox", model=Surv(time,
                                                                                                                                                                        status==0)~1), "cr2"=list(fit="cox", model=Surv(time,
                                                                                                                                                                                                                        status==2)~chemo+stage+hgb+age) ), treat.model=chemo~stage,
            treat.effect="ate", no.small.steps=500,
            sl.models=list(mod1=list(Surv(time, status==1)~chemo+stage+hgb+age, t0
                                     = (1:50)/2000)), output.km=TRUE, output.bhaz=TRUE, V=3,
            lambda.cvs=seq(0.1, 0.03, length=10), maxit=1e5, penalize.time=FALSE,
            verbose=TRUE, iterative=TRUE, tau=20, target=1)


#######################################################################################

bhazs <- bhaz.sl[[1]]

bhazs[, chaz1 := cumsum(dhaz1*exp1), by = "chemo"]
bhazs[, chaz2 := cumsum(dhaz2*exp2), by = "chemo"] 
bhazs[, chaz0 := cumsum(dhaz0*exp0), by = "chemo"]
 
#######################################################################################
#
# Cause one events (add changepoint)

log.t0.1 <- -2#-0.5
log.t1.1 <- -0.5#0.75
log.t2.1 <- 2.5#2
log.t3.1 <- 3.15
log.t4.1 <- 3.5#3

kmin.1.t1.1 <- min((1:nrow(bhazs[chaz1>0 & chemo==1]))[log(bhazs[chaz1>0 & chemo==1][["time"]])>log.t0.1])
kmax.1.t1.1 <- max((1:nrow(bhazs[chaz1>0 & chemo==1]))[log(bhazs[chaz1>0 & chemo==1][["time"]])<log.t1.1])
kmin.1.t1.0 <- min((1:nrow(bhazs[chaz1>0 & chemo==0]))[log(bhazs[chaz1>0 & chemo==0][["time"]])>log.t0.1])
kmax.1.t1.0 <- max((1:nrow(bhazs[chaz1>0 & chemo==0]))[log(bhazs[chaz1>0 & chemo==0][["time"]])<log.t1.1])

kmin.1.t2.1 <- min((1:nrow(bhazs[chaz1>0 & chemo==1]))[log(bhazs[chaz1>0 & chemo==1][["time"]])>log.t1.1])
kmax.1.t2.1 <- max((1:nrow(bhazs[chaz1>0 & chemo==1]))[log(bhazs[chaz1>0 & chemo==1][["time"]])<log.t2.1])
kmin.1.t2.0 <- min((1:nrow(bhazs[chaz1>0 & chemo==0]))[log(bhazs[chaz1>0 & chemo==0][["time"]])>log.t1.1])
kmax.1.t2.0 <- max((1:nrow(bhazs[chaz1>0 & chemo==0]))[log(bhazs[chaz1>0 & chemo==0][["time"]])<log.t2.1])

kmin.1.t3.1 <- min((1:nrow(bhazs[chaz0>0 & chemo==1]))[log(bhazs[chaz0>0 & chemo==1][["time"]])>log.t2.1])
kmax.1.t3.1 <- max((1:nrow(bhazs[chaz0>0 & chemo==1]))[log(bhazs[chaz0>0 & chemo==1][["time"]])<log.t3.1])
kmin.1.t3.0 <- min((1:nrow(bhazs[chaz0>0 & chemo==0]))[log(bhazs[chaz0>0 & chemo==0][["time"]])>log.t2.1])
kmax.1.t3.0 <- max((1:nrow(bhazs[chaz0>0 & chemo==0]))[log(bhazs[chaz0>0 & chemo==0][["time"]])<log.t3.1])

kmin.1.t4.1 <- min((1:nrow(bhazs[chaz0>0 & chemo==1]))[log(bhazs[chaz0>0 & chemo==1][["time"]])>log.t3.1])
kmax.1.t4.1 <- max((1:nrow(bhazs[chaz0>0 & chemo==1]))[log(bhazs[chaz0>0 & chemo==1][["time"]])<log.t4.1])
kmin.1.t4.0 <- min((1:nrow(bhazs[chaz0>0 & chemo==0]))[log(bhazs[chaz0>0 & chemo==0][["time"]])>log.t3.1])
kmax.1.t4.0 <- max((1:nrow(bhazs[chaz0>0 & chemo==0]))[log(bhazs[chaz0>0 & chemo==0][["time"]])<log.t4.1])

#plot(log(bhazs[chaz1>0 & chemo==1][["time"]])[kmin.1.t1.1:kmax.1.t1.1],log(bhazs[chaz1>0 & chemo==1][["chaz1"]][kmin.1.t1.1:kmax.1.t1.1]))
fit.status1.t1.1 <- lm(log(chaz1)~log(time), data=bhazs[chaz1>0 & chemo==1][kmin.1.t1.1:kmax.1.t1.1])
#abline(a = coef(fit.status1.t1.1)[1], b = coef(fit.status1.t1.1)[2], col = "red")
(gamma.status1.t1.1 <- coef(fit.status1.t1.1)[2])
(lambda.status1.t1.1 <- exp(coef(fit.status1.t1.1)[1]/gamma.status1.t1.1))
#plot(log(bhazs[chaz1>0 & chemo==0][["time"]])[kmin.1.t1.0:kmax.1.t1.0],log(bhazs[chaz1>0 & chemo==0][["chaz1"]][kmin.1.t1.0:kmax.1.t1.0]))
fit.status1.t1.0 <- lm(log(chaz1)~log(time), data=bhazs[chaz1>0 & chemo==0][kmin.1.t1.0:kmax.1.t1.0])
#abline(a = coef(fit.status1.t1.0)[1], b = coef(fit.status1.t1.0)[2], col = "red")
(gamma.status1.t1.0 <- coef(fit.status1.t1.0)[2])
(lambda.status1.t1.0 <- exp(coef(fit.status1.t1.0)[1]/gamma.status1.t1.0))

# plot(log(bhazs[chaz1>0 & chemo==1][["time"]])[kmin.1.t2.1:kmax.1.t2.1],log(bhazs[chaz1>0 & chemo==1][["chaz1"]][kmin.1.t2.1:kmax.1.t2.1]))
fit.status1.t2.1 <- lm(log(chaz1)~log(time), data=bhazs[chaz1>0 & chemo==1][kmin.1.t2.1:kmax.1.t2.1])
#abline(a = coef(fit.status1.t2.1)[1], b = coef(fit.status1.t2.1)[2], col = "red")
(gamma.status1.t2.1 <- coef(fit.status1.t2.1)[2])
(lambda.status1.t2.1 <- exp(coef(fit.status1.t2.1)[1]/gamma.status1.t2.1))
#plot(log(bhazs[chaz1>0 & chemo==0][["time"]])[kmin.1.t2.0:kmax.1.t2.0],log(bhazs[chaz1>0 & chemo==0][["chaz1"]][kmin.1.t2.0:kmax.1.t2.0]))
fit.status1.t2.0 <- lm(log(chaz1)~log(time), data=bhazs[chaz1>0 & chemo==0][kmin.1.t2.0:kmax.1.t2.0])
#abline(a = coef(fit.status1.t2.0)[1], b = coef(fit.status1.t2.0)[2], col = "red")
(gamma.status1.t2.0 <- coef(fit.status1.t2.0)[2])
(lambda.status1.t2.0 <- exp(coef(fit.status1.t2.0)[1]/gamma.status1.t2.0))

# plot(log(bhazs[chaz1>0 & chemo==1][["time"]])[kmin.1.t3.1:kmax.1.t3.1],log(bhazs[chaz1>0 & chemo==1][["chaz1"]][kmin.1.t3.1:kmax.1.t3.1]))
fit.status1.t3.1 <- lm(log(chaz1)~log(time), data=bhazs[chaz1>0 & chemo==1][kmin.1.t3.1:kmax.1.t3.1])
#abline(a = coef(fit.status1.t3.1)[1], b = coef(fit.status1.t3.1)[2], col = "red")
(gamma.status1.t3.1 <- coef(fit.status1.t3.1)[2])
(lambda.status1.t3.1 <- exp(coef(fit.status1.t3.1)[1]/gamma.status1.t3.1))
#plot(log(bhazs[chaz1>0 & chemo==0][["time"]])[kmin.1.t3.0:kmax.1.t3.0],log(bhazs[chaz1>0 & chemo==0][["chaz1"]][kmin.1.t3.0:kmax.1.t3.0]))
fit.status1.t3.0 <- lm(log(chaz1)~log(time), data=bhazs[chaz1>0 & chemo==0][kmin.1.t3.0:kmax.1.t3.0])
#abline(a = coef(fit.status1.t3.0)[1], b = coef(fit.status1.t3.0)[2], col = "red")
(gamma.status1.t3.0 <- coef(fit.status1.t3.0)[2])
(lambda.status1.t3.0 <- exp(coef(fit.status1.t3.0)[1]/gamma.status1.t3.0))

# plot(log(bhazs[chaz1>0 & chemo==1][["time"]])[kmin.1.t4.1:kmax.1.t4.1],log(bhazs[chaz1>0 & chemo==1][["chaz1"]][kmin.1.t4.1:kmax.1.t4.1]))
fit.status1.t4.1 <- lm(log(chaz1)~log(time), data=bhazs[chaz1>0 & chemo==1][kmin.1.t4.1:kmax.1.t4.1])
#abline(a = coef(fit.status1.t4.1)[1], b = coef(fit.status1.t4.1)[2], col = "red")
(gamma.status1.t4.1 <- coef(fit.status1.t4.1)[2])
(lambda.status1.t4.1 <- exp(coef(fit.status1.t4.1)[1]/gamma.status1.t4.1))
#plot(log(bhazs[chaz1>0 & chemo==0][["time"]])[kmin.1.t4.0:kmax.1.t4.0],log(bhazs[chaz1>0 & chemo==0][["chaz1"]][kmin.1.t4.0:kmax.1.t4.0]))
fit.status1.t4.0 <- lm(log(chaz1)~log(time), data=bhazs[chaz1>0 & chemo==0][kmin.1.t4.0:kmax.1.t4.0])
#abline(a = coef(fit.status1.t4.0)[1], b = coef(fit.status1.t4.0)[2], col = "red")
(gamma.status1.t4.0 <- coef(fit.status1.t4.0)[2])
(lambda.status1.t4.0 <- exp(coef(fit.status1.t4.0)[1]/gamma.status1.t4.0))

#######################################################################################
#
# Cause two events (add changepoint)

log.t0.2 <- 0 
log.t1.2 <- 1.5 
log.t2.2 <- 3
log.t3.2 <- 3.6  

kmin.2.t1.1 <- min((1:nrow(bhazs[chaz2>0 & chemo==1]))[log(bhazs[chaz2>0 & chemo==1][["time"]])>log.t0.2])
kmax.2.t1.1 <- max((1:nrow(bhazs[chaz2>0 & chemo==1]))[log(bhazs[chaz2>0 & chemo==1][["time"]])<log.t1.2])
kmin.2.t1.0 <- min((1:nrow(bhazs[chaz2>0 & chemo==0]))[log(bhazs[chaz2>0 & chemo==0][["time"]])>log.t0.2])
kmax.2.t1.0 <- max((1:nrow(bhazs[chaz2>0 & chemo==0]))[log(bhazs[chaz2>0 & chemo==0][["time"]])<log.t1.2])

kmin.2.t2.1 <- min((1:nrow(bhazs[chaz2>0 & chemo==1]))[log(bhazs[chaz2>0 & chemo==1][["time"]])>log.t1.2])
kmax.2.t2.1 <- max((1:nrow(bhazs[chaz2>0 & chemo==1]))[log(bhazs[chaz2>0 & chemo==1][["time"]])<log.t2.2])
kmin.2.t2.0 <- min((1:nrow(bhazs[chaz2>0 & chemo==0]))[log(bhazs[chaz2>0 & chemo==0][["time"]])>log.t1.2])
kmax.2.t2.0 <- max((1:nrow(bhazs[chaz2>0 & chemo==0]))[log(bhazs[chaz2>0 & chemo==0][["time"]])<log.t2.2])

kmin.2.t3.1 <- min((1:nrow(bhazs[chaz2>0 & chemo==1]))[log(bhazs[chaz2>0 & chemo==1][["time"]])>log.t2.2])
kmax.2.t3.1 <- max((1:nrow(bhazs[chaz2>0 & chemo==1]))[log(bhazs[chaz2>0 & chemo==1][["time"]])<log.t3.2])
kmin.2.t3.0 <- min((1:nrow(bhazs[chaz2>0 & chemo==0]))[log(bhazs[chaz2>0 & chemo==0][["time"]])>log.t2.2])
kmax.2.t3.0 <- max((1:nrow(bhazs[chaz2>0 & chemo==0]))[log(bhazs[chaz2>0 & chemo==0][["time"]])<log.t3.2])

#plot(log(bhazs[chaz2>0 & chemo==1][["time"]])[kmin.2.t1.1:kmax.2.t1.1],log(bhazs[chaz2>0 & chemo==1][["chaz2"]][kmin.2.t1.1:kmax.2.t1.1]))
fit.status2.t1.1 <- lm(log(chaz2)~log(time), data=bhazs[chaz2>0 & chemo==1][kmin.2.t1.1:kmax.2.t1.1])
#abline(a = coef(fit.status2.t1.1)[1], b = coef(fit.status2.t1.1)[2], col = "red")
(gamma.status2.t1.1 <- coef(fit.status2.t1.1)[2])
(lambda.status2.t1.1 <- exp(coef(fit.status2.t1.1)[1]/gamma.status2.t1.1))
#plot(log(bhazs[chaz2>0 & chemo==0][["time"]])[kmin.2.t1.0:kmax.2.t1.0],log(bhazs[chaz2>0 & chemo==0][["chaz2"]][kmin.2.t1.0:kmax.2.t1.0]))
fit.status2.t1.0 <- lm(log(chaz2)~log(time), data=bhazs[chaz2>0 & chemo==0][kmin.2.t1.0:kmax.2.t1.0])
#abline(a = coef(fit.status2.t1.0)[1], b = coef(fit.status2.t1.0)[2], col = "red")
(gamma.status2.t1.0 <- coef(fit.status2.t1.0)[2])
(lambda.status2.t1.0 <- exp(coef(fit.status2.t1.0)[1]/gamma.status2.t1.0))

# plot(log(bhazs[chaz2>0 & chemo==1][["time"]])[kmin.2.t2.1:kmax.2.t2.1],log(bhazs[chaz2>0 & chemo==1][["chaz2"]][kmin.2.t2.1:kmax.2.t2.1]))
fit.status2.t2.1 <- lm(log(chaz2)~log(time), data=bhazs[chaz2>0 & chemo==1][kmin.2.t2.1:kmax.2.t2.1])
#abline(a = coef(fit.status2.t2.1)[1], b = coef(fit.status2.t2.1)[2], col = "red")
(gamma.status2.t2.1 <- coef(fit.status2.t2.1)[2])
(lambda.status2.t2.1 <- exp(coef(fit.status2.t2.1)[1]/gamma.status2.t2.1))
#plot(log(bhazs[chaz2>0 & chemo==0][["time"]])[kmin.2.t2.0:kmax.2.t2.0],log(bhazs[chaz2>0 & chemo==0][["chaz2"]][kmin.2.t2.0:kmax.2.t2.0]))
fit.status2.t2.0 <- lm(log(chaz2)~log(time), data=bhazs[chaz2>0 & chemo==0][kmin.2.t2.0:kmax.2.t2.0])
#abline(a = coef(fit.status2.t2.0)[1], b = coef(fit.status2.t2.0)[2], col = "red")
(gamma.status2.t2.0 <- coef(fit.status2.t2.0)[2])
(lambda.status2.t2.0 <- exp(coef(fit.status2.t2.0)[1]/gamma.status2.t2.0))

# plot(log(bhazs[chaz2>0 & chemo==1][["time"]])[kmin.2.t3.1:kmax.2.t3.1],log(bhazs[chaz2>0 & chemo==1][["chaz2"]][kmin.2.t3.1:kmax.2.t3.1]))
fit.status2.t3.1 <- lm(log(chaz2)~log(time), data=bhazs[chaz2>0 & chemo==1][kmin.2.t3.1:kmax.2.t3.1])
#abline(a = coef(fit.status2.t3.1)[1], b = coef(fit.status2.t3.1)[2], col = "red")
(gamma.status2.t3.1 <- coef(fit.status2.t3.1)[2])
(lambda.status2.t3.1 <- exp(coef(fit.status2.t3.1)[1]/gamma.status2.t3.1))
#plot(log(bhazs[chaz2>0 & chemo==0][["time"]])[kmin.2.t3.0:kmax.2.t3.0],log(bhazs[chaz2>0 & chemo==0][["chaz2"]][kmin.2.t3.0:kmax.2.t3.0]))
fit.status2.t3.0 <- lm(log(chaz2)~log(time), data=bhazs[chaz2>0 & chemo==0][kmin.2.t3.0:kmax.2.t3.0])
#abline(a = coef(fit.status2.t3.0)[1], b = coef(fit.status2.t3.0)[2], col = "red")
(gamma.status2.t3.0 <- coef(fit.status2.t3.0)[2])
(lambda.status2.t3.0 <- exp(coef(fit.status2.t3.0)[1]/gamma.status2.t3.0))

#######################################################################################
#
# Censoring events (covariate dependent, add changepoint)

log.t0.0 <- 1.2
log.t1.0 <- 1.9   
log.t2.0 <- 3.2
log.t3.0 <- 3.5

kmin.0.t1.1 <- min((1:nrow(bhazs[chaz0>0 & chemo==1]))[log(bhazs[chaz0>0 & chemo==1][["time"]])>log.t0.0])
kmax.0.t1.1 <- max((1:nrow(bhazs[chaz0>0 & chemo==1]))[log(bhazs[chaz0>0 & chemo==1][["time"]])<log.t1.0])
kmin.0.t1.0 <- min((1:nrow(bhazs[chaz0>0 & chemo==0]))[log(bhazs[chaz0>0 & chemo==0][["time"]])>log.t0.0])
kmax.0.t1.0 <- max((1:nrow(bhazs[chaz0>0 & chemo==0]))[log(bhazs[chaz0>0 & chemo==0][["time"]])<log.t1.0])

kmin.0.t2.1 <- min((1:nrow(bhazs[chaz0>0 & chemo==1]))[log(bhazs[chaz0>0 & chemo==1][["time"]])>log.t1.0])
kmax.0.t2.1 <- max((1:nrow(bhazs[chaz0>0 & chemo==1]))[log(bhazs[chaz0>0 & chemo==1][["time"]])<log.t2.0])
kmin.0.t2.0 <- min((1:nrow(bhazs[chaz0>0 & chemo==0]))[log(bhazs[chaz0>0 & chemo==0][["time"]])>log.t1.0])
kmax.0.t2.0 <- max((1:nrow(bhazs[chaz0>0 & chemo==0]))[log(bhazs[chaz0>0 & chemo==0][["time"]])<log.t2.0])

kmin.0.t3.1 <- min((1:nrow(bhazs[chaz0>0 & chemo==1]))[log(bhazs[chaz0>0 & chemo==1][["time"]])>log.t2.0])
kmax.0.t3.1 <- max((1:nrow(bhazs[chaz0>0 & chemo==1]))[log(bhazs[chaz0>0 & chemo==1][["time"]])<log.t3.0])
kmin.0.t3.0 <- min((1:nrow(bhazs[chaz0>0 & chemo==0]))[log(bhazs[chaz0>0 & chemo==0][["time"]])>log.t2.0])
kmax.0.t3.0 <- max((1:nrow(bhazs[chaz0>0 & chemo==0]))[log(bhazs[chaz0>0 & chemo==0][["time"]])<log.t3.0])

#plot(log(bhazs[chaz0>0 & chemo==1][["time"]])[kmin.0.t1.1:kmax.0.t1.1],log(bhazs[chaz0>0 & chemo==1][["chaz0"]][kmin.0.t1.1:kmax.0.t1.1]))
fit.status0.t1.1 <- lm(log(chaz0)~log(time), data=bhazs[chaz0>0 & chemo==1][kmin.0.t1.1:kmax.0.t1.1])
#abline(a = coef(fit.status0.t1.1)[1], b = coef(fit.status0.t1.1)[2], col = "red")
(gamma.status0.t1.1 <- coef(fit.status0.t1.1)[2])
(lambda.status0.t1.1 <- exp(coef(fit.status0.t1.1)[1]/gamma.status0.t1.1))
#plot(log(bhazs[chaz0>0 & chemo==0][["time"]])[kmin.0.t1.0:kmax.0.t1.0],log(bhazs[chaz0>0 & chemo==0][["chaz0"]][kmin.0.t1.0:kmax.0.t1.0]))
fit.status0.t1.0 <- lm(log(chaz0)~log(time), data=bhazs[chaz0>0 & chemo==0][kmin.0.t1.0:kmax.0.t1.0])
#abline(a = coef(fit.status0.t1.0)[1], b = coef(fit.status0.t1.0)[2], col = "red")
(gamma.status0.t1.0 <- coef(fit.status0.t1.0)[2])
(lambda.status0.t1.0 <- exp(coef(fit.status0.t1.0)[1]/gamma.status0.t1.0))

# plot(log(bhazs[chaz0>0 & chemo==1][["time"]])[kmin.0.t2.1:kmax.0.t2.1],log(bhazs[chaz0>0 & chemo==1][["chaz0"]][kmin.0.t2.1:kmax.0.t2.1]))
fit.status0.t2.1 <- lm(log(chaz0)~log(time), data=bhazs[chaz0>0 & chemo==1][kmin.0.t2.1:kmax.0.t2.1])
#abline(a = coef(fit.status0.t2.1)[1], b = coef(fit.status0.t2.1)[2], col = "red")
(gamma.status0.t2.1 <- coef(fit.status0.t2.1)[2])
(lambda.status0.t2.1 <- exp(coef(fit.status0.t2.1)[1]/gamma.status0.t2.1))
#plot(log(bhazs[chaz0>0 & chemo==0][["time"]])[kmin.0.t2.0:kmax.0.t2.0],log(bhazs[chaz0>0 & chemo==0][["chaz0"]][kmin.0.t2.0:kmax.0.t2.0]))
fit.status0.t2.0 <- lm(log(chaz0)~log(time), data=bhazs[chaz0>0 & chemo==0][kmin.0.t2.0:kmax.0.t2.0])
#abline(a = coef(fit.status0.t2.0)[1], b = coef(fit.status0.t2.0)[2], col = "red")
(gamma.status0.t2.0 <- coef(fit.status0.t2.0)[2])
(lambda.status0.t2.0 <- exp(coef(fit.status0.t2.0)[1]/gamma.status0.t2.0))

# plot(log(bhazs[chaz0>0 & chemo==1][["time"]])[kmin.0.t3.1:kmax.0.t3.1],log(bhazs[chaz0>0 & chemo==1][["chaz0"]][kmin.0.t3.1:kmax.0.t3.1]))
fit.status0.t3.1 <- lm(log(chaz0)~log(time), data=bhazs[chaz0>0 & chemo==1][kmin.0.t3.1:kmax.0.t3.1])
#abline(a = coef(fit.status0.t3.1)[1], b = coef(fit.status0.t3.1)[2], col = "red")
(gamma.status0.t3.1 <- coef(fit.status0.t3.1)[2])
(lambda.status0.t3.1 <- exp(coef(fit.status0.t3.1)[1]/gamma.status0.t3.1))
#plot(log(bhazs[chaz0>0 & chemo==0][["time"]])[kmin.0.t3.0:kmax.0.t3.0],log(bhazs[chaz0>0 & chemo==0][["chaz0"]][kmin.0.t3.0:kmax.0.t3.0]))
fit.status0.t3.0 <- lm(log(chaz0)~log(time), data=bhazs[chaz0>0 & chemo==0][kmin.0.t3.0:kmax.0.t3.0])
#abline(a = coef(fit.status0.t3.0)[1], b = coef(fit.status0.t3.0)[2], col = "red")
(gamma.status0.t3.0 <- coef(fit.status0.t3.0)[2])
(lambda.status0.t3.0 <- exp(coef(fit.status0.t3.0)[1]/gamma.status0.t3.0))

#######################################################################################
#
# Censoring events (independent, add changepoint)

bhazs.uninformative.cens <- bhaz.uninformative.cens[[1]]
bhazs.uninformative.cens[, chaz0 := cumsum(dhaz0*exp0), by = "chemo"] 

bhazs.uninformative.cens.long <- melt(bhazs.uninformative.cens, id.vars=c("chemo", "time")) 
bhazs.uninformative.cens.long[, variable2:=substr(variable,1,4)] 
bhazs.uninformative.cens.long <- bhazs.uninformative.cens.long[variable2=="chaz"][, status:=paste0("status = ", gsub("chaz", "", variable))]
bhazs.uninformative.cens.long[, chemo:=paste0("chemo = ", chemo)]

log.t0.0 <- 1.2
log.t1.0 <- 1.9   
log.t2.0 <- 3.2
log.t3.0 <- 3.5

kmin.0.t1.1 <- min((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==1]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["time"]])>log.t0.0])
kmax.0.t1.1 <- max((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==1]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["time"]])<log.t1.0])
kmin.0.t1.0 <- min((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==0]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["time"]])>log.t0.0])
kmax.0.t1.0 <- max((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==0]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["time"]])<log.t1.0])

kmin.0.t2.1 <- min((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==1]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["time"]])>log.t1.0])
kmax.0.t2.1 <- max((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==1]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["time"]])<log.t2.0])
kmin.0.t2.0 <- min((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==0]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["time"]])>log.t1.0])
kmax.0.t2.0 <- max((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==0]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["time"]])<log.t2.0])

kmin.0.t3.1 <- min((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==1]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["time"]])>log.t2.0])
kmax.0.t3.1 <- max((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==1]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["time"]])<log.t3.0])
kmin.0.t3.0 <- min((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==0]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["time"]])>log.t2.0])
kmax.0.t3.0 <- max((1:nrow(bhazs.uninformative.cens[chaz0>0 & chemo==0]))[log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["time"]])<log.t3.0])

#plot(log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["time"]])[kmin.0.t1.1:kmax.0.t1.1],log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["chaz0"]][kmin.0.t1.1:kmax.0.t1.1]))
fit.status0.independent.t1.1 <- lm(log(chaz0)~log(time), data=bhazs.uninformative.cens[chaz0>0 & chemo==1][kmin.0.t1.1:kmax.0.t1.1])
#abline(a = coef(fit.status0.independent.t1.1)[1], b = coef(fit.status0.independent.t1.1)[2], col = "red")
(gamma.status0.independent.t1.1 <- coef(fit.status0.independent.t1.1)[2])
(lambda.status0.independent.t1.1 <- exp(coef(fit.status0.independent.t1.1)[1]/gamma.status0.independent.t1.1))
#plot(log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["time"]])[kmin.0.t1.0:kmax.0.t1.0],log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["chaz0"]][kmin.0.t1.0:kmax.0.t1.0]))
fit.status0.independent.t1.0 <- lm(log(chaz0)~log(time), data=bhazs.uninformative.cens[chaz0>0 & chemo==0][kmin.0.t1.0:kmax.0.t1.0])
#abline(a = coef(fit.status0.independent.t1.0)[1], b = coef(fit.status0.independent.t1.0)[2], col = "red")
(gamma.status0.independent.t1.0 <- coef(fit.status0.independent.t1.0)[2])
(lambda.status0.independent.t1.0 <- exp(coef(fit.status0.independent.t1.0)[1]/gamma.status0.independent.t1.0))

# plot(log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["time"]])[kmin.0.t2.1:kmax.0.t2.1],log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["chaz0"]][kmin.0.t2.1:kmax.0.t2.1]))
fit.status0.independent.t2.1 <- lm(log(chaz0)~log(time), data=bhazs.uninformative.cens[chaz0>0 & chemo==1][kmin.0.t2.1:kmax.0.t2.1])
#abline(a = coef(fit.status0.independent.t2.1)[1], b = coef(fit.status0.independent.t2.1)[2], col = "red")
(gamma.status0.independent.t2.1 <- coef(fit.status0.independent.t2.1)[2])
(lambda.status0.independent.t2.1 <- exp(coef(fit.status0.independent.t2.1)[1]/gamma.status0.independent.t2.1))
#plot(log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["time"]])[kmin.0.t2.0:kmax.0.t2.0],log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["chaz0"]][kmin.0.t2.0:kmax.0.t2.0]))
fit.status0.independent.t2.0 <- lm(log(chaz0)~log(time), data=bhazs.uninformative.cens[chaz0>0 & chemo==0][kmin.0.t2.0:kmax.0.t2.0])
#abline(a = coef(fit.status0.independent.t2.0)[1], b = coef(fit.status0.independent.t2.0)[2], col = "red")
(gamma.status0.independent.t2.0 <- coef(fit.status0.independent.t2.0)[2])
(lambda.status0.independent.t2.0 <- exp(coef(fit.status0.independent.t2.0)[1]/gamma.status0.independent.t2.0))

# plot(log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["time"]])[kmin.0.t3.1:kmax.0.t3.1],log(bhazs.uninformative.cens[chaz0>0 & chemo==1][["chaz0"]][kmin.0.t3.1:kmax.0.t3.1]))
fit.status0.independent.t3.1 <- lm(log(chaz0)~log(time), data=bhazs.uninformative.cens[chaz0>0 & chemo==1][kmin.0.t3.1:kmax.0.t3.1])
#abline(a = coef(fit.status0.independent.t3.1)[1], b = coef(fit.status0.independent.t3.1)[2], col = "red")
(gamma.status0.independent.t3.1 <- coef(fit.status0.independent.t3.1)[2])
(lambda.status0.independent.t3.1 <- exp(coef(fit.status0.independent.t3.1)[1]/gamma.status0.independent.t3.1))
#plot(log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["time"]])[kmin.0.t3.0:kmax.0.t3.0],log(bhazs.uninformative.cens[chaz0>0 & chemo==0][["chaz0"]][kmin.0.t3.0:kmax.0.t3.0]))
fit.status0.independent.t3.0 <- lm(log(chaz0)~log(time), data=bhazs.uninformative.cens[chaz0>0 & chemo==0][kmin.0.t3.0:kmax.0.t3.0])
#abline(a = coef(fit.status0.independent.t3.0)[1], b = coef(fit.status0.independent.t3.0)[2], col = "red")
(gamma.status0.independent.t3.0 <- coef(fit.status0.independent.t3.0)[2])
(lambda.status0.independent.t3.0 <- exp(coef(fit.status0.independent.t3.0)[1]/gamma.status0.independent.t3.0))

######################################################################

#--- for simulating baseline covariates
#

p.stage <- mean(follic[["stage"]] == 1)
p.age <- fitdistr(follic[["age"]], "normal")
p.hgb <- fitdistr(follic[["hgb"]], "normal")
p.chemo <- mean(follic[["chemo"]] == 1)

######################################################################
### estimate.weibulls
