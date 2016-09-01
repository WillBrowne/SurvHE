rm(list=ls())
#setwd("~/Msc/Diss/")

library(flexsurv)
library(INLA)
library(data.table)
library(R2OpenBUGS)
WINE="/opt/local/bin/wine"
WINEPATH="/opt/local/bin/winepath"
OpenBUGS.pgm="/Users/wilbrowne/.wine_new/drive_c/Program\ Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe"

source("survHE.R")

# Reads data in
dat <- read.table("http://www.statistica.it/gianluca/survHE/data.txt",header=TRUE)
# Adds some fictious covariates - for testing
dat$sex <- rbinom(dim(dat)[1],1,.5)
dat$age <- rpois(dim(dat)[1],32)
dat$imd <- cut(runif(dim(dat)[1],0,100),breaks=seq(0,100,20)); levels(dat$imd)=1:5
dat$ethnic <- cut(runif(dim(dat)[1],0,100),breaks=seq(0,100,20)); levels(dat$ethnic)=1:5
formula=Surv(time,censored)~as.factor(arm)+as.factor(imd)

 test <- attributes(terms(formula))$term.labels
  ncovs <- length(test)
  formula.temp <- as.formula(gsub("inla.surv","Surv",deparse(formula)))
  time <- all.vars(formula.temp,data)[1]
  event <- all.vars(formula.temp,data)[2]
if (length(ncovs)>1) {
      Xraw <- model.frame(formula.temp,data)
      w <- (which(sapply(Xraw,is.factor)==1))-1
      if (length(w)>=1) {
          factors <- gsub("as.factor[( )]","",test) 
          factors <- gsub("[( )]","",factors)
          covs <- test[-w]
          if (length(covs)==0) {
              covs <- NULL
          }
      } else {
          factors <- NULL
          covs <- test
      }
  }

    if(!is.null(covs)) {
      X <- data[,pmatch(covs,colnames(data))]
      K <- ifelse(!is.null(dim(X)[2]),dim(X)[2],1)
  }
  # If there are categorical covariates (in the vector 'factors'), makes sure they have the right form
  if(!is.null(factors)) {
      cols <- pmatch(factors,colnames(data))
      H <- length(cols)
      D <- numeric()
      for (i in 1:H) {
          data[,cols[i]] <- as.factor(data[,cols[i]])
          nlevs <- length(levels(data[,cols[i]]))
          # Now if the method is MCMC recodes the levels of the factors so that they can be used in BUGS
          if (method=="mcmc") {
              levels(data[,cols[i]]) <- 1:nlevs
          }
          D[i] <- nlevs
      }
  }
        print("get position")
        position <- pmatch(distr,availables.mcmc)
        print("write data")
        dataBugs <- write.data(position)
        print("get params")
        params <- write.params(position)
        print("write model")
        model.file <- write.model(position)$name
model <- R2OpenBUGS::bugs(data=dataBugs,inits=inits,parameters.to.save=params,
                                  model.file=model.file,n.chains=2,n.iter=10000,
                                  n.burnin=5000,n.thin=0,debug=TRUE,
                                  OpenBUGS.pgm=OpenBUGS.pgm,useWINE=TRUE,
                                  WINEPATH=WINEPATH,WINE=WINE)

x <- rgenf()
dfmod <- R2OpenBUGS::bugs(data=datadf,inits=list(inits1,inits2),parameters.to.save=params,
                                  model.file="~/Msc/Diss/genf.txt",n.chains=2,n.iter=10000,
                                  n.burnin=5000,n.thin=0,debug=TRUE,
                                  OpenBUGS.pgm=OpenBUGS.pgm,useWINE=TRUE,
                                  WINEPATH=WINEPATH,WINE=WINE)


mods.mle <- c("weibull","exp","gamma","lnorm","llogis","gengamma")   
mods.inla <- c("exp","weibull","lognormal","loglogistic")

fit.splines <- fit.models(formula=formula,data=data,knots = 4,method="splines")


fit.mle <- 
fit.mcmc <- fit.models(formula=formula,data=dat,distr="gompertz",method="mcmc",n.iter=10000,useWINE = TRUE,WINEPATH=WINEPATH,OpenBUGS.pgm = OpenBUGS.pgm,WINE = WINE)

fit.mle <- fit.models(formula=formula,data=dat,distr="genf",method="mle")

fit.mcmc <- fit.models(formula=formula,data=dat,distr="genf",method="mcmc",n.iter=10000,useWINE = TRUE,WINEPATH=WINEPATH,OpenBUGS.pgm = OpenBUGS.pgm,WINE = WINE)

fit.mle <- fit.models(formula=formula,data=dat,distr="exp",method="mle")
fit.mcmc <- fit.models(formula=formula,data=dat,distr="exp",method="mcmc",n.iter=10000,useWINE = TRUE,WINEPATH=WINEPATH,OpenBUGS.pgm = OpenBUGS.pgm,WINE = WINE)
fit.mcmc.dum <- fit.models(formula=formula,data=dat,distr="exp_d",method="mcmc",n.iter=10000,useWINE = TRUE,WINEPATH=WINEPATH,OpenBUGS.pgm = OpenBUGS.pgm,WINE = WINE)


fit.mle <- fit.models(formula=formula,data=dat,distr=mods.mle)
fit.inla <- fit.models(formula=formula,data=dat,distr=mods.inla,
                       method="inla",control.family=list(lognormal=list(initial=0)))

fit.inla <- fit.models(formula=formula,data=dat,distr="exp", method="inla")
"exp"
fit.mcmc <- fit.models(formula=formula,data=dat,distr="genf",method="mcmc",n.iter=10000,useWINE = TRUE,WINEPATH=WINEPATH,OpenBUGS.pgm = OpenBUGS.pgm,WINE = WINE)
fit.mcmc <- fit.models(formula=formula,data=dat,distr="gompz",method="mcmc",n.iter=10000,useWINE = TRUE,WINEPATH=WINEPATH,OpenBUGS.pgm = OpenBUGS.pgm,WINE = WINE)

tmp <- fit.models(formula=formula,data=dat,distr="loglogistic",method="mcmc",n.iter=10000,
  useWINE = TRUE,WINEPATH=WINEPATH,OpenBUGS.pgm = OpenBUGS.pgm,WINE = WINE)


psa <- make.surv(fit.inla,nsim = 100)

psa.plot(psa)

print.survHE(fit.mle) ## specify mod as vector

### this seems really ugly at the moment.
plot.survHE(fit.mle)

model.fit.plot(fit.mle) # seems to work well for all models.

## Create a smaller version of the data to see the influence of the priors
n.red <- 15
w <- sample(dim(dat)[1],n.red) # is this meant to be dat?

## replace dat to dat_short
dat_short <- dat[w,] # similarly

distr="exponential"

# Runs flexsurv (MLE analysis)
library(flexsurv)
formula <- as.formula("Surv(time,censored)~as.factor(arm)")
tic <- proc.time()
mle <- flexsurvreg(formula=formula,data=data,dist=distr)
toc <- proc.time()-tic
mle$time2run <- toc[3]


# Runs INLA
library(INLA)
data$treat <- as.factor(data$arm)
formula2 <- inla.surv(time,censored) ~ treat
inla <- inla(formula2,family=distr,data=data,
          # These are needed to then simulate from the joint posterior of all parameters
          control.compute=list(config=TRUE),
          control.inla=list(int.strategy="grid",dz=.1,diff.logdens=5),
          # These controls the settings for the priors 
          control.fixed=list(mean=list(treat=0),mean.intercept=0,
                             prec=list(treat=0.001),prec.intercept=0.0000001),
          control.family=list(hyper=list(theta=list(prior="loggamma",param=c(25,25))))
)

# Now rescales the coefficients to the AF parametrisation
jpost <- suppressWarnings(inla.posterior.sample(n=1000,inla))
shape <- unlist(lapply(jpost, function(x) x$hyperpar))	
scale <- exp(unlist(lapply(jpost, function(x) x$latent["(Intercept).1",])))^(1/-shape)
effect <- (exp(unlist(lapply(jpost, function(x) x$latent["treat1.1",]))))^(1/-shape)
tmp <- cbind(shape,scale,effect=log(effect)); rownames(tmp) <- NULL



# Runs JAGS (models the Weibull as a special case of Generalised Gamma distribution --- AF parameterisation)
library(R2jags)
model.string <- "model {
for (i in 1:n) {
is.censored[i] ~ dinterval(t[i],censored[i])
t[i] ~ dgen.gamma(1,lambda[i],shape)
log(lambda[i]) <- beta0 + beta1[(arm[i]+1)]
}
# Priors (same as default in INLA)
beta0 ~ dnorm(0,prec.int.inla)  
beta1[1] <- 0
beta1[2] ~ dnorm(1,prec.fe.inla)
shape ~ dgamma(shape.inla,rate.inla)
scale <- exp(-beta0)
effect <- exp(-beta1[2])
}"

dataJags <- with(data,
		list("is.censored"=(censored==0),
                     "censored"=time,
                     "t"=ifelse(censored==1,time,NA),
                     "arm"=as.numeric(arm),"n"=length(arm),
		# These induce the same priors as in the default INLA??
		      "shape.inla"=25,"rate.inla"=25,
		      "prec.fe.inla"=0.001,"prec.int.inla"=0.000001
		)
)

inits <- function(){
	list(t=ifelse(dataJags$is.censored,data$time+1,NA),
	shape=runif(1),beta0=rnorm(1),beta1=c(NA,rnorm(1)))
}

n.iter <- 10000
n.burnin <- 5000
n.samples <- 1000
n.chains <- 2
n.thin <- floor((n.iter-n.burnin)/(n.samples/n.chains))
params <- c("shape","beta1","beta0","effect","scale")
tic <- proc.time()
mcmc <- jags(data=dataJags,inits=inits,parameters.to.save=params,
              model.file=textConnection(model.string),n.chains=n.chains,
              n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin)
toc <- proc.time()-tic
mcmc$time2run <- toc[3]
tmp2 <- cbind(mcmc$BUGSoutput$sims.list$shape,mcmc$BUGSoutput$sims.list$scale,
              log(mcmc$BUGSoutput$sims.list$effect)); colnames(tmp2) <- c("shape","scale","effect")


#### Runs JAGS (directly using Weibull distribution --- PH parameterisation)
library(R2jags)
model.string <- "model {
for (i in 1:n) {
is.censored[i] ~ dinterval(t[i],censored[i])
t[i] ~ dweib(shape,lambda[i])
log(lambda[i]) <- beta0 + beta1[(arm[i]+1)]
}
# Priors (same as default in INLA)
beta0 ~ dnorm(0,prec.int.inla)  
beta1[1] <- 0
beta1[2] ~ dnorm(1,prec.fe.inla)
shape ~ dgamma(shape.inla,rate.inla)
scale <- pow(exp(beta0),1/-shape)
effect <- pow(exp(beta1[2]),1/-shape)
}"

dataJags <- with(data,
                 list("is.censored"=(censored==0),
                      "censored"=time,
                      "t"=ifelse(censored==1,time,NA),
                      "arm"=as.numeric(arm),"n"=length(arm),
                      # These induce the same priors as in the default INLA??
                      "shape.inla"=25,"rate.inla"=25,
                      "prec.fe.inla"=0.001,"prec.int.inla"=0.000001
                 )
)

inits <- function(){
  list(t=ifelse(dataJags$is.censored,data$time+1,NA),
       shape=runif(1),beta0=rnorm(1),beta1=c(NA,rnorm(1)))
}

n.iter <- 10000
n.burnin <- 5000
n.samples <- 1000
n.chains <- 2
n.thin <- floor((n.iter-n.burnin)/(n.samples/n.chains))
params <- c("shape","beta1","beta0","effect","scale")
tic <- proc.time()
mcmc2 <- jags(data=dataJags,inits=inits,parameters.to.save=params,
             model.file=textConnection(model.string),n.chains=n.chains,
             n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin)
toc <- proc.time()-tic
mcmc2$time2run <- toc[3]
tmp3 <- cbind(mcmc2$BUGSoutput$sims.list$shape,mcmc2$BUGSoutput$sims.list$scale,
              log(mcmc2$BUGSoutput$sims.list$effect)); colnames(tmp3) <- c("shape","scale","effect")



## Compares the results
# MLE
mle$res

# INLA
t(apply(tmp,2,function(x) c("mean"=mean(x),quantile(x,.025),quantile(x,.975),"sd"=sd(x))))

# JAGS (AF parameterisation)
t(apply(tmp2,2,function(x) c("mean"=mean(x),quantile(x,.025),quantile(x,.975),"sd"=sd(x))))

# JAGS (PH parameterisation)
t(apply(tmp3,2,function(x) c("mean"=mean(x),quantile(x,.025),quantile(x,.975),"sd"=sd(x))))




### Checks stuff in INLA
distr <- c("weibull","exponential","loglogistic")
control.family <- replicate(length(distr),list(inla.set.control.family.default()))
names(control.family) <- distr

# Defines the control.family parameters for just two out of the three distributions specified in distr
cf=list(weibull=list(param=c(0.1,.01)),loglogistic=list(intials=2))
pres <- names(cf)
pos.pres <- pmatch(pres,distr)
mis <- distr[-pmatch(pres,distr)]
pos.mis <- pmatch(mis,distr)
control.family <- replicate(length(distr),list())
control.family[[pos.pres]] <- cf
control.family[[mis.pres]] <- cf[[mis.pres]]



### CHECKS FLEXIBLE SPLINE MODEL IN INLA & MLE
sp1 <- flexsurvspline(Surv(recyrs,censrec)~group,data=bc,k=1,scale="hazard")
sp2 <- inla(inla.surv(recyrs,censrec)~group+f(log(bc$recyrs),model="rw2",hyper = list(prec = list(prior="loggamma",param=c(1,0.01)))),data=bc,family="weibull")
 


#### FOR KATRIN
formula <- Surv(time,censored)~arm
weibull <- fit.models(formula=formula,data=dat,distr="weibull",method="mcmc",n.iter=0)




setwd("~/Dropbox/HE/CERM/Roche/Work/RModel/Markov/DigitiseIT/Arm covariate/")
library(flexsurv)
library(INLA)
library(R2OpenBUGS)
source("~/Dropbox/UCL/Mapi/Projects/Survival/survHE.R")
data <- make.ipd(list("IPD_PFS_Benda.txt","IPD_PFS_GBenda.txt"),ctr=1,var.labs=c("Patient","time","event","arm"))
data <- data[,-5]
data$arm[data$arm==1] <- 2; data$arm[data$arm==0] <- 1; data$arm[data$arm==2] <- 0
data.inla <- data; data.inla$time <- data.inla$time/10

formula <- Surv(time,event)~as.factor(arm)
fit.mle <- fit.models(formula=formula,data=data,distr=c("weibull","exp"))
fit.mcmc <- fit.models(formula=formula,data=data,distr="weibull",method="mcmc",n.iter=6000)
fit.inla <- fit.models(formula=formula,data=data.inla,distr="weibull",method="inla")

data$time <- data$time/10
