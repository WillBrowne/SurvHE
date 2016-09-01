## SET OF UTILITY FUNCTIONS TO INCLUDE SURVIVAL ANALYSIS RESULTS INTO A HEALTH ECONOMIC MODEL
## Gianluca Baio, 4 Mar 2016

fit.models <- function(formula=NULL,data,distr=NULL,method="mle",...) {
  ## Main function - runs the survival analysis with several useful options
  ## formula = a formula specifying the model to be used, in the form
  ##           Surv(time,event)~treatment[+covariates] for flexsurv
  ##           inla.surv(time,event)~treatment[+covariates] for INLA
  ## data = a data frame containing the data to be used
  ## distr = a (vector of) string(s) containing the name(s) of the model(s) to be fitted
  ## method = a string specifying the inferential method ("mle", "inla" or "mcmc) 
  ## 
  ## ... additional options (mainly to do with INLA & MCMC)
  ##
  ## **INLA** specific options 
  ## dz = defines the step length for the grid search over the hyperparameters space (default = 0.1)
  ## diff.logdens = defines the difference in the log-density for the hyperparameters to stop integration (default = 5)
  ## control.fixed = defines the default for the priors, unless specified by the user. Default values are
  ##                 prior mean = 0 for *all* fixed effects
  ##                 prior var = 1000 for *all* fixed effects
  ##                 prior mean = 0 for the intercept
  ##                 prior prec -> 0 for the intercept 
  ## control.family = a list of options. If distr is a vector, then can be provided as a named
  ##                  list of options, for example something like this: 
  ##                  control.family=list(weibull=list(param=c(.1,.1)),lognormal=list(initial=2))
  ##                  the names of the elements of the list need to be the same as those given
  ##                  in the vector distr
  ##
  ## **MCMC** specific options
  ## n.iter = total number of iterations (default = 10000). If n.iter=0, then only
  ##          creates the data, parameters and model file to then run BUGS separately
  ## n.samples = number of iterations to be retained for summary (default = 1000)
  ## n.chains = number of chains (default = 2)
  ## n.burnin = number of iterations to be discarded (default = n.iter-n.samples/n.chains)
  ## n.thin = thinning (default = 1)
  ## debug.mode = whether the OpenBUGS window should open (default = TRUE for Windows and FALSE for Linux/Mac OS)
  ## OpenBUGS.pgm = full path to the OpenBUGS launcher
  ## useWINE = whether should use WINE (under Linux or Mac OS) - default = FALSE
  ## newWINE = whether should use newer version of WINE (under Linux or Mac OS) - default = TRUE
  ## WINEPATH = - default = NULL
  ## WINE = - default = NULL
  ## max_splines 
  
  ##### NB: COULD SPECIFY DIFFERENT control.family OPTIONS FOR THE INLA METHOD WITH MANY DISTRIBUTIONS
  
  # Lists all the additional inputs
  exArgs <- list(...)
  # Avoids the 'no visible binding for global variable' error, when compiling
  model <- NULL
  
  # Needs to specify either the formula or the list of variables!
  if(is.null(formula)) {
    stop("You need to specify a model 'formula', e.g. 'formula=Surv(time,event)~treat'")
  }
  # ensures method is lower case
  method <- tolower(method)
  # ensire method is one of "mle","inla" or "mcmc"
  if(!method %in% c("mcmc","inla","mle","splines")) {
    stop("Methods available for use are mcmc, inla, mle or splines")
  }
 if(is.null(distr) & method != "splines") {
    stop("You need to specify a distribution")
  }
  # INLA can only do a limited set of models (for now) so if user has selected
  # one that is not available, then falls back on MLE analysis
  availables.mle <- c("genf", "genf.orig", "gengamma", "gengamma.orig", "exp", 
                      "weibull", "weibullPH", "lnorm", "gamma", "gompertz", 
                      "llogis", "exponential", "lognormal")
  availables.inla <- c("exponential","weibull","lognormal","loglogistic")
  availables.mcmc <- c("weibull","exponential","gengamma","lognormal","gamma","loglogistic","genf","gompertz","dexp","weibullPH","dloglogis")

  
  # Standardises labels for model names
  labs <- distr
  labs[pmatch("weibull",labs)] <- "Weibull"
  labs[pmatch("weibullPH",labs)] <- "Weibull"
  labs[pmatch("exp",labs)] <- "Exponential"
  labs[pmatch("exponential",labs)] <- "Exponential"
  labs[pmatch("gamma",labs)] <- "Gamma"
  labs[pmatch("lnorm",labs)] <- "log-Normal"
  labs[pmatch("lognormal",labs)] <- "log-Normal"
  labs[pmatch("llogis",labs)] <- "log-Logistic"
  labs[pmatch("loglogistic",labs)] <- "log-Logistic"
  labs[pmatch("gengamma",labs)] <- "Gen. Gamma"
  labs[pmatch("genf",labs)] <- "Gen. F"
  labs[pmatch("gompz",labs)] <- "Gompertz"
  labs[pmatch("dexp",labs)] <- "Bastardised exp"
  labs[pmatch("dloglogis",labs)] <- "log-Logistic"

  if(method=="inla") {
    # Checks that the distribution name(s) are consistent with INLA
    user.distr <- distr
    distr[pmatch("llogis",user.distr)] <- "loglogistic"
    distr[pmatch("exp",user.distr)] <- "exponential"
    distr[pmatch("lnorm",user.distr)] <- "lognormal"
    # But if there still are some that are just not available then falls back on MLE
    if (any(is.na(pmatch(distr,availables.inla)))) {
      method <- "mle"
      cat("NB: INLA can only fit Exponential, Weibull, log-Logistic or log-Normal parametric survival models. \nFalling back on MLE analysis")
     
    }
  }
  if(method=="mcmc") {
      # Checks that the distribution name(s) are consistent with MCMC
      user.distr <- distr
      distr[pmatch("llogis",user.distr)] <- "loglogistic"
      distr[pmatch("exp",user.distr)] <- "exponential"
      distr[pmatch("lnorm",user.distr)] <- "lognormal"
  }
  # Reconstructs the vars list based on the formula
  test <- attributes(terms(formula))$term.labels
  ncovs <- length(test)
  formula.temp <- as.formula(gsub("inla.surv","Surv",deparse(formula)))
  time <- all.vars(formula.temp,data)[1]
  event <- all.vars(formula.temp,data)[2]
  if (ncovs>0) {
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
  } else {
      covs <- factors <- NULL
  }
  # If there are covariates, creates a matrix and sets the dimension
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
  vars <- list(time=time,event=event,factors=factors,covs=covs,nlevs=D)
  
  # Need to create a formula for the KM that is in the right format (flexsurv-like)
  chk <- is.na(pmatch("Surv",attributes(terms(formula))$variables[[2]][1])) # If TRUE, then it's in inla.surv() terms
                                                                            # If FALSE, then formula is in Surv() terms
  
  # If it's in inla.surv() needs to create a Surv() formula for the KM + check that the method isn't mle
  if(chk) {
    tmp <- deparse(formula)
    km.formula <- as.formula(gsub("inla.surv","Surv",tmp))
    # If the method was originally INLA but has been modified or if the method is mle but the formula is in the wrong
    # format, then need to change the formula to be consistent with the Surv() notation
    if (method=="mle") {
      formula <- km.formula
    }
  } else {
    km.formula <- formula
    # If the method isn't inla but the formula is in Surv() terms, then change the formula
    if (method=="inla") {
      tmp <- deparse(formula)
      formula <- as.formula(gsub("Surv","inla.surv",tmp))
    }
  }
  # Computes the Kaplan Meier curve using the package "rms"
  ObjSurvfit=rms::npsurv(        # Uses the function "npsurv" from the package "rms"
    formula=km.formula,          # to fit the model specified in the "formula" object
    data=data                    # to the dataset named "data"
  )

  # If method = MLE, then fits the model(s) using flexsurvreg
  if (method=="mle") {
    # Checks that the distribution name(s) are consistent with flexsurv
    # The only problem here is if the user has specified a log-Logistic in INLA terminology
    user.distr <- distr
    distr[pmatch("loglogistic",user.distr)] <- "llogis"
    # Then run the model(s)
    runMLE <- function(distr) {
      tic <- proc.time()
      model <- flexsurv::flexsurvreg(formula=formula,data=data,dist=distr)
      toc <- proc.time()-tic
      time2run=toc[3]
      list(model=model,time2run=time2run)
    }
    output <- lapply(distr,function(x) runMLE(x))
    mod <- lapply(output, function(x) x$model)
    time2run <- unlist(lapply(output, function(x) x$time2run)); names(time2run) <- labs
    aic <- unlist(lapply(mod,function(x) x$AIC))
    bic <- unlist(lapply(mod,function(x) -2*x$loglik+length(x$coefficients)*log(x$N)))
    ## This could be calculated?
    dic <- NULL
  }
  
  # If method = INLA, then fits model(s) using inla
  if (method=="inla") {
    # If INLA is not installed, then asks for it
    if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
     stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
    }
    
    # Set up optional parameters to default values if the user hasn't done it themselves
    # 1. defines the step length for the grid search over the hyperparameters space
    if(exists("dz",where=exArgs)) {dz <- exArgs$dz} else {dz <- 0.1}
    # 2. defines the difference in the log-density for the hyperparameters to stop integration
    if(exists("diff.logdens",where=exArgs)) {diff.logdens <- exArgs$diff.logdens} else {diff.logdens <- 5}
    # 3. defines the default for the priors, unless specified by the user
    if(exists("control.fixed",where=exArgs)) {
      control.fixed <- exArgs$control.fixed
    } else {
      control.fixed <- INLA::inla.set.control.fixed.default()
      # prior mean = 0 for *all* fixed effects
      # prior var = 1000 for *all* fixed effects
      # prior mean = 0 for the intercept
      # prior prec -> 0 for the intercept 
    }
    if(exists("control.family",where=exArgs)) {
      control.family <- replicate(length(distr),list(INLA::inla.set.control.family.default()))
      names(control.family) <- distr
      cf <- exArgs$control.family
      news <- pmatch(names(cf),names(control.family))
      for (i in 1:length(news)) {
        control.family[[news[i]]] <- cf[[i]]
        names(control.family[[news[i]]]) <- names(cf[[i]])
      }
    } else {
      # If not specified, uses the default values in INLA (depending on the model selected)
      control.family <- replicate(length(distr),list(INLA::inla.set.control.family.default()))
      # And sets some more sensible/general values for some of the controls
#       if(!is.na(pmatch("weibull",distr))) {
#         control.family[[pmatch("weibull",distr)]]$param <- c(.1,.1)
#       }
#        for (i in 1:length(control.family)) {
#          control.family[[i]]$initial=.1
#        }
    }

    # 4. Finally runs INLA
    mod <- lapply(1:length(distr), function(x) {
      INLA::inla(formula,family=distr[x],data=data,control.compute=list(config=TRUE,dic=TRUE),
                 control.inla=list(int.strategy="grid",dz=dz,diff.logdens=diff.logdens),
                 control.fixed=control.fixed,control.family=control.family[[x]]
      )
    })
    time2run <- unlist(lapply(mod,function(x) x$cpu.used["Total"])); names(time2run) <- labs
    # NB Internally, DIC = model$dic$mean.deviance+model$dic$p.eff
    dic <- unlist(lapply(mod,function(x) x$dic$dic))
    # NB But to estimate the AIC & BIC is probably best to use model$dic$deviance.mean!
    aic <- unlist(lapply(1:length(mod), function(i) 2*mod[[i]]$dic$p.eff+mod[[i]]$dic$deviance.mean)); names(aic) <- NULL
    bic <- unlist(lapply(1:length(mod),function(i) 
      mod[[i]]$dic$deviance.mean+mod[[i]]$dic$p.eff*log(mod[[i]]$size.linear.predictor$n))); names(bic) <- NULL
    for (i in 1:length(distr)) {mod[[i]]$dlist$name <- distr[i]}
  }
  
  # If method="MCMC" then fits model(s) using BUGS
  if(method=="mcmc") {
    if(!isTRUE(requireNamespace("R2OpenBUGS",quietly=TRUE))) {
     stop("You need to install the package 'R2OpenBUGS'. Please run in your R terminal:\n install.packages('R2OpenBUGS')")
    }
    print("writing model")
    # 1. Creates a little function that writes the relevant model file
    write.model <- function(position){
      # Associates models with indices (just to read the code more neatly)
      weib <- 1; expo <- 2; ggam <- 3; lnor <- 4; gamm <- 5; llog <- 6 ; genf <- 7 ; gompz <- 8 ; dexp <- 9 ; weibPH <-10 ; dloglogis <- 11
      
      # Start model
      start.mod <- "model {"
      # Start loop
      start.loop <- "for (i in 1:n) {"
      # Data model (OpenBUGS notation)
      mod.data <- character()
      mod.data[weib] <- "t[i] ~ dggamma(1,lambda[i],shape)C(cens[i],)"              # weibull AF
      mod.data[expo] <- "t[i] ~ dexp(lambda[i])C(cens[i],)"                         # exponential
      mod.data[ggam] <- "t[i] ~ dggamma(q,lambda[i],shape)C(cens[i],)"              # generalised gamma
      mod.data[lnor] <- "t[i] ~ dlnorm(lambda[i],prec)C(cens[i],)"                  # log-normal
      mod.data[gamm] <- "t[i] ~ dgamma(shape,lambda[i])C(cens[i],)"                 # gamma
      mod.data[llog] <- "t.log[i] ~ dlogis(lambda[i],tau)C(cens.log[i],)"           # log-logistic
      mod.data[genf] <- "t[i] ~ df(df1,df2,lambda[i],shape)C(cens[i],)"        # generalised F
      mod.data[gompz] <- "dummy[i] <- 0\n dummy[i] ~ dloglik(logLike[i]) \n logLike[i] <- -lambda[i]/shape * (exp(shape * t[i]) -1) + death[i]*(log(lambda[i])+shape*t[i])"        # Gompertz dist doing the loglik trick 
      mod.data[dexp] <- "dummy[i] <- 0 \n dummy[i] ~ dloglik(logLike[i]) \n logLike[i] <- death[i]*log(lambda[i])-lambda[i]*t[i]"
      mod.data[weibPH] <- "t[i] ~ dweib(shape,lambda[i])C(cens[i],)" # weibullPH
      mod.data[dloglogis] <- "dummy[i] <- 0 \n dummy[i] ~ dloglik(logLike[i]) \n logLike[i] <- death[i]*(log(shape) - log(lambda[i]) + (shape - 1)*(log(t[i]) - log(lambda[i])) - log(1 + (pow(t[i]/lambda[i],shape)))) -log(1 + (pow(t[i]/lambda[i],shape))) " # Loglogistic

      # ... all other models I want to implement
      # Linear predictor 
      linpred <- "beta"
      if(!is.null(covs)) {
        if (K==1) {
          linpred <- paste0(linpred,"+",paste0("gamma*X[i]",collapse="+"))
        } else {
          linpred <- paste0(linpred,"+",paste0(paste0(paste0("gamma[",1:K,"]"),"*"),
                            "X[i,",1:K,"]",collapse="+"))
          print(linpred)
        }
      }
      if(!is.null(factors)) {
        if (length(factors)==1) {
          linpred <- paste0(linpred,"+",paste0("delta[",factors,"[i]]"),collapse="+")
        } else {
          linpred <- paste0(linpred,"+",paste0(unlist(lapply(1:H,function(i) paste0("delta",i,"[",factors[i],"[i]]"))),collapse="+"))
          print(linpred)
        }
      }
      ##################################################################################
      ## Would be needed to remove the last '+', but got rid of this problem already  ##
      ## linpred <- gsub(".{1}$","",paste0("log(lambda[i]) <- ",linpred))             ##
      ##################################################################################
      
      # For the lognormal and the loglogistic don't need to log the linear predictor! 
      if (all(position!=c(4,6))) {
        # So if position is not 4 or 6, needs to use log(lambda[i])
        linpred <- paste0("log(lambda[i]) <- ",linpred)
      } else {
        # But if it is  either 4 or 6 then don't need to log!
        linpred <- paste0("lambda[i] <- ",linpred)
      }
      # End loop
      end.loop <- "}"
      # Priors - regression coefficients
      intercept <- "beta ~ dnorm(mu.beta,tau.beta)\n"
      if(is.null(covs)) {
        cov.eff <- ""
      }
      if(!is.null(covs)) {
        if (K==1) {
          cov.eff <- "gamma ~ dnorm(mu.gamma,tau.gamma)\n"
        } else {
          cov.eff <- "gamma[1:K] ~ dmnorm(mu.gamma[],Q.gamma[,])\n"
        }
      }
      if(is.null(factors)) {
        fact <- ""
      }
      if(!is.null(factors)) {
        if (H==1) {
          fact <- "delta[1] <- 0 \nfor (d in 2:D) {\ndelta[d]~dnorm(mu.delta,tau.delta)\n}\n"
        } else {
          fact <- paste0(unlist(lapply(1:H,function(i) paste0("delta",i,"[1] <-0 \nfor(d in 2:D[",i,"]){\ndelta",i,"[d]~dnorm(mu.delta,tau.delta)\n}\n"))),collapse="")
        }
      }
      # Priors - distribution-specific parameters
      parameters <- character()
      parameters[weib] <- "shape ~ dgamma(a.shape,b.shape)\n"
      parameters[expo] <- " "
      parameters[ggam] <- "shape ~ dgamma(a.shape,b.shape) \nq ~ dgamma(a.q,b.q)\n"
      parameters[lnor] <- "sd ~ dunif(0,10) \nprec <- pow(sd,-2)\n"
      parameters[gamm] <- parameters[1]
      parameters[llog] <- "tau ~ dgamma(a.shape,b.shape)\n"
      parameters[genf] <- "shape ~ dgamma(a.shape,b.shape) \n df1 ~ dgamma(a.n,b.n) \n df2 ~ dgamma(a.m,b.m)\n "
      parameters[gompz] <- "shape ~ dgamma(a.shape,b.shape)\n "
      parameters[dexp] <- " "
      parameters[weibPH] <- "shape ~ dgamma(a.shape,b.shape)\n "
      parameters[dloglogis] <- "shape ~ dgamma(a.shape,b.shape)\n "

      # ... priors for all the other models I will implement

      # End model
      end.mod <- "}"
      
      # Human-readable labels to specify the distributional assumption for the observed data
      labels <- character()
      labels[weib] <- paste0("# Weibull model (AF parameterisation) - written on: ",Sys.time())
      labels[expo] <- paste0("# Exponential model - written on: ",Sys.time())
      labels[ggam] <- paste0("# Generalised Gamma model - written on: ",Sys.time())
      labels[lnor] <- paste0("# log-Normal model - written on: ",Sys.time())
      labels[gamm] <- paste0("# Gamma model - written on: ",Sys.time())
      labels[llog] <- paste0("# log-Logistic model - written on: ",Sys.time())
      labels[genf] <- paste0("# Generalised F model - written on: ",Sys.time())
      labels[gompz] <- paste0("# Gompertz model - written on: ",Sys.time())
      labels[dexp] <- paste0("# Exponential dummy model - written on: ",Sys.time())
      labels[weibPH] <- paste0("# Weibull model (PH parameterisation) - written on: ",Sys.time())
      labels[dloglogis] <- paste0("# Loglogistic model  - written on: ",Sys.time())
      
      # Determines which data model has been used and selects the relevant text
      model.string <- paste0(labels[position],"\n",start.mod,"\n",start.loop,"\n",
                           mod.data[position],"\n",linpred,"\n",end.loop,"\n",
                           intercept,cov.eff,fact,parameters[position],end.mod)
      
      filein <- paste0(getwd(),"/model.txt")
      writeLines(model.string, con=filein)
      list(name=filein,string=model.string)
    }
    
        # 2. Creates a little function that writes the relevant data to a list (to be passed to OpenBUGS)
    write.data <- function(position) {
      # Associates models with indices (just to read the code more neatly)
      weib <- 1; expo <- 2; ggam <- 3; lnor <- 4; gamm <- 5; llog <- 6; genf <- 7; gompz <-8; dexp <- 9; weibPH <- 10 ; dloglogis <- 11

      # Basic data (observed variables -- irrespective of the distribution)
      # Needs to define defaults for the parameters of the intercept & trt.effect
      mu.beta=mu.delta=0
      tau.beta=tau.delta=.0001
      ### CHECK: DO I NEED TO USE eval(parse(text=time))???
      if(position==6) {
        dataBugs <- list(
          "cens.log"=ifelse(eval(parse(text=paste0("data$",event)))==0,eval(parse(text=paste0("log(data$",time,")")))+1,0),
          "t.log"=ifelse(eval(parse(text=paste0("data$",event)))==1,eval(parse(text=paste0("log(data$",time,")")))+1,NA),
          mu.beta=mu.beta,tau.beta=tau.beta,n=dim(data)[1]
        )
#         dataBugs <- with(data,
#                          list("cens.log"=ifelse(eval(parse(text=event))==0,log(time),-20),
#                               "t.log"=ifelse(eval(parse(text=event))==1,log(time),NA),
#                               mu.beta=mu.beta,tau.beta=tau.beta,n=dim(data)[1]
#                          )
#         )
      } else if (position %in% c(7,8,9,11)){
        dataBugs <- list(
          "death"= ifelse(eval(parse(text=paste0("data$",event)))==1,1,0),
          "t"=eval(parse(text=paste0("data$",time))),
          mu.beta=mu.beta,tau.beta=tau.beta,n=dim(data)[1]
          )
      } else {
        dataBugs <- list(
          "cens"=ifelse(eval(parse(text=paste0("data$",event)))==0,eval(parse(text=paste0("data$",time))),0),
          "t"=ifelse(eval(parse(text=paste0("data$",event)))==1,eval(parse(text=paste0("data$",time))),NA),
          mu.beta=mu.beta,tau.beta=tau.beta,n=dim(data)[1]
        )
#         dataBugs <- with(data,
#                          list("cens"=ifelse(eval(parse(text=event))==0,time,0),
#                               "t"=ifelse(eval(parse(text=event))==1,time,NA),
#                               mu.beta=mu.beta,tau.beta=tau.beta,n=dim(data)[1]
#                          )
#         )
      }
      
      # Now need to add the other data
      if(!is.null(covs)) {
        if (K==1) {
          mu.gamma <- 0
          tau.gamma <- .0001
          regr2 <- "dataBugs$mu.gamma=mu.gamma;dataBugs$tau.gamma=tau.gamma;dataBugs$X=X"
        } else {
          mu.gamma <- rep(0,K)
          Q.gamma <- .00001*diag(K)
          regr2 <- "dataBugs$mu.gamma=mu.gamma;dataBugs$Q.gamma=Q.gamma;dataBugs$K=K;dataBugs$X=as.matrix(X)"
        }
        eval(parse(text=regr2))
      }
      if(!is.null(factors)) {
        n.temp <- length(dataBugs)
        for (i in 1:H) {
          dataBugs[[(n.temp+i)]] <- as.numeric(data[,cols[i]])
          names(dataBugs)[(n.temp+i)] <- factors[i]
        }
        regr3 <- "dataBugs$mu.delta=mu.delta;dataBugs$tau.delta=tau.delta;dataBugs$D=D"
        eval(parse(text=regr3))
      }
      hyperpars <- character()
      hyperpars[weib] <- "dataBugs$a.shape=0.01;dataBugs$b.shape=0.01" 
      hyperpars[expo] <- ""
      hyperpars[ggam] <- "dataBugs$a.shape=0.01;dataBugs$b.shape=0.01;dataBugs$a.q=0.01;dataBugs$b.q=0.01"
      hyperpars[lnor] <- ""
      hyperpars[gamm] <- hyperpars[1]
      hyperpars[llog] <- hyperpars[1] 
      hyperpars[genf] <-  "dataBugs$a.shape=0;dataBugs$b.shape=1000;"
      hyperpars[gompz] <-  "dataBugs$a.shape=0.001;dataBugs$b.shape=0.001;"
      hyperpars[dexp] <- ""
      hyperpars[weibPH] <- hyperpars[1] 
      hyperpars[dloglogis] <-  "dataBugs$a.shape=0.001;dataBugs$b.shape=0.001;"

      eval(parse(text=hyperpars[position]))
      return(dataBugs)
    }
    # 3. Creates a little function that writes the inits (to be passed to OpenBUGS)
    # Inits
  inits <- function(){
      # Associates models with indices (just to read the code more neatly)
      weib <- 1; expo <- 2; ggam <- 3; lnor <- 4; gamm <- 5; llog <- 6 ; genf <- 7; gompz <- 8 ; dexp <- 9 ; weibPH <- 10 ; dloglogis <- 11
      position <- pmatch(distr,availables.mcmc)
      dataBugs <- write.data(position)
      
      inits.list <- list()
      if (position==6) {
        inits.t <- eval(parse(text="inits.list$t.log=ifelse(dataBugs$cens.log!=-20,log(exp(dataBugs$cens.log)+1),NA)"))
      } else if (position %in% c(7,8,9,11)){
        #inits.t <- eval(parse(text="inits.list$t=dataBugs$t+1"))
        print("t is observed")

      }else {
        inits.t <- eval(parse(text="inits.list$t=ifelse(dataBugs$cens!=0,(dataBugs$cens)+1,NA)"))
      }
      inits.par <- character()
      inits.par[weib] <- paste0("inits.list$shape=runif(1)")
      inits.par[expo] <- ""
      inits.par[ggam] <- paste0("inits.list$shape=runif(1); inits.list$q=runif(1)")
      inits.par[lnor] <- paste0("inits.list$sd=runif(1)")
      inits.par[gamm] <- paste0("inits.list$shape=runif(1)")
      inits.par[llog] <- paste0("inits.list$tau=runif(1)")
      inits.par[genf] <- paste0("inits.list$shape=runif(1);inits.list$df1=runif(1);inits.list$df2=runif(1)")
      inits.par[gompz] <- paste0("inits.list$shape=runif(1)")
      inits.par[dexp] <- ""
      inits.par[weibPH] <- paste0("inits.list$shape=runif(1)")
      inits.par[dloglogis] <- paste0("inits.list$shape=runif(1)")


      eval(parse(text=inits.par[position]))
      inits.coef <- paste0("inits.list$beta=rnorm(1)")
      eval(parse(text=inits.coef))
      if(!is.null(covs)) {
        if(K==1){
          inits.coef3 <- "inits.list$gamma=rnorm(1)"
        } else {
          inits.coef3 <- "inits.list$gamma=rnorm(K,0,1)"
        }
        eval(parse(text=inits.coef3))
      }
      if(!is.null(factors)) {
        if(H==1) {
          inits.coef4 <- "inits.list$delta=c(NA,rnorm((D-1),0,1))"
        } else {
          inits.coef4 <- paste0("inits.list$delta",1:H,"=c(NA,rnorm((D[",1:H,"]-1),0,1))",collapse=";")
        }
        eval(parse(text=inits.coef4))
      }
      return(inits.list)
    }


    
    # 4. Checks that optional parameters to be passed to OpenBUGS are not given as extra arguments and if not
    #    sets the default values
    if(exists("n.chains",where=exArgs)) {n.chains <- exArgs$n.chains} else {n.chains <- 2}
    if(exists("n.iter",where=exArgs)) {n.iter <- exArgs$n.iter} else {n.iter <- 10000}
    if(exists("n.samples",where=exArgs)) {n.samples <- exArgs$n.samples} else {n.samples <- 1000}
    if(exists("n.burnin",where=exArgs)) {n.burnin <- exArgs$n.burnin} else {n.burnin <- n.iter-n.samples/n.chains}
    if(exists("n.thin",where=exArgs)) {n.thin <- exArgs$n.thin} else {n.thin <- 1}
    if(exists("debug.mode",where=exArgs)) {debug.mode <- exArgs$debug.mode} else {
      # For Windows machines, makes debug=T the default. If not Windows, then debug=F by default
      if(Sys.info()["sysname"] == "Windows") {debug.mode <- TRUE} else {debug.mode <- FALSE}
    }
    if(exists("OpenBUGS.pgm",where=exArgs)) {OpenBUGS.pgm <- exArgs$OpenBUGS.pgm} else {OpenBUGS.pgm=NULL}
    if(exists("useWINE",where=exArgs)) {useWINE <- exArgs$useWINE} else {useWINE=FALSE}
    if(exists("newWINE",where=exArgs)) {newWINE <- exArgs$newWINE} else {newWINE=TRUE}
    if(exists("WINEPATH",where=exArgs)) {WINEPATH <- exArgs$WINEPATH} else {WINEPATH=NULL}
    if(exists("WINE",where=exArgs)) {WINE <- exArgs$WINE} else {WINE=NULL}
    
    # 5. Creates a little function that defines the parameters to be monitored
    write.params <- function(position) {
      # Associates models with indices (just to read the code more neatly)
      weib <- 1; expo <- 2; ggam <- 3; lnor <- 4; gamm <- 5; llog <- 6 ; genf <-7 ; gompz <- 8 ; dexp <- 9 ; weibPH <- 10; ; dloglogis <- 11
      
      params <- character()
      params[weib] <- paste0("c('beta','shape'")
      params[expo] <- paste0("c('beta'")
      params[ggam] <- paste0("c('beta'")
      params[lnor] <- paste0("c('beta','sd'")
      params[gamm] <- paste0("c('beta','shape'")
      params[llog] <- paste0("c('beta','tau',")
      params[genf] <- paste0("c('beta','shape','df1','df2'")
      params[gompz] <- paste0("c('beta','shape'")
      params[dexp] <- paste0("c('beta'")
      params[weibPH] <- paste0("c('beta','shape'")
      params[dloglogis] <- paste0("c('beta','shape'")

      if (!is.null(covs)) {
        params <- paste(params,"'gamma'",sep=",")
      }
      if (!is.null(factors)) {
        if (H==1) {
          params <- paste(params,"'delta'",sep=",")
        } else {
          params <- paste(params,paste0("'delta",1:H,"'",collapse=","),sep=",")
        }
      }
      params <- paste0(params,")")
      eval(parse(text=params[position]))
    }
    # If n.iter is set to 0, then will only create model file, data list and inits list
    if(n.iter==0) {
      mod <- list()
      position <- pmatch(distr,availables.mcmc)
      dataBugs <- mod$dataBugs <- write.data(position) 
      mod$param <- write.params(position)
      mod$model.file <- write.model(position)
      mod$inits <- lapply(1:n.chains,function(i) inits())
      aic <- bic <- dic <- time2run <- NULL
    } else {
      runBUGS <- function(distr) {        
        print("get position")
        position <- pmatch(distr,availables.mcmc)
        print(position)
        print(distr)
        print("write data")
        dataBugs <- write.data(position)
        print("get params")
        params <- write.params(position)
        print("write model")
        model.file <- write.model(position)$name
        # Determines which model assumption is selected
        tic <- proc.time()
        model <- R2OpenBUGS::bugs(data=dataBugs,inits=inits,parameters.to.save=params,
                                  model.file=model.file,n.chains=n.chains,n.iter=n.iter,
                                  n.burnin=n.burnin,n.thin=n.thin,debug=debug.mode,
                                  OpenBUGS.pgm=OpenBUGS.pgm,useWINE=useWINE,newWINE=newWINE,
                                  WINEPATH=WINEPATH,WINE=WINE)
        toc <- proc.time()-tic
        time2run <- toc[3]
        list(model=model,time2run=time2run)
      }
      output <- lapply(distr,function(x) runBUGS(x))
      mod <- lapply(output, function(x) x$model)
      time2run <- unlist(lapply(output, function(x) x$time2run)); names(time2run) <- labs
      for (i in 1:length(distr)) {mod[[i]]$dlist$name <- distr[i]}
     for (i in 1:length(distr)){
      # NEED TO RECOMPUTE THE DEVIANCE FOR THE log-logistic MODEL WHICH IS FITTED TO *log* DATA!
      if(mod[[i]]$dlist$name=="loglogistic") {
        print("recomputing log")
        position <- pmatch(mod[[i]]$dlist$name,availables.mcmc)
        dataBugs <- mod[[i]]$dataBugs <- write.data(position) 
        # Computes the density for the log-logistic model
        f <- function(t,shape,scale){
          num <- (shape/scale)*(t/scale)^(shape-1)
          den <- (1+(t/scale)^shape)^2
          f <- num/den
          f <- ifelse(is.na(f),1,f)
          return(f)
        }
        # Computes the survival function based on the log-logistic model
        S <- function(t,shape,scale){
          num <- 1
          den <- 1+(t/scale)^shape
          S <- num/den
          S <- ifelse(is.na(S),1,S)
          return(S)
        }
        
        # Defines shape & scale based on the estimated parameters
        shape=mod[[i]]$sims.list$tau
        scale=exp(mod[[i]]$sims.list$lambda)
        scale1 <- scale[,which(!is.na(dataBugs$t.log))]
        scale2 <- scale[,which(is.na(dataBugs$t.log))]
        
        # Defines the relevant times (on the NATURAL scale)
        obs.time <- exp(dataBugs$t.log[!is.na(dataBugs$t.log)])
        cens.time <- exp(dataBugs$cens.log[is.na(dataBugs$t.log)])
#       all.time <- ifelse(is.na(dataBugs$t.log),exp(dataBugs$cens.log),exp(dataBugs$t.log))
        # Computes the log-likelihood (including censoring) & the deviance on the NATURAL scale
        log.lik1 <- matrix(0,mod[[i]]$n.sims,length(obs.time))
        log.lik2 <- matrix(0,mod[[i]]$n.sims,length(cens.time))
        # log.lik3 <- matrix(0,mod[[1]]$n.sims,length(all.time))
        

        log.lik1 <- matrix(unlist(lapply(1:dim(log.lik1)[1],function(i) {
          lapply(1:length(obs.time),function(j) {
            sum(log(f(obs.time[j],shape[i],scale1[i,j])))
          })
        })),nrow=mod[[i]]$n.sims,ncol=length(obs.time))
        
        log.lik2 <- matrix(unlist(lapply(1:dim(log.lik2)[1],function(i) {
          lapply(1:length(cens.time),function(j) {
            sum(log(S(cens.time[j],shape[i],scale2[i,j])))
          })
        })),nrow=mod[[i]]$n.sims,ncol=length(cens.time))
        
#         log.lik3 <- matrix(unlist(lapply(1:dim(log.lik3)[1],function(i) {
#           lapply(1:length(all.time),function(j) {
#             sum(log(S(all.time[j],shape[i],scale[i])))
#           })
#         })),nrow=model$n.sims,ncol=length(all.time))
        
        log.lik <- apply(log.lik1,1,sum)+apply(log.lik2,1,sum)###-apply(log.lik3,1,sum)
        deviance <- -2*log.lik
        mod[[i]]$sims.list$deviance.log <- mod[[i]]$sims.list$deviance
        mod[[i]]$sims.list$deviance <- deviance
        mod[[i]]$summary["deviance",1:7] <- c(mean(deviance),sd(deviance),quantile(deviance,.025),quantile(deviance,.25),
                                           median(deviance),quantile(deviance,.75),quantile(deviance,.975))
        mod[[i]]$DIC <- mod[[i]]$summary["deviance",1]+mod[[i]]$pD
      }
      }
      aic <- unlist(lapply(mod, function(x) 2*x$pD+x$summary["deviance",1]))
      bic <- unlist(lapply(mod, function(x) x$summary["deviance",1]+x$pD*log(sum(ObjSurvfit$n))))
      dic <- unlist(lapply(mod, function(x) x$DIC))
    }
  }
  if (method == "splines"){
     if(!exists("knots",where=exArgs)) {
    knots <- 5
  } else { knots <- exArgs$knots}
  if (knots > 10){
    stop("More than ten knots may overfit the data, please choose fewer knots")
  }
spline_parameter_tuning <- function(knot,scale,data){
       model <- flexsurvspline(formula, data=data, k=knot, scale=scale)
       return(model$AIC)
      }
spline_dataframe <- expand.grid(seq(1,knots), c("hazard","odds","normal"))
spline_dataframe$AIC <- apply(spline_dataframe,1,function(x) spline_parameter_tuning(as.numeric(x[1]),x[2],data))
      ### Idea for testing. Take all non-censored data in the 
print(spline_dataframe)
k <- spline_dataframe[spline_dataframe$AIC == min(spline_dataframe$AIC),]$Var1
scale <- as.character(spline_dataframe[spline_dataframe$AIC == min(spline_dataframe$AIC),]$Var2)
tic <- proc.time()
### add knots and scale used in final model
mod <- list(flexsurvspline(formula, data=data, k=k, scale=scale))
toc <- proc.time()-tic

time2run=toc[3]
aic <- mod[[1]]$AIC
bic <- NULL # Need to try and make this work
dic <- NULL # Need to try and make this work

    }
  # Now defines the outputs of the function
  model.fitting <- list(aic=aic,bic=bic,dic=dic)
  misc <- list(time2run=time2run,formula=formula,km=ObjSurvfit,data=data)
  if(method=="mcmc") {misc$vars <- vars}
  res <- list(models=mod,model.fitting=model.fitting,method=method,misc=misc)
  class(res) <- "survHE"
  return(res)
}


make.surv <- function(fit,mod=1,t=NULL,newdata=NULL,nsim=1,...) {
  ## Creates the survival curves for the fitted model(s)
  # fit = the result of the call to the fit.models function, containing the model fitting (and other relevant information)
  # mod = the index of the model. Default value is 1, but the user can choose which model fit to visualise, 
  #     if the call to fit.models has a vector argument for distr (so many models are fitted & stored in the same object)
  # t = the time framework to be used for the estimation of the survival curve
  # newdata = a list (of lists), specifiying the values of the covariates at which the computation is performed. For example
  #           'list(list(arm=0),list(arm=1))' will create two survival curves, one obtained by setting the covariate 'arm'
  #           to the value 0 and the other by setting it to the value 1. In line with 'flexsurv' notation, the user needs
  #           to either specify the value for *all* the covariates or for none (in which case, 'newdata=NULL', which is the
  #           default). If some value is specified and at least one of the covariates is continuous, then a single survival
  #           curve will be computed in correspondence of the average values of all the covariates (including the factors, 
  #           which in this case are expanded into indicators). 
  #           THE ORDER OF THE VARIABLES IS list(covs,factors)
  # nsim = the number of simulations from the distribution of the survival curves. Default at nsim=1, in which case
  #          uses the point estimate for the relevant distributional parameters and computes the resulting survival curve
  # ... = additional options
  
  exArgs <- list(...)
  if(is.null(t)) {
    t <- sort(unique(fit$misc$km$time))
  }
  
  m <- fit$models[[mod]]                  # Extracts the model object from the survHE output
  dist <- fit$models[[mod]]$dlist$name    # Extracts the name of the distribution fitted
  n.elements <- ifelse(is.null(newdata),0,length(newdata))
  n.provided <- unlist(lapply(newdata,function(x) length(x)))
  
  # If no newdata are provided then see what to do!
  data <- fit$misc$data
  test <- attributes(terms(fit$misc$formula))$term.labels
  ncovs <- length(test)
  formula.temp <- as.formula(gsub("inla.surv","Surv",deparse(fit$misc$formula)))
  Xraw <- model.frame(formula.temp,data=data)
  is.fac <- sapply(Xraw, is.factor)[-1]
  w <- (which(sapply(Xraw,is.factor)==1))-1
  X <- matrix(colMeans(model.matrix(formula.temp,data=data)), nrow = 1)
  if(fit$method=="inla") {
     colnames(X) <- rownames(m$summary.fixed)
   } else {
     colnames(X) <- colnames(model.matrix(formula.temp,data=data)) #c("Intercept",test)
   }
  # newdata is not given (ie = NULL); this implies n.provided=NULL
  if (n.elements==0) {
    # If all the covariates are factors and mode_factor = False, then get survival curves for all the combinations
    if(all(is.fac) & length(is.fac)>0 ) {
      X <- unique(model.matrix(formula.temp,data=data))
      nam <- as.matrix(unique(X))
      for (i in 2:ncol(nam)) {
        nam[, i] <- paste(colnames(nam)[i],nam[, i], sep = "=")
      }
      rownames(X) <- apply(nam, 1, paste, collapse = ",")
    }
    else if (all(is.fac) & length(is.fac)>0 ){
      count <- table(apply(model.matrix(formula.temp,data=data), 1, paste, collapse=" ")) 
       col_splits <-strsplit(names(count[which.max(count)])," ")
       mode_col <- as.matrix(sapply(col_splits,as.numeric),by = row)
       X <- unique(model.matrix(formula.temp,data=data))
       mode_factor <- which(sapply(seq(1:dim(X)[1]), function(i) all(as.matrix(X[i,]) == mode_col) ))
       print(mode_factor)
       ## Is the line below redundant?
        nam <- as.matrix(unique(X))
        for (i in 2:ncol(nam)) {
        nam[, i] <- paste(colnames(nam)[i],nam[, i], sep = "=")
      }
      rownames(X) <- apply(nam, 1, paste, collapse = ",")
      X <- X[mode_factor,]
    }
  }
  # newdata is a list with a single list inside (only one scenario given by the user)
  if (n.elements==1) {
    if (n.provided!=ncovs) {
      stop("You need to provide data for *all* the covariates specified in the model, in the list 'newdata'")
    } else {
      if(all(is.fac)) {
       X <- unique(model.matrix(formula.temp,data=data))

       ## Is the line below redundant?
      nam <- as.matrix(unique(X))
      for (i in 2:ncol(nam)) {
        nam[, i] <- paste(colnames(nam)[i],nam[, i], sep = "=")
      }
      rownames(X) <- apply(nam, 1, paste, collapse = ",")
      #X <- X[mode_factor,]

        n.elements <- ifelse(dim(X)[1]>1,2,1)
      }
    }
  }
  if (n.elements>1) {
    if (!all(n.provided==ncovs)) {
      stop("You need to provide data for *all* the covariates specified in the model, in the list 'newdata'")
    } else {
      X <- matrix(rep(X,n.elements),nrow=n.elements,byrow=T)
      print(X)
      if(fit$method=="inla") {
        colnames(X) <- rownames(m$summary.fixed)
      } else {
        colnames(X) <- c("Intercept",test)
      }
      # Just like flexsurv, if you provide values for the covariates, you have to do so for *all*!
      names <- unique(unlist(lapply(newdata,function(x) names(x))))
      positions <- lapply(1:length(names),function(i) which(grepl(names[i],colnames(X))))
      print(positions)
      temp <- matrix(unlist(newdata),nrow=length(newdata),byrow=T)
      colnames(temp) <- names
      # Could change the value in X with the value in temp[-w] for the continuous variables
      contin <- (1:length(names))[-w]
      # Do this only if there're some continuous covariates!
      if (length(contin>0)) {
        for (i in 1:length(contin)) {
          for (j in 1:n.elements) {
            X[j,positions[[contin[i]]]] <- temp[j,contin[i]]
          }
        }
      }
      # And then change the value in X with the factor expansion for the categorical variables
      for (i in 1:length(w)) {
        for (j in 1:n.elements) {
          check <- eval(parse(text=paste0("temp[j,w[i]]==as.numeric(levels(as.factor(data$",names[w[i]],")))")))
          X[j,positions[[w[i]]]] <- check[-1]
        }
      }
    }
  }
  
  # If the original model(s) have been fitted using MLE & flexsurvreg, then use flexsurv::summary to compute the survival curves
  if(fit$method %in% c("mle","splines")) {
    dist <- ifelse(dist=="weibull.quiet","weibull",dist)
    S <- sim <-list()
    if(nsim==1) {
        if(is.null(newdata)) {
    df <- data.table(data)
        test <- attributes(terms(formula))$term.labels
        ncovs <- length(test)
        factors <- gsub("as.factor[( )]","",test) 
        factors <- gsub("[( )]","",factors)
        nd <- df[,.N,by=factors][N==max(N),]
        if (ncol(nd)==1) { 
          S[[1]] <- summary(m,t=t)
         } else {
      S[[1]] <- summary(m,t=t,newdata = nd) }

      sim <- NULL
  } else {
    newdatalist <- lapply(1:length(newdata),function(x) data.frame(newdata[[x]]))
    nd <- do.call("rbind", newdatalist)
    S[[1]] <- summary(m,t=t,newdata = nd)
      sim <- NULL
  }
        
    } else {
      if (is.null(newdata)){
          sim <- list(flexsurv::normboot.flexsurvreg(m,B=nsim,newdata=X))
      } else {
        sim <- lapply(1:n.elements,function(i) flexsurv::normboot.flexsurvreg(m,B=nsim,newdata=newdata[[i]]))
      }
      txt1 <- paste("x[",1:dim(sim[[1]])[2],"]",sep="",collapse=",")
      tmp <- lapply(1:length(sim), function(i) eval(parse(text=paste0("t(apply(sim[[",i,"]],1,function(x) d",dist,"(t,",txt1,")/h",dist,"(t,",txt1,")))"))))  
      # if(dist=="exp") {
      #    tmp <- lapply(1:length(sim), function(i) t(apply(sim[[i]],1,function(x) dexp(t,x[1])/hexp(t,x[1]))))
      # } else {
      #   tmp <- lapply(1:length(sim), function(i) eval(parse(text=paste0("t(apply(sim[[",i,"]],1,function(x) d",dist,"(t,x[1],x[2])/h",dist,"(t,x[1],x[2])))"))))  
      # }
      S <- list(list())
      S <- lapply(1:nsim,function(i) {
          lapply(1:length(sim),function(j) {
              cbind(t,tmp[[j]][i,])
          })
      })
    }
  } 
  
  # If the original model(s) have been fitted using INLA, then use the (summaries of the) posterior distributions to compute the survival curves
  if(fit$method=="inla") {
    # A function to rescale the parameters of a given model and then computes the survival curve
    rescale.inla <- function(m,linpred) {
      if (m$dlist$name=="weibull") {
        shape <- m$summary.hyperpar[1,1]
        scale <- exp(linpred)^(1/-shape)
        S <- lapply(1:length(scale), function(x) cbind(t,dweibull(t,shape,scale[x])/hweibull(t,shape,scale[x]))) 
      }
      if (m$dlist$name=="exponential") {
        rate <- exp(linpred)
        S <- lapply(1:length(rate), function(x) cbind(t,dexp(t,rate[x])/hexp(t,rate[x])))  
      }
      if (m$dlist$name=="loglogistic") {
        shape <- m$summary.hyperpar[1,1]
        scale <- exp(linpred)
        S <- lapply(1:length(scale), function(x) cbind(t,dllogis(t,shape,scale[x])/hllogis(t,shape,scale[x]))) 
      }
      if (m$dlist$name=="lognormal") {
        mulog <- linpred
        sdlog <- INLA::inla.contrib.sd(m)$hyper[1,1]
        S <- lapply(1:length(mulog), function(x) cbind(t,dlnorm(t,mulog[x],sdlog)/hlnorm(t,mulog[x],sdlog))) 
      }
      return(S)
    }
    
    # Now computes the survival curves for the relevant case
    if (nsim==1) {
      S <- list()
      linpred <- apply(m$summary.fixed[,1]*t(X),2,sum)
      S[[1]] <- rescale.inla(m,linpred)
      sim <- NULL
    } else {
      jpost <- suppressWarnings(INLA::inla.posterior.sample(n=nsim,m))
      pos <- pmatch(rownames(m$summary.fixed),rownames(jpost[[1]]$latent))
      sim1 <- matrix(unlist(lapply(jpost,function(x) x$latent[pos,])),ncol=length(pos),byrow=T)
      colnames(sim1) <- m$names.fixed
      if (m$nhyper>0) {
        sim2 <- matrix(unlist(lapply(jpost,function(x) x$hyperpar)),ncol=m$nhyper,byrow=T)
        sim <- cbind(sim2,sim1)
        colnames(sim) <- c(paste0("hyperpar",1:m$nhyper),m$names.fixed)
      } else {
        sim <- sim1
      }
      linpred <- matrix(unlist(lapply(1:nsim,function(i) apply(sim1[i,]*t(X),2,sum))),nrow=nsim,byrow=T)
      if(m$dlist$name=="weibull") {
        shape <- sim[,1]
        scale <- exp(linpred)^(1/-shape)
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(scale)[2],function(j) {
            cbind(t,dweibull(t,shape[i],scale[i,j])/hweibull(t,shape[i],scale[i,j]))
          })
        })
      }
      if(m$dlist$name=="exponential") {
        rate <- exp(linpred)
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(rate)[2],function(j) {
            cbind(t,dexp(t,rate[i,j])/hexp(t,rate[i,j]))
          })
        })
      }
      if (m$dlist$name=="loglogistic") {
        shape <- sim[,1]
        scale <- exp(linpred)
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(scale)[2],function(j) {
            cbind(t,dlogis(log(t),scale[i,j],1/shape[i])/hllogis(log(t),scale[i,j],helpshape[i]))
          })
        })
      }
      if (m$dlist$name=="lognormal") {
        mulog <- linpred
        sdlog <- sim[,1]
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(scale)[2],function(j) {
            cbind(t,dlnorm(t,mulog[i,j],sdlog[i])/hlnorm(t,mulog[i,j],sdlog[i]))
          })
        })
      }
    }
  }
  
  if(fit$method=="mcmc") {

    # Now computes the survival curves for the relevant case
    if (nsim==1) {
      S <- list()
      rels <- c(grep("beta",rownames(m$summary)),grep("delta",rownames(m$summary)),grep("gamma",rownames(m$summary)))
      linpred <- apply(m$summary[rels,1]*t(X),2,sum)
      
      if (m$dlist$name=="weibull") {
        shape <- mean(m$sims.matrix[,"shape"])
        scale <- exp(-linpred)
        S[[1]] <- lapply(1:length(scale),function(j) {
            cbind(t,dweibull(t,shape,scale[j])/hweibull(t,shape,scale[j]))
        })
      }
      if (m$dlist$name=="exponential") {
        rate <- exp(linpred)
        S[[1]] <- lapply(1:length(rate),function(j) {
          cbind(t,dexp(t,rate[j])/hexp(t,rate[j]))
        })
      }
      if (m$dlist$name=="loglogistic") {
        shape <- mean(m$sims.matrix[,"tau"])
        scale <- exp(linpred)
        S[[1]] <- lapply(1:length(scale),function(j) {
          cbind(t,dllogis(t,shape,scale[j])/hllogis(t,shape,scale[j]))
        })
      }
      if (m$dlist$name=="dloglogis") {
        shape <- mean(m$sims.matrix[,"shape"])
        scale <- exp(linpred)
        S[[1]] <- lapply(1:length(scale),function(j) {
          cbind(t,dllogis(t,shape,scale[j])/hllogis(t,shape,scale[j]))
        })
      }
      if (m$dlist$name=="lognormal") {
        sdlog <- mean(m$sims.matrix[,"sd"])
        mulog <- linpred
        S[[1]] <- lapply(1:length(mulog),function(j) {
          cbind(t,dlnorm(t,mulog[j],sdlog)/hlnorm(t,mulog[j],sdlog))
        })
      }
      if (m$dlist$name=="weibullPH") {
        shape <- mean(m$sims.matrix[,"shape"])
        scale  <- exp(linpred)^(1/(-shape))
        S[[1]] <- lapply(1:length(scale),function(j) {
            cbind(t,dweibull(t,shape,scale[j])/hweibull(t,shape,scale[j]))
        })
      }
      if (m$dlist$name=="gamma") {
        shape <- mean(m$sims.matrix[,"shape"])
        scale  <- exp(linpred)
        S[[1]] <- lapply(1:length(scale),function(j) {
            cbind(t,dgamma(t,shape,scale[j])/hgamma(t,shape,scale[j]))
        })
      }
        if (m$dlist$name=="gompertz") {
        shape <- mean(m$sims.matrix[,"shape"])
        scale <- exp(linpred)
        S[[1]] <- lapply(1:length(scale),function(j) {
            cbind(t,dgompertz(t,shape,1/scale[j])/hgompertz(t,shape,1/scale[j]))
        })
      }
      ### ALL THE OTHER CASES!
      sim <- NULL
    } else {
      rels <- c(grep("beta",rownames(m$summary)),grep("delta",rownames(m$summary)),grep("gamma",rownames(m$summary)))
      sim <- m$sims.matrix[,rels]
      if (rels ==1){ 
        linpred <- matrix(unlist(lapply(1:nsim,function(i) apply(sim[i]*t(X),2,sum))),nrow=nsim,byrow=T)
      }
      else{
      linpred <- matrix(unlist(lapply(1:nsim,function(i) apply(sim[i,]*t(X),2,sum))),nrow=nsim,byrow=T)
    }
      if (m$dlist$name=="weibull") {
        shape <- m$sims.matrix[,"shape"]
        scale <- exp(-linpred)
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(scale)[2],function(j) {
            cbind(t,dweibull(t,shape[i],scale[i,j])/hweibull(t,shape[i],scale[i,j]))
          })
        })
      }
      if (m$dlist$name=="exponential") {
        rate <- exp(linpred)
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(rate)[2],function(j) {
            cbind(t,dexp(t,rate[i,j])/hexp(t,rate[i,j]))
          })
        })
      }
      if (m$dlist$name=="loglogistic") {
        shape <- m$sims.matrix[,"tau"]
        scale <- exp(linpred)
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(mulog)[2],function(j) {
            cbind(t,dllogis(t,shape[i],scale[i,j])/hllogis(t,shape[i],scale[i,j]))
          })
        })
      }
      if (m$dlist$name=="dloglogis") {
        shape <- mean(m$sims.matrix[,"shape"])
        scale <- exp(linpred)
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(scale)[2],function(j) {
            cbind(t,dllogis(t,shape[i],scale[i,j])/hllogis(t,shape[i],scale[i,j]))
          })
        })
      }
      if (m$dlist$name=="lognormal") {
        sdlog <- m$sims.matrix[,"sd"]
        mulog <- linpred
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(mulog)[2],function(j) {
            cbind(t,dlnorm(t,mulog[i,j],sdlog[i])/hlnorm(t,mulog[i,j],sdlog[i]))
          })
        })
      }
      if (m$dlist$name=="weibullPH") {
        shape <- m$sims.matrix[,"shape"]
        scale  <- exp(linpred)^(1/(-shape[1:nsim]))
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(scale)[2],function(j) {
            cbind(t,dweibull(t,shape[i],scale[i,j])/hweibull(t,shape[i],scale[i,j]))
          })
        })
      }
      ### needs to be fixed
      if (m$dlist$name=="gamma") {
        shape <- m$sims.matrix[,"shape"]
        scale  <- exp(linpred)
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(scale)[2],function(j) {
            cbind(t,1-pgamma(t,shape[i],scale[i,j]))
          })
        })
      }
      if (m$dlist$name=="gompertz") {
        shape <- m$sims.matrix[,"shape"]
        scale  <- exp(linpred)
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(scale)[2],function(j) {
            cbind(t,dgompertz(t,shape[i],scale[i,j])/hgompertz(t,shape[i],scale[i,j]))
          })
        })
      }
    }
      ## ALL OTHER CASES!!!

  }
  n.elements <- length(S[[1]]) 
  mat <- lapply(1:n.elements,function(j) matrix(unlist(lapply(1:nsim,function(i) S[[i]][[j]][,2])),nrow=nsim,byrow=T))
  # Now defines the output of the function
  list(S=S,sim=sim,nsim=nsim,mat=mat)
}



print.survHE <- function(x,mod=1,...) {
  # Creates a print method for the objects in the class survHE
  # x is the survHE object (the output of the call to fit.models)
  # mod is the index of the model. Default value is 1, but the user can choose which model fit to visualise, 
  #     if the call to fit.models has a vector argument for distr (so many models are fitted & stored in the same object)
  # ... optional arguments
  # digits = number of significant digits to be shown in the summary table (default = 6)
  # nsim = number of simulations from the joint posterior for INLA (default = 100)
  # original = a flag to say whether the *original* table from either INLA or MCMC should be printed
  
  exArgs <- list(...)
  if(exists("original",where=exArgs)) {original=exArgs$original} else {original=FALSE}
  # Can select the number of digits to be printed in the output table
  if(!exists("digits",where=exArgs)){digits=6} else {digits=exArgs$digits}

  if(x$method %in% c("mle","splines")) {
    res <- x$models[[mod]]$res[,c(1,4,2,3)]
    if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  }
  if(x$method=="inla" & original==FALSE) {
    # Rescales the parameters to make the estimates comparable with flexsurvreg
    if(!exists("nsim",where=exArgs)){nsim <- 100} else {nsim=exArgs$nsim} 

    # This is a rescaling function for the built-in models (that INLA can do by default)
    rescale.print.inla <- function(x,mod,nsim) {
      # Simulates from the joint posterior of *all* parameters & hyperparameters
      jpost <- suppressWarnings(INLA::inla.posterior.sample(n=nsim,x$models[[mod]]))
      # This finds the position of the hyperparameters in the simulations from the joint posterior
      pos <- pmatch(rownames(x$models[[mod]]$summary.fixed),rownames(jpost[[1]]$latent))
      
      if(x$models[[mod]]$dlist=="weibull") {
        shape <- unlist(lapply(jpost,function(x) x$hyperpar)); names(shape) <- NULL
        scale <- exp(unlist(lapply(jpost,function(x) x$latent[pos[1],])))^(1/-shape)
        effects <- matrix(NA,nrow=(length(pos)-1),ncol=nsim)
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          for (j in 2:length(pos)) {
            effects[(j-1),] <- log(exp(unlist(lapply(jpost,function(x) x$latent[pos[j],])))^(1/-shape))
          }
          rownames(effects) <- x$models[[mod]]$names.fixed[-1]
        }
        tab <- rbind(shape,scale,effects)
      }
      if(x$models[[mod]]$dlist=="exponential") {
        rate <- exp(unlist(lapply(jpost,function(x) x$latent[pos[1],])))
        effects <- matrix(NA,nrow=(length(pos)-1),ncol=nsim)
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          for (j in 2:length(pos)) {
            effects[(j-1),] <- unlist(lapply(jpost,function(x) x$latent[pos[j],]))
          }
          rownames(effects) <- x$models[[mod]]$names.fixed[-1]
        }
        tab <- rbind(rate,effects)
      }
      if(x$models[[mod]]$dlist=="lognormal") {
        prec <- unlist(lapply(jpost,function(x) x$hyperpar)); names(prec) <- NULL
        sdlog <- 1/sqrt(prec)
        meanlog <- unlist(lapply(jpost,function(x) x$latent[pos[1],]))
        effects <- matrix(NA,nrow=(length(pos)-1),ncol=nsim)
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          for (j in 2:length(pos)) {
            effects[(j-1),] <- unlist(lapply(jpost,function(x) x$latent[pos[j],]))
          }
          rownames(effects) <- x$models[[mod]]$names.fixed[-1]
        }
        tab <- rbind(meanlog,sdlog,effects)
      }
      if(x$models[[mod]]$dlist=="loglogistic") {
        shape <- unlist(lapply(jpost,function(x) x$hyperpar)); names(shape) <- NULL
        scale <- exp(unlist(lapply(jpost,function(x) x$latent[pos[1],])))
        effects <- matrix(NA,nrow=(length(pos)-1),ncol=nsim)
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          for (j in 2:length(pos)) {
            effects[(j-1),] <- unlist(lapply(jpost,function(x) x$latent[pos[j],]))
          }
          rownames(effects) <- x$models[[mod]]$names.fixed[-1]
        }
        tab <- rbind(shape,scale,effects)
      }
      return(tab)
    }
    # The user could specify a rescale.print function for their own specific model and that would be used instead
    if(exists("rescale.print",where=exArgs)) {
      func <- exArgs$rescale.print
      if(exists("inputs",where=exArgs)) {
        inputs=exArgs$inputs
      } else {
        inputs=list()
      }
    } else {
        func <- rescale.print.inla
        inputs <- list(x,mod,nsim)
    }
    tab <- do.call(what=func,args=inputs)
    
    res <- t(apply(tab,1,function(x) c(mean(x),sd(x),quantile(x,.025),quantile(x,.975))))
    if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  }

  if(x$method=="mcmc" & original==FALSE) {
    ### CREATE A FUNCTION AND THEN GENERALISE, JUST LIKE FOR INLA!!
    ### ALSO CHECK WHAT HAPPENS WHEN THERE'RE NO COVARIATES!!!!
    # Sets nsim to the number of available MCMC runs
    nsim <- x$models[[mod]]$n.sims
    D <- x$misc$vars$nlevs

    if (x$models[[mod]]$dlist$name=="weibull") {
      nodes <- rownames(x$models[[1]]$summary)
      shape <- x$models[[mod]]$sims.list$shape
      scale <- exp(-x$models[[mod]]$sims.list$beta)
      # Searches for the factors & covariates
      effects <- matrix(nrow=0,ncol=nsim)
      which.facts <- grep("delta",colnames(x$models[[1]]$sims.matrix))
      if (length(which.facts)>0) {
        names.facts <- x$misc$vars$factors
        effects.facts <- t(log(exp(-x$models[[mod]]$sims.matrix[,which.facts])))
        rownames(effects.facts) <- unlist(lapply(1:length(names.facts),function(i) paste0("as.factor(",names.facts[i],")",(2:D[i]))))
        effects <- rbind(effects,effects.facts)
      }
      which.covs <- grep("gamma",colnames(x$models[[1]]$sims.matrix))
      if (length(which.covs)>0) {
        names.covs <- x$misc$vars$covs
        effects.covs <- t(log(exp(-x$models[[mod]]$sims.matrix[,which.covs])))
        rownames(effects.covs) <- names.covs
        effects <- rbind(effects,effects.covs)
      }
      tab <- rbind(shape,scale,effects)
    }
    if (x$models[[mod]]$dlist$name=="exponential") {
      nodes <- rownames(x$models[[1]]$summary)
      rate <- exp(x$models[[mod]]$sims.list$beta)
      # Searches for the factors & covariates
      effects <- matrix(nrow=0,ncol=nsim)
      which.facts <- grep("delta",colnames(x$models[[1]]$sims.matrix))
      if (length(which.facts)>0) {
        names.facts <- x$misc$vars$factors
        effects.facts <- t(log(exp(x$models[[mod]]$sims.matrix[,which.facts])))
        rownames(effects.facts) <- unlist(lapply(1:length(names.facts),function(i) paste0("as.factor(",names.facts[i],")",(2:D[i]))))
        effects <- rbind(effects,effects.facts)
      }
      which.covs <- grep("gamma",colnames(x$models[[1]]$sims.matrix))
      if (length(which.covs)>0) {
        names.covs <- x$misc$vars$covs
        effects.covs <- t(log(exp(x$models[[mod]]$sims.matrix[,which.covs])))
        rownames(effects.covs) <- names.covs
        effects <- rbind(effects,effects.covs)
      }
      tab <- rbind(rate,effects)
    }
    if (x$models[[mod]]$dlist$name=="loglogistic") {
      nodes <- rownames(x$models[[1]]$summary)
      shape <- (x$models[[mod]]$sims.list$tau)
      scale <- exp(x$models[[mod]]$sims.list$beta)
      # Searches for the factors & covariates
      effects <- matrix(nrow=0,ncol=nsim)
      which.facts <- grep("delta",colnames(x$models[[1]]$sims.matrix))
      if (length(which.facts)>0) {
        names.facts <- x$misc$vars$factors
        effects.facts <- t(log(exp(x$models[[mod]]$sims.matrix[,which.facts])))
        rownames(effects.facts) <- unlist(lapply(1:length(names.facts),function(i) paste0("as.factor(",names.facts[i],")",(2:D[i]))))
        effects <- rbind(effects,effects.facts)
      }
      which.covs <- grep("gamma",colnames(x$models[[1]]$sims.matrix))
      if (length(which.covs)>0) {
        names.covs <- x$misc$vars$covs
        effects.covs <- t(log(exp(x$models[[mod]]$sims.matrix[,which.covs])))
        rownames(effects.covs) <- names.covs
        effects <- rbind(effects,effects.covs)
      }
      tab <- rbind(shape,scale,effects)
    }
    if (x$models[[mod]]$dlist$name=="lognormal") {
      nodes <- rownames(x$models[[1]]$summary)
      sdlog <- x$models[[mod]]$sims.list$sd
      meanlog <- x$models[[mod]]$sims.list$beta
      # Searches for the factors & covariates
      effects <- matrix(nrow=0,ncol=nsim)
      which.facts <- grep("delta",colnames(x$models[[1]]$sims.matrix))
      if (length(which.facts)>0) {
        names.facts <- x$misc$vars$factors
        effects.facts <- t(log(exp(x$models[[mod]]$sims.matrix[,which.facts])))
        rownames(effects.facts) <- unlist(lapply(1:length(names.facts),function(i) paste0("as.factor(",names.facts[i],")",(2:D[i]))))
        effects <- rbind(effects,effects.facts)
      }
      which.covs <- grep("gamma",colnames(x$models[[1]]$sims.matrix))
      if (length(which.covs)>0) {
        names.covs <- x$misc$vars$covs
        effects.covs <- t(log(exp(x$models[[mod]]$sims.matrix[,which.covs])))
        rownames(effects.covs) <- names.covs
        effects <- rbind(effects,effects.covs)
      }
      tab <- rbind(meanlog,sdlog,effects)
    }
    res <- t(apply(tab,1,function(x) c(mean(x),sd(x),quantile(x,.025),quantile(x,.975))))
    if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  }
  
  # Finally creates the table
  
  # Original formatting of the tables from INLA & BUGS
  if(original==TRUE) {
    if (x$method=="mle") {
      print(x$models[[mod]])
    }
    if (x$method=="inla") {
      print(summary(x$models[[mod]]))
    }
    if (x$method=="mcmc") {
      if (x$models[[1]]$dlist$name=="loglogistic") {
        x$models[[mod]]$summary <- x$models[[mod]]$summary[-grep("lambda",rownames(x$models[[mod]]$summary)),]
      }
      print(x$models[[mod]],digits=digits)
    }
  } else {
    # FORMATS THE TABLE
    # Now recodes the model name to a standardised string
    if(x$models[[mod]]$dlist$name=="exp" | x$models[[mod]]$dlist$name=="exponential") {label="Exponential"}
    if(x$models[[mod]]$dlist$name=="gamma") {label="Gamma"}
    if(x$models[[mod]]$dlist$name=="lognormal" | x$models[[mod]]$dlist$name=="lnorm") {label="log-Normal"}
    if(x$models[[mod]]$dlist$name=="llogis" | x$models[[mod]]$dlist$name=="loglogistic") {label="log-Logistic"}
    if(x$models[[mod]]$dlist$name=="gengamma") {label="Generalised Gamma"}
    if(x$models[[mod]]$dlist$name=="weibull" | x$models[[mod]]$dlist$name=="weibull.quiet" | x$models[[mod]]$dlist$name=="weibullPH") {label="Weibull"}
    if(x$models[[mod]]$dlist$name=="genf") {label="Generalised F"}
    if(x$models[[mod]]$dlist$name=="gompertz") {label="Gompertz"}
    if(x$method == "splines") {label = "flexible parametric distribution"}
    if(x$method=="mle") {label.met="Flexsurvreg \n(Maximum Likelihood Estimate)"}
    if(x$method=="inla") {label.met="INLA (Bayesian inference via \nIntegrated Nested Laplace Approximation)"}
    if(x$method=="mcmc") {label.met="OpenBUGS (Bayesian inference via \nMarkov Chain Monte Carlo)"}
    if(x$method=="splines") {label.met="Flexsurvreg \n(Maximum Likelihood Estimate)"}

    cat("\n")
    cat(paste0("Model fit for the ",label," model, obtained using ",label.met,". Running time: ",
               format(x$misc$time2run[[mod]],digits=5,nsmall=3)," seconds"))
    cat("\n\n")
    print(res,quote=F,digits=digits,justify="center")
    cat("\n")
    cat("Model fitting summaries\n")
    cat(paste0("Akaike Information Criterion (AIC)....: ",format(x$model.fitting$aic[[mod]],digits=6,nsmall=3)))
    cat("\n")
    cat(paste0("Bayesian Information Criterion (BIC)..: ",format(x$model.fitting$bic[[mod]],digits=6,nsmall=3)))
    if(x$method=="inla" | x$method=="mcmc") {
      cat("\n")
      cat(paste0("Deviance Information Criterion (DIC)..: ",format(x$model.fitting$dic[[mod]],digits=6,nsmall=3)))
    }
    cat("\n\n")
  }
}


psa.plot <- function(psa,...) {
  # Plots the survival curves for all the PSA simulations
  # psa = the result of the call to the function make.surv
  # ... = additional arguments
  # xlab = label for the x-axis
  # ylab = label for the y-axis
  # col = vector of colours with which to plot the curves
  # alpha = parameter to determine the transparency (default = 0.1)
  
  n.elements <- length(psa$S[[1]]) 
  times <- psa$S[[1]][[1]][,1]
  exArgs <- list(...)
  if(!exists("xlab",where=exArgs)) {xlab="Time"} else {xlab=exArgs$xlab}
  if(!exists("ylab",where=exArgs)) {ylab="Survival"} else {ylab=exArgs$ylab}
  if(!exists("col",where=exArgs)) {col=sample(colors(),n.elements)} else {col=exArgs$col}
  if(!exists("alpha",where=exArgs)) {alpha=0.1} else {alpha=exArgs$alpha}
  if(!exists("main",where=exArgs)) {main=""} else {main=exArgs$main}
  
  # If there's only the average value for the survival curve, simpler plot
  if (psa$nsim==1) {
    alpha <- 1
    plot(psa$S[[1]][[1]][,1:2],t="l",xlab=xlab,ylab=ylab,col=adjustcolor(col[1],alpha),ylim=c(0,1),xlim=range(pretty(times)),
         main=main,axes=F)
    if (n.elements>1) {
        pts2 <- lapply(2:n.elements,function(i) points(psa$S[[1]][[i]],t="l",col=adjustcolor(col[i],alpha)))
    }
  }
  
  # If there are nsim simulations from the survival curves, then more complex plot
  if (psa$nsim>1) {
    tmp <- lapply(1:n.elements,function(j) matrix(unlist(lapply(1:psa$nsim,function(i) psa$S[[i]][[j]][,2])),nrow=psa$nsim,byrow=T))
    q025 <- lapply(1:n.elements, function(j) apply(tmp[[j]],2,function(x) quantile(x,.025)))
    q500 <- lapply(1:n.elements, function(j) apply(tmp[[j]],2,function(x) quantile(x,.5))) 
    q975 <- lapply(1:n.elements, function(j) apply(tmp[[j]],2,function(x) quantile(x,.975))) 
    print(psa$S[[1]][[1]][,1])
    print(q500[[1]])
    print(adjustcolor(col[1],1))
    #plot(psa$S[[1]][[1]][,1],q500[[1]],col=adjustcolor(col[1],1),t="l",xlab=xlab,ylab=ylab,ylim=c(0,1),xlim=range(pretty(times)),lwd=2,main=main,axes=F)
    points(psa$S[[1]][[1]][,1],q500[[1]],col=adjustcolor(col[1],1),t="l",lwd=2,main=main)

    polygon(c(psa$S[[1]][[1]][,1],rev(psa$S[[1]][[1]][,1])),c(q975[[1]],rev(q025[[1]])),col=adjustcolor(col[1],alpha),border=NA)
    if (n.elements>1) {
        lapply(2:n.elements, function(i) {
            pts1 <- points(psa$S[[1]][[i]][,1],q500[[i]],col=adjustcolor(col[i],1),t="l",lwd=2) 
            pts2 <- polygon(c(psa$S[[1]][[i]][,1],rev(psa$S[[1]][[i]][,1])),c(q975[[i]],rev(q025[[i]])),col=adjustcolor(col[i],alpha),border=NA)
        })
    }
  }
  axis(1);axis(2)
}


plot.survHE <- function(x,...) {
  ## Plots the KM + the results of the model fitted by fit.models()
  ## Uses different commands, depending on which method has been used to fit the models
  #
  # x = the result of the call to the fit.model function
  #
  # ... = additional (mainly graphical) options
  # xlab
  # ylab
  # lab.trt
  # cex.trt
  # n.risk
  # xlim
  # colors
  # labs
  # add.km = TRUE (whether to also add the Kaplan Meier estimates of the data)
  
  exArgs <- list(...) 		# Lists all the additional inputs
  
  mod <- x$models
  nmodels <- length(mod)  # Number of models fitted by fit.models()
  
  # Checks that extra options are specified
  if (!exists("t",where=exArgs)) {t=sort(unique(x$misc$km$time))} else {t=exArgs$t}
  if (!exists("xlab",where=exArgs)) {xl="time"} else {xl=exArgs$xlab}
  if (!exists("ylab",where=exArgs)) {yl="Survival"} else {yl=exArgs$ylab}
  if (!exists("lab.trt",where=exArgs)) {lab.trt=names(x$misc$km$strata)} else {lab.trt=names(x$km$strata)<-exArgs$lab.trt}
  if (!exists("cex.trt",where=exArgs)) {cex.trt=.8} else {cex.trt=exArgs$cex.trt}
  if (!exists("n.risk",where=exArgs)) {nrisk=FALSE} else {nrisk=exArgs$n.risk}
  if (!exists("newdata",where=exArgs)) {newdata = NULL} else {newdata=exArgs$newdata}
  if (!exists("xlim",where=exArgs) & !exists("t",where=exArgs)) {
    xlm=range(pretty(x$misc$km$time))
  } 
  if (!exists("xlim",where=exArgs) & exists("t",where=exArgs)) {
    xlm=range(pretty(t))
  }
  if (exists("xlim",where=exArgs) & !exists("t",where=exArgs)) {
    xlm <- exArgs$xlim
  }
  if (exists("xlim",where=exArgs) & exists("t",where=exArgs)) {
    xlm <- exArgs$xlim
  }
  if (!exists("colors",where=exArgs)) {
    if (nmodels>1) {colors=(2:(nmodels+1))} else {colors=2}
  } else {colors=exArgs$colors}
  if (!exists("labs",where=exArgs)) {
    labs <- sapply(mod,function(m) m$dlist$name)
    labs[pmatch("weibull.quiet",labs)] <- "Weibull"
    labs[pmatch("weibull",labs)] <- "Weibull"
    labs[pmatch("exp",labs)] <- "Exponential"
    labs[pmatch("exponential",labs)] <- "Exponential"
    labs[pmatch("gamma",labs)] <- "Gamma"
    labs[pmatch("lnorm",labs)] <- "log-Normal"
    labs[pmatch("lognormal",labs)] <- "log-Normal"
    labs[pmatch("llogis",labs)] <- "log-Logistic"
    labs[pmatch("loglogistic",labs)] <- "log-Logistic"
    labs[pmatch("gengamma",labs)] <- "Gen. Gamma"
    labs[pmatch("dloglogis",labs)] <- "log-Logistic"

  } else {labs=exArgs$labs}
  labs <- c("Kaplan Meier",labs)
  if (!exists("add.km",where=exArgs)) {add.km=TRUE} else {add.km=exArgs$add.km}
  
  # Now plots the KM curve using "rms" if add.km is set to TRUE
  if (add.km==TRUE) {
    rms::survplot(x$misc$km,                                     # Specialised plot from "rms" 
                  xlab=xl,ylab=yl,		                         # x- and y- labels
                  label.curves=list(labels=lab.trt,cex=cex.trt), # specifies curve labels
                  n.risk=nrisk,   	                             # tells R to show number at risk 
                  lwd=2,xlim=xlm 	                             # defines the size of the lines (2 pts)
    )
  } else {
    labs <- labs[-1]
    colors <- colors-1
    plot(0,0,col="white",xlab=xl,ylab=yl,axes=F,xlim=xlm,ylim=c(0,1));axis(1);axis(2);
  }
  res <- lapply(1:nmodels,function(i) make.surv(x,nsim=1,t=t,mod=i,newdata = newdata))
  for (i in 1:nmodels) {
    pts <- lapply(res[[i]]$S[[1]],function(m) cbind(m[,1],m[,2]))
    lapply(1:length(pts), function(x) points(pts[[x]],t="l",col=colors[i],lty=x))
  }
  legend(x="topright",legend=labs,lwd = 2, bty = "n",col = c("black",colors),cex=.8)
}


bugs2survHE <- function(model,distr,data,formula) {
  # Utility function to format a BUGS object (eg a model run independently of survHE) to an object that can 
  # be post-processed using the functionalities of survHE (ie to compute, plot and store survival curves etc)
  #
  # model = a (list of) BUGS object(s), including the output of the call to BUGS
  # distr = a (vector of) string(s) with the name(s) of the distribution(s) implemented in the BUGS model
  # data = a *named* data frame containing the original data fitted using BUGS
  # formula = a flexsurv-like formula to describe the model fitted using BUGS
  
  # If the user only specifies one single model, then turn this into a list
  if(class(model)!="list") {
    model <- list(model)
  }
  # Sets 'method' to 'mcmc'
  method <- "mcmc"
  
  # Defines the survHE object and starts filling it
  obj <- list(list())
  obj$models <- lapply(1:length(model), function(x) model[[x]])
  
  # Needs to determine the distribution implemented in 'model' and write it to the element $models[[1]]$dlist$name
  for (i in 1:length(model)) {
    obj$models[[i]]$dlist$name <- distr[i]
  }
  
  # Creates DIC & AIC for all the models in the list
  obj$model.fitting$dic <- unlist(lapply(1:length(model),function(x) model[[x]]$DIC))
  obj$model.fitting$aic <- unlist(lapply(1:length(model),function(x) {model[[x]]$pD+model[[x]]$summary["deviance",1]}))

  # Writes the method to "mcmc"
  obj$method <- method
  
  # The running time is not provided by the standard BUGS output, so set to NA
  obj$misc$time2run <- rep(NA,length(model))
  obj$misc$formula <- formula
  obj$misc$data <- data
  
  # Reconstructs the vars list based on the formula
  test <- attributes(terms(formula))$term.labels
  ncovs <- length(test)
  formula.temp <- as.formula(gsub("inla.surv","Surv",deparse(formula)))
  time <- all.vars(formula.temp,data)[1]
  event <- all.vars(formula.temp,data)[2]
  if (length(ncovs>1)) {
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
  } else {
    covs <- factors <- NULL
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
      if (method=="MCMC" | method=="mcmc") {
        levels(data[,cols[i]]) <- 1:nlevs
      }
      D[i] <- nlevs
    }
  }
  vars <- list(time=time,event=event,factors=factors,covs=covs,nlevs=D)
  
  # Need to create a formula for the KM that is in the right format (flexsurv-like)
  chk <- is.na(pmatch("Surv",attributes(terms(formula))$variables[[2]][1])) # If TRUE, then it's in inla.surv() terms
                                                                            # If FALSE, then formula is in Surv() terms
  # If it's in inla.surv() needs to create a Surv() formula for the KM + check that the method isn't mle
  if(chk) {
    tmp <- deparse(formula)
    km.formula <- as.formula(gsub("inla.surv","Surv",tmp))
    # If the method was originally INLA but has been modified or if the method is mle but the formula is in the wrong
    # format, then need to change the formula to be consistent with the Surv() notation
    if (method=="mle" | method=="MLE") {
      formula <- km.formula
    }
  } else {
    km.formula <- formula
    # If the method isn't inla but the formula is in Surv() terms, then change the formula
    if (method=="inla" | method=="INLA") {
      tmp <- deparse(formula)
      formula <- as.formula(gsub("Surv","inla.surv",tmp))
    }
  }
  # Computes the Kaplan Meier curve using the package "rms"
  ObjSurvfit=rms::npsurv(        # Uses the function "npsurv" from the package "rms"
    formula=km.formula,          # to fit the model specified in the "formula" object
    data=data                    # to the dataset named "data"
  )
  
  obj$misc$vars <- vars
  obj$misc$km <- ObjSurvfit
  
  # NB The BIC needs the value of n to be computed: 
  obj$model.fitting$bic <- unlist(lapply(1:length(model),function(x) {
    model[[x]]$summary["deviance",1]+model[[x]]$pD*log(sum(ObjSurvfit$n))
  }))
  
  # Formats the resulting object as an element of the class 'survHE'
  class(obj) <- "survHE"
  return(obj)
}



model.fit.plot <- function(fit,type="aic",...) {
  ## Plots a summary of the model fit for all the models tested by flexsurvreg
  # fit a flexsurv object containing the results of the analysis for several models
  # type = should the AIC or the BIC plotted? (values = "aic" or "bic")
  # ... additional (mainly graphical) options
  
  # Checks that the models to be used have a label and if not creates one
  exArgs <- list(...)
  if (!exists("models",where=exArgs)) {
    mod <- fit$models
    models <- sapply(mod,function(m) m$dlist$name)
    models[pmatch("weibull.quiet",models)] <- "Weibull"
    models[pmatch("weibull",models)] <- "Weibull"
    models[pmatch("exp",models)] <- "Exponential"
    models[pmatch("exponential",models)] <- "Exponential"
    models[pmatch("gamma",models)] <- "Gamma"
    models[pmatch("lnorm",models)] <- "log-Normal"
    models[pmatch("lognormal",models)] <- "log-Normal"
    models[pmatch("llogis",models)] <- "log-Logistic"
    models[pmatch("loglogistic",models)] <- "log-Logistic"
    models[pmatch("gengamma",models)] <- "Gen. Gamma"
  } else {models=exArgs$models}
  
  # Defines the data to be plotted
  if (type=="aic" | type=="AIC" | type=="a" | type=="A") {
    mf <- data.frame(model=models,AIC=fit$model.fitting$aic)
    lab.type <- "AIC"
  } else if (type=="bic" | type=="BIC" | type=="b" | type=="B") {
    mf <- data.frame(model=models,BIC=fit$model.fitting$bic)
    lab.type <- "BIC"
  } else if (type=="dic" | type=="DIC" | type=="d" | type=="D") {
    mf <- data.frame(model=models,DIC=fit$model.fitting$dic)
    lab.type <- "DIC"
  }
  
  # Finally do the plot
  if (!exists("xlim",where=exArgs)) {xlm=range(pretty(mf[,2]))} else {xlm=exArgs$xlim}
  if (!exists("digits",where=exArgs)) {digits=7} else {digits=exArgs$digits}
  if (!exists("nsmall",where=exArgs)) {nsmall=3} else {nsmall=exArgs$nsmall}
  if (!exists("main",where=exArgs)) {main=paste0("Model comparison based on ",lab.type)} else {main=exArgs$main}
  if (!exists("mar",where=exArgs)) {mar=c(4,6,3,1.3)} else {mar=exArgs$mar}
  if (!exists("cex.names",where=exArgs)) {cex.names=0.8} else {cex.names=exArgs$cex.names}
  par(mar=mar)                                           # Bottom,left,top & right margins
  b <- barplot(                                          # Function to draw a barplot (see BMS NICE submission)
    mf[,2],  	                                           # Makes a barplot using the values of the AIC or BIC
    names.arg=mf$model,	                                 # Names of the models (can be formatted differently)
    xlab=lab.type,                                       # Label for the x-axis
    xlim=xlm,
    xpd=F,                                               # Graphical parameter to clip at the lowest end of the range
    horiz=T,                                             # Plots the graph horizontally (better readability)
    las=1,                                               # Rotates the labels on the y-axis (better readability)
    cex.names=cex.names                                # Rescales the labels on the y-axis to 80% of normal size
    #main=main
  )
  # And then adds the actual value of the AIC/BIC for each of the models
  text(mf[,2],		                                       # Position of the text on the x-axis
       b,                                                # Position of the text on the y-axis
       format(mf[,2],digits=digits,nsmall=nsmall),       # Formats the values of the AICs/BICs/DICs, using 3 dp
       pos=4,                                            # Puts the text to the right of the bars
       cex=.8                                            # Rescales the labels on the y-axis to 80% of normal size
  )
}

test.linear.assumptions <- function(fit,mod=1,coxph=TRUE,label = FALSE){
  m <- fit$models[[mod]]
  dist <- m$dlist$name
  split_vector <- c(1)
  for(i in 2:length(fit$misc$km$time)){
  if(fit$misc$km$time[i]<fit$misc$km$time[i-1]){
  split_vector <- c(split_vector,i-1,i)
  }
  }
  split_vector <- c(split_vector,length(fit$misc$km$time))
  split_mat <- matrix(split_vector,length(split_vector)/2,2,byrow = T)
  times <- lapply(1:dim(split_mat)[1],function(x) fit$misc$km$time[split_mat[x,][1]:split_mat[x,][2]])
  survs <- lapply(1:dim(split_mat)[1],function(x) fit$misc$km$surv[split_mat[x,][1]:split_mat[x,][2]])
  if (dist %in% c("Exponential","exp","exponential")){
  plot(0,0,col="white",xlab='time',ylab='log(S(t))',axes=F,xlim=range(pretty(fit$misc$km$time)))
  axis(1)
  axis(2)
  pts <- lapply(1:dim(split_mat)[1],function(m) cbind(times[[m]],log(survs[[m]])))
    lapply(1:length(pts), function(x) points(pts[[x]],t="l",lty=x))
    if (label){legend('topright','Exponential distributional assumption')}

  #text(max(pts[[1]][,1]),max(pts[[1]][,2]), 'Exponential linear assumption', cex=0.6, col="red")

  }
  if (dist %in% c("weibull","weibullPH","weibull.quiet")){
  plot(0,0,col="white",xlab='log(time)',ylab='log(-log(S(t)))',axes=F,xlim=range(pretty(log(fit$misc$km$time))), ylim =range(pretty(log(-log(survs[[1]])))))
  axis(1)
  axis(2)
  pts <- lapply(1:dim(split_mat)[1],function(m) cbind(log(times[[m]]),log(-log(survs[[m]]))))
    lapply(1:length(pts), function(x) points(pts[[x]],t="l",lty=x))
        if (label){legend('topright','Weibull distributional assumption')}

  }
  if (dist == "llogis"){
  plot(0,0,col="white",xlab='time',ylab='log(S(t)/(1-S(t)))',axes=F,xlim=range(pretty(log(fit$misc$km$time))), ylim =range(pretty(log(survs[[1]]/(1-survs[[1]])))))
  axis(1)
  axis(2)
  pts <- lapply(1:dim(split_mat)[1],function(m) cbind(log(times[[m]]),log(survs[[m]]/(1-survs[[m]]))))
    lapply(1:length(pts), function(x) points(pts[[x]],t="l",lty=x))
    if (label){legend('topright','log-Logistic distributional assumption')}

    }
    if (dist %in% c("lognormal","lnorm")){
      ### add log normal 
  plot(0,0,col="white",xlab='time',ylab='log(S(t))',axes=F,xlim=range(pretty(log(fit$misc$km$time)))   )#, ylim =range(qnorm(1-survs[[1]])))
  axis(1)
  axis(2)
  pts <- lapply(1:dim(split_mat)[1],function(m) cbind(times[[m]],qnorm(1-survs[[m]])))
  print(pts)
    lapply(1:length(pts), function(x) points(pts[[2]],t="l",lty=x))
    if (label){legend('topright','lognormal distributional assumption')}

  }
  if (dist == "gompertz"){
  estimate.h <- function(s,t){
  denom <- t-c(t[-1],max(t)+1)
  print(denom)
  numerator <- log(s) - log(c(s[-1],0))
  print(numerator)
  return(-numerator/denom)
}
  plot(0,0,col="white",xlab='log(time)',ylab='h(t)',axes=F,xlim=range(pretty(fit$misc$km$time)), ylim =range(pretty(estimate.h(survs[[1]],times[[1]]))))
  axis(1)
  axis(2)
  pts <- lapply(1:dim(split_mat)[1],function(m) data.table(cbind(times[[m]],estimate.h(survs[[m]],times[[m]])))[V2!=0,])
    lapply(1:length(pts), function(x) points(pts[[x]],t="l",lty=x))
        if (label){legend('topright','Gompertz distributional assumption')}

  }
}

# ################################################################################################################################
# ## DO WE STILL NEED THIS????
# extrapolate.plot <- function(object,ci=FALSE,...) {
#   ## Plots the extrapolated curves 
#   # object = the output of the extrapolate.surv function
#   # ci = should confidence bands be plotted too?
#   # ... = additional options (eg graphical parameters)
#   
#   if(is.element("flexsurv",installed.packages()[,1])==F){install.packages("flexsurv")}
#   
#   exArgs <- list(...)
#   if (!exists("xlab",where=exArgs)) {xl="time"} else {xl=exArgs$xlab}
#   if (!exists("ylab",where=exArgs)) {yl="Survival"} else {yl=exArgs$ylab}
#   if (!exists("lab.trt",where=exArgs)) {lab.trt=paste0("Intervention ",1:length(object))} else {lab.trt<-exArgs$lab.trt}
#   if (!exists("colors",where=exArgs)) {
#     if (length(object)>1) {colors=seq(1:length(object))} else {colors=exArgs$colors}
#   }
#   if (!exists("main",where=exArgs)) {main="Extrapolated survival curve"} else {main=exArgs$main}
#   
#   if (ci==FALSE || ci==F) {
#     plot(object[[1]]$time,object[[1]]$est,t="l",xlab=xl,ylab=yl,col=colors[1],main=main)
#     if (length(object)>1) {
#       for (i in 1:length(object)) {
#         points(object[[i]]$time,object[[i]]$est,t="l",col=colors[i])
#       }
#       legend(x="topright",legend=lab.trt,lwd=2,bty="n",col=colors,cex=.8)
#     }
#   }
#   if (ci==TRUE || ci==T) {
#     plot(object[[1]]$time,object[[1]]$est,t="l",xlab=xl,ylab=yl,col=colors[1],lwd=2,main=main)
#     polygon(c(rev(object[[1]]$time),object[[1]]$time),c(rev(object[[1]]$lcl),object[[1]]$ucl),col="grey80",border=NA)
#     points(object[[1]]$time,object[[1]]$est,t="l",xlab=xl,ylab=yl,col=colors[1],lwd=2)
#     if (length(object)>1) {
#       for (i in 2:length(object)) {
#         polygon(c(rev(object[[i]]$time),object[[i]]$time),c(rev(object[[i]]$lcl),object[[i]]$ucl),col="grey80",border=NA)
#         points(object[[i]]$time,object[[i]]$est,t="l",col=colors[1],lwd=2,lty=i)
#       }
#       legend(x="topright",legend=lab.trt,lwd=2,bty="n",lty=1:length(object),cex=.8)
#     }
#   }
# }
################################################################################################################################


write.surv <- function(object,file,sheet) {
  # Writes the survival summary to an excel file (helpful to then call the values in the Markov model)
  # object = a summary.flexsurvreg object containing the survival curves (with times, estimates and interval limits)
  # file = a string with the full path to the file name to be saved
  # sheet = a string with the name of the sheet to be created
  
  if(is.element("xlsx",installed.packages()[,1])==F){install.packages("xlsx")}
  if(is.element("flexsurv",installed.packages()[,1])==F){install.packages("flexsurv")}
  if(is.element("tools",installed.packages()[,1])==F){install.packages("tools")}
  
  # If it already exists, we need to append the data to a different sheet
  if (file.exists(file)) {
    wb <- loadWorkbook(file)
    # If worksheet already exists needs to replace it & overwrite it
    if (sheet %in% names(getSheets(wb))) {removeSheet(wb,sheetName=sheet)}
    sheet <- createSheet(wb,sheet)
    nobjs <- length(names(object))
    sc <- seq(from=1,by=5,length.out=nobjs)
    for (i in 1:nobjs) {
      addDataFrame(object[[i]],sheet=sheet,startRow=1,startColumn=sc[i],row.names=F)
    }
    saveWorkbook(wb,file)
  }
  
  # But if file does not exist, then create it
  if (!file.exists(file)) {
    exts <- file_ext(file)
    ## Should put some restriction as to what file extensions we want?
    wb <- createWorkbook(type=exts)
    sheet <- createSheet(wb,sheet)
    nobjs <- length(names(object))
    sc <- seq(from=1,by=5,length.out=nobjs)
    for (i in 1:nobjs) {
      addDataFrame(object[[i]],sheet=sheet,startRow=1,startColumn=sc[i],row.names=F)
    }
    saveWorkbook(wb,file)
  }
}


digitise <- function(surv_inp,nrisk_inp,km_output="KMdata.txt",ipd_output="IPDdata.txt") {
  # Post-process the data obtained by DigitizeIT to obtain the KM data and the individual level data
  # surv_inp = a txt file obtained by DigitizeIT and containing the input survival times from graph reading
  # nrisk_inp = a txt file obtained by DigitizeIT and containing the reported number at risk
  # km_output = the name of the file to which the KM data will be written
  # ipd_output = the name of the file to which the individual level data data will be written
  # Adapted from Patricia Guyot (2012)
  
  # Defines the working directory (same as the one where the DigitizeIT data are)
  working.dir <- dirname(surv_inp)
  ####  setwd(working.dir); working.dir <- paste0(getwd(),"/")
  tot.events<-"NA"  #tot.events = total no. of events reported. If not reported, then tot.events="NA"
  arm.id<-1         #arm indicator
  
  #Read in survival times read by digizeit
  digizeit <- read.table(surv_inp,header=TRUE,row.names=NULL)
  t.S<-digizeit[,2]     # times recorded from DigitizeIT
  S<-digizeit[,3]       # survival from DigitizeIT
  
  #Read in published numbers at risk, n.risk, at time, t.risk, lower and upper indexes for time interval
  pub.risk<-read.table(nrisk_inp,header=TRUE,row.names=NULL)
  ## Needs to get rid of possible time intervals with no digitised observations
  pub.risk <- pub.risk[pub.risk[,4]>0,]
  ## Needs to recode the first ever occurrence to 1??
  if (!(pub.risk[1,3]==1)) {pub.risk[1,3] <- 1}
  
  # Defines the variables needed for the algorithm
  t.risk<-pub.risk[,2]
  lower<-pub.risk[,3]
  upper<-pub.risk[,4]
  n.risk<-pub.risk[,5]
  n.int<-length(n.risk)
  n.t<- upper[n.int]
  
  #Initialise vectors
  arm <- rep(arm.id,n.risk[1])
  n.censor <- rep(0,(n.int-1))
  n.hat <- rep(n.risk[1]+1,n.t)
  cen <- d <- rep(0,n.t)
  KM.hat <- rep(1,n.t)
  last.i <- rep(1,n.int)
  sumdL <- 0
  
  # Executes Patricia's algorithm to determine censoring
  if (n.int > 1){
    #Time intervals 1,...,(n.int-1)
    for (i in 1:(n.int-1)){
      #First approximation of no. censored on interval i
      n.censor[i]<- round(n.risk[i]*S[lower[i+1]]/S[lower[i]]- n.risk[i+1])
      #Adjust tot. no. censored until n.hat = n.risk at start of interval (i+1)
      while((n.hat[lower[i+1]]>n.risk[i+1])||((n.hat[lower[i+1]]<n.risk[i+1])&&(n.censor[i]>0))){
        if (n.censor[i]<=0){
          cen[lower[i]:upper[i]]<-0
          n.censor[i]<-0
        }
        if (n.censor[i]>0){
          cen.t<-rep(0,n.censor[i])
          for (j in 1:n.censor[i]){
            cen.t[j]<- t.S[lower[i]] +
              j*(t.S[lower[(i+1)]]-t.S[lower[i]])/(n.censor[i]+1)
          }
          #Distribute censored observations evenly over time. Find no. censored on each time interval.
          cen[lower[i]:upper[i]]<-hist(cen.t,breaks=t.S[lower[i]:lower[(i+1)]],plot=F)$counts
        }
        #Find no. events and no. at risk on each interval to agree with K-M estimates read from curves
        n.hat[lower[i]]<-n.risk[i]
        last<-last.i[i]
        for (k in lower[i]:upper[i]){
          if (i==1 & k==lower[i]){
            d[k]<-0
            KM.hat[k]<-1
          }
          else {
            d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
            KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          }
          n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
          if (d[k] != 0) last<-k
        }
        n.censor[i]<- n.censor[i]+(n.hat[lower[i+1]]-n.risk[i+1])
      }
      if (n.hat[lower[i+1]]<n.risk[i+1]) n.risk[i+1]<-n.hat[lower[i+1]]
      last.i[(i+1)]<-last
    }
  }
  #Time interval n.int.
  if (n.int>1){
    #Assume same censor rate as average over previous time intervals.
    n.censor[n.int]<- min(round(sum(n.censor[1:(n.int-1)])*(t.S[upper[n.int]]-
                                                              t.S[lower[n.int]])/(t.S[upper[(n.int-1)]]-t.S[lower[1]])), n.risk[n.int])
  }
  if (n.int==1){n.censor[n.int]<-0}
  if (n.censor[n.int] <= 0){
    cen[lower[n.int]:(upper[n.int]-1)]<-0
    n.censor[n.int]<-0
  }
  if (n.censor[n.int]>0){
    cen.t<-rep(0,n.censor[n.int])
    for (j in 1:n.censor[n.int]){
      cen.t[j]<- t.S[lower[n.int]] +
        j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
    }
    cen[lower[n.int]:(upper[n.int]-1)]<-hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],plot=F)$counts
  }
  #Find no. events and no. at risk on each interval to agree with K-M estimates read from curves
  n.hat[lower[n.int]]<-n.risk[n.int]
  last<-last.i[n.int]
  for (k in lower[n.int]:upper[n.int]){
    if(KM.hat[last] !=0){
      d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))} else {d[k]<-0}
    KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
    n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
    #No. at risk cannot be negative
    if (n.hat[k+1] < 0) {
      n.hat[k+1]<-0
      cen[k]<-n.hat[k] - d[k]
    }
    if (d[k] != 0) last<-k
  }
  #If total no. of events reported, adjust no. censored so that total no. of events agrees.
  if (tot.events != "NA"){
    if (n.int>1){
      sumdL<-sum(d[1:upper[(n.int-1)]])
      #If total no. events already too big, then set events and censoring = 0 on all further time intervals
      if (sumdL >= tot.events){
        d[lower[n.int]:upper[n.int]]<- rep(0,(upper[n.int]-lower[n.int]+1))
        cen[lower[n.int]:(upper[n.int]-1)]<- rep(0,(upper[n.int]-lower[n.int]))
        n.hat[(lower[n.int]+1):(upper[n.int]+1)]<- rep(n.risk[n.int],(upper[n.int]+1-lower[n.int]))
      }
    }
    #Otherwise adjust no. censored to give correct total no. events
    if ((sumdL < tot.events)|| (n.int==1)){
      sumd<-sum(d[1:upper[n.int]])
      while ((sumd > tot.events)||((sumd< tot.events)&&(n.censor[n.int]>0))){
        n.censor[n.int]<- n.censor[n.int] + (sumd - tot.events)
        if (n.censor[n.int]<=0){
          cen[lower[n.int]:(upper[n.int]-1)]<-0
          n.censor[n.int]<-0
        }
        if (n.censor[n.int]>0){
          cen.t<-rep(0,n.censor[n.int])
          for (j in 1:n.censor[n.int]){
            cen.t[j]<- t.S[lower[n.int]] +
              j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
          }
          cen[lower[n.int]:(upper[n.int]-1)]<-hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],plot=F)$counts
        }
        n.hat[lower[n.int]]<-n.risk[n.int]
        last<-last.i[n.int]
        for (k in lower[n.int]:upper[n.int]){
          d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
          KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          if (k != upper[n.int]){
            n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
            #No. at risk cannot be negative
            if (n.hat[k+1] < 0) {
              n.hat[k+1]<-0
              cen[k]<-n.hat[k] - d[k]
            }
          }
          if (d[k] != 0) last<-k
        }
        sumd<- sum(d[1:upper[n.int]])
      }
    }
  }
  
  # Now writes the results to the output files
  KMdata <- data.frame(time=t.S,n.risk=n.hat[1:n.t],n.event=d,n.censored=cen)
  write.table(KMdata,km_output,sep="\t",row.names=FALSE,col.names=TRUE)
  
  # And forms IPD data
  #Initialise vectors
  t.IPD <- rep(t.S[n.t],n.risk[1])
  event.IPD <- rep(0,n.risk[1])
  #Write event time and event indicator (=1) for each event, as separate row in t.IPD and event.IPD
  k=1
  for (j in 1:n.t){
    if(d[j]!=0){
      t.IPD[k:(k+d[j]-1)]<- rep(t.S[j],d[j])
      event.IPD[k:(k+d[j]-1)]<- rep(1,d[j])
      k<-k+d[j]
    }
  }
  #Write censor time and event indicator (=0) for each censor, as separate row in t.IPD and event.IPD
  for (j in 1:(n.t-1)){
    if(cen[j]!=0){
      t.IPD[k:(k+cen[j]-1)]<- rep(((t.S[j]+t.S[j+1])/2),cen[j])
      event.IPD[k:(k+cen[j]-1)]<- rep(0,cen[j])
      k<-k+cen[j]
    }
  }
  #Output IPD
  IPD <- data.frame(time=t.IPD,event=event.IPD,arm)
  write.table(IPD,ipd_output,sep="\t",row.names=FALSE,col.names=TRUE)
  
  if (dirname(km_output)==".") {
    cat("\n")
    cat(paste0("Kaplan Meier data written to file: ",working.dir,km_output))    
  } else {
    cat("\n")
    cat(paste0("Kaplan Meier data written to file: ",km_output))    
  }
  if (dirname(ipd_output)==".") {
    cat("\n")
    cat(paste0("IPD data written to file: ",working.dir,ipd_output))    
  } else {
    cat("\n")
    cat(paste0("IPD data written to file: ",ipd_output))
  }
}



make.ipd <- function(ipd_files,ctr=1,var.labs=c("time","event","arm")) {
  ## Piles in the simulated IPD resulting from running digitise for more than one treatment arm  
  ## ipd_files = a list including the names of the IPD files created as output of digitise
  ## ctr = the index of the file associated with the control arm (default, the first file).
  ##       This will be coded as 0
  ## var.labs = a vector of labels for the column of the resulting data matrix. NB these
  ##            should match the arguments to the formula specified for fit.models. The
  ##            user can specify values. These should be 4 elements (ID, TIME, EVENT, ARM)
  
  # Identifies the number of arms (= number of IPD files)
  n_arms <- length(ipd_files)
  index <- 1:n_arms
  active <- index[-ctr]
  data <- read.table(ipd_files[[ctr]],header=TRUE,row.names=NULL)
  data[,"arm"] <- 0 # sets the value of "arm" to 0, for the control group
  arm.ind <- 1
  for (i in active) {
    tmp <- read.table(ipd_files[[index[i]]],header=TRUE,row.names=NULL)
    tmp[,"arm"] <- arm.ind
    data <- rbind(data,tmp)
    arm.ind <- arm.ind+1
  }
  colnames(data) <- var.labs
  return(data)
}

