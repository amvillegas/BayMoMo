#' Create a Bayesian Lee-Carter model
#'
#' Initialise a BayMoMo object representing a Bayesuab
#' Lee-Carter model
#'
#'
blc <- function(staticAgePrior = c("Gamma", "Normal", "RW", "RW2"),
                periodAgePrior = c("Normal", "RW", "RW2")){

  #Observation equation
  obsModel <- '
    #Observational model
    for (x in 1:k){
      for (t in 1:n){
        Dxt[x,t] ~ dpois(lambda.hat[x,t])
        logmu.hat[x,t] <- (ax[x] + bx[x]*kt[t])
        lambda.hat[x,t] <- Ext[x,t]*exp(logmu.hat[x,t])
      }
    }
    '
  tsModel <- '
    #Time series model
    kt[1] <- 0
    for (t in 2:(n+h)){
      kt[t] ~ dnorm(d + kt[t-1], 1/sigma.kt^2)
    }

    #Predictions
    for (x in 1:k){
      for (t in (n+1):(n+h)){
        logmu.hat[x,t] <- (ax[x] + bx[x]*kt[t])
      }
    }
    '
  #kappa prior
  kPrior <- '
    #kt priors
    d ~ dnorm(d0, 1/sigma.d^2)
    sigma.kt ~ dunif(0,A.sigma.kt)
    '

  #alpha prior
  staticAgePrior <- match.arg(staticAgePrior)
  aPrior <- switch(staticAgePrior,
         Gamma = '
    #ax prior
    for (x in 1:k){
      ax[x] <- log(ex[x])
      ex[x] ~ dgamma(a.ex[x],b.ex[x])
    }
    ',
         Normal = '
    #ax prior
    for (x in 1:k){
      ax[x] ~ dnorm(nu.ax,1/sigma.ax^2)
    }
    sigma.ax ~ dunif(0,A.sigma.ax)
   ',
         RW = '
    #ax prior
    ax[1] ~ dnorm(0, 1/sigma.ax1^2)
    for (x in 2:k){
      ax[x] ~ dnorm(ax[x-1],1/tau.ax^2)
    }
    tau.ax ~ dunif(0,A.tau.ax)
   ',
         RW2 = '
    #ax prior
    ax[1] ~ dnorm(0, 1/sigma.ax1^2)
    ax[2] ~ dnorm(0, 1/sigma.ax2^2)
    for (x in 3:k){
      ax[x] ~ dnorm(2*ax[x-1]-ax[x-2],1/tau.ax^2)
    }
    tau.ax ~ dunif(0,A.tau.ax)
   '
         )

  #beta prior
  periodAgePrior <- match.arg(periodAgePrior)
  bPrior <- switch(periodAgePrior,
         Normal = '
    #bx prior
    for (x in 1:(k-1)){
      bx[x] ~ dnorm(1/k,1/sigma.bx^2)
    }
    sigma.bx ~ dunif(0,A.sigma.bx)
    bx[k] <- 1- sum(bx[1:(k-1)])
   ',
         RW = '
    #bx prior
    bx[1] ~ dnorm(1/k, 1/sigma.bx1^2)
    for (x in 2:(k-1)){
      bx[x] ~ dnorm(bx[x-1],1/tau.bx^2)
    }
    tau.bx ~ dunif(0,A.tau.bx)
    bx[k] <- 1- sum(bx[1:(k-1)])
   ',
         RW2 = '
    #bx prior
    bx[1] ~ dnorm(1/k, 1/sigma.bx1^2)
    bx[2] ~ dnorm(1/k, 1/sigma.bx2^2)
    for (x in 3:(k-1)){
      bx[x] ~ dnorm(2*bx[x-1]-bx[x-2],1/tau.bx^2)
    }
    tau.bx ~ dunif(0,A.tau.bx)
    bx[k] <- 1- sum(bx[1:(k-1)])
   '
         )

  #Construct the model
  BUGSmodel <- paste("model {", obsModel, tsModel, kPrior, aPrior, bPrior, "}", sep = "\n")

  out <- list(staticAgePrior = staticAgePrior,
              periodAgePrior = periodAgePrior,
              BUGSmodel = BUGSmodel,
              model = lc(const = "first"))
  class(out) <- "BayMoMo"
  out
}
