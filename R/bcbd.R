#' Create a Bayesian CBD model
#'
#' Initialise a BayMoMo object representing a Bayesian
#' CBD model
#'
#'
bcbd <- function(){
  #Construct the model
  BUGSmodel <- '
    model {
    #Observational model
    for (x in 1:k) {
      for (t in 1:n) {
        Dxt[x, t] ~ dpois(lambda.hat[x, t])
        logmu.hat[x, t] <- (kt[1, t] + (x-xbar) * kt[2, t])
        lambda.hat[x, t] <- Ext[x,t] * exp(logmu.hat[x, t])
      }
    }

    #Time series model
    for (t in 2:(n+h)) {
      kt[1:2, t] ~ dmnorm(d[] + kt[1:2, t - 1], sigma.inv.kt[,])
    }

    #Predictions
    for (x in 1:k){
      for (t in (n+1):(n+h)){
        logmu.hat[x, t] <- (kt[1, t] + (x-xbar)*kt[2,t])
      }
    }

    #Priors
    #kt priors
    kt[1, 1] ~ dnorm(k10[1], 1/sigma.k1[1]^2)
    kt[2, 1] ~ dnorm(k10[2], 1/sigma.k1[2]^2)
    d[1] ~ dnorm(d0[1], 1/sigma.d[1]^2)
    d[2] ~ dnorm(d0[2], 1/sigma.d[2]^2)
    sigma.inv.kt[1:2, 1:2] ~ dwish(0.5*V[,], 2)
    }
    '

  out <- list(BUGSmodel = BUGSmodel,
              model = cbd(link = "log"))
  class(out) <- "BayMoMo"
  out
}
