InfoCriteria <-
function(criteria,lmresult,nObs,sigma_sqr){
  p <- lmresult$rank
  nY <- 1
  lm_res <- residuals(lmresult)
  SSEp <- t(lm_res) %*% lm_res      #SSEp: Sum Squared error with p X included
  #SSEpSq <- t(SSEp)%*%SSEp         #SSEpSq: Sum Squared error with p X included
  #SSESqdet <- abs(det(SSEpSq))
  SSEdet <- abs(det(SSEp))
  
  if(criteria=="AIC"){
    Fitstat <- nObs*log(SSEdet/nObs)+(2*p*nY*nObs+nY*(nY+1))/nObs-2/nObs+nObs+2         # documents in SAS9.3
  }
  if(criteria=="AICc"){
    Fitstat <- nObs*log(SSEdet/nObs)+nObs*(nObs+p)*nY/(nObs-p-nY-1)
  }
  if(criteria=="BIC"){
    Fitstat <- nObs*log(SSEdet/nObs)+2*(2+p)*(nObs*sigma_sqr/SSEdet)-2*(nObs*sigma_sqr/SSEdet)^2   ## maybe some error in 
  }
  if(criteria=="HQ"){
    Fitstat <- log(SSEdet)+2*log(log(nObs))*p*nY/nObs
  }
  #if(criteria=="HQc"){
  #  Fitstat <- log(SSESqdet)+2*log(log(nObs))*p*nY/(nObs-p-nY-1)
  #}
  if(criteria=="SBC"){
    Fitstat <- nObs*log(SSEdet/nObs)+log(nObs)*p
  }
  return(Fitstat)
}
