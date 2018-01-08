gtWAS <- function(Tdata,alldata,independent,selection='stepwise',select="SL",Choose="SBC",vecThr=c(0.001,0.001,0.05),
                  correct="Bonferroni"){
  ## initialize p values
  nT <- ncol(Tdata)
  if(independent=="E(B)" | independent=="B(E)"){
    nbase <- ncol(alldata)/2
  }else if(independent=="E" | independent=="B"){
    nbase <- ncol(alldata)
  }else{
    stop("independent must be B or E or B(E) or E(B)")
  }
  
  pvalue <- matrix(NA,nbase,nT*2)
  pvalue <- as.data.frame(pvalue)
  colnames(pvalue) <- sort(c(paste(colnames(Tdata),"_P",sep=""),paste(colnames(Tdata),"_Sig",sep="")))
  vecT <- colnames(Tdata)
  if(independent=="E"){  # just for expression data
    rownames(pvalue) <- paste("exp",c(1:nbase),sep="")
  }else if(independent=="B"){ # just for base data
    rownames(pvalue) <- paste("base",c(1:nbase),sep="")
  }else if(independent=="E(B)"){ # just for expression and base data
    rownames(pvalue) <- paste("exp",c(1:nbase),":base",c(1:nbase),sep="")
  }else if(independent=="B(E)"){
    rownames(pvalue) <- paste("base",c(1:nbase),":exp",c(1:nbase),sep="")
  }else{
    stop("independent must be B or E or B(E) or E(B)")
  }
  
  ## compute all independent variables' p value based on expression or base or both data
  result <- list()
  if(!all(vecThr > 0 && vecThr < 1)){
    stop("significant levels must be 0~1")
  }
  sle=vecThr[1]
  sls=vecThr[2]

  #n = 1
  for(n in 1:nT){
    ## stepwise regression
    if(independent == "E"){
      TMdata <- cbind(Tdata[colnames(Tdata)[n]],alldata)
    }else if(independent == "B"){
      TMdata <- cbind(Tdata[colnames(Tdata)[n]],alldata)
      TMdata <- sapply(TMdata, as.numeric)
      TMdata <- as.data.frame(TMdata)
    }else if(independent == "E(B)" | independent == "B(E)"){
      TMdata <- cbind(Tdata[colnames(Tdata)[n]],alldata)
    }else{
      stop("independent must be B or E or B(E) or E(B)")
    }
    
    CF_RCpp <- stepwise(TMdata,independent,selection,select,sle,sls,Choose)
    Xrep <- CF_RCpp[which(duplicated(CF_RCpp[,"name"])),"name"]
    if(length(Xrep)!=0){
      Xomit <- NULL
      for(m in 1:length(Xrep)){
        #m=2
        if(sum(CF_RCpp[,"name"] %in% Xrep[m])%%2==0){
          Xomit <- append(Xomit,Xrep[m])
        }
      }
      tempCF_P1 <- unique(setdiff(CF_RCpp$name,Xomit))
    }else{
      tempCF_P1 <- CF_RCpp[,"name"]
    }
    CF <- tempCF_P1[-1]       #delete "intercept"
    if(length(CF)==0){
      Mx <- 1
      cat("There is no cofactors for ",colnames(Tdata)[n]," in stepwise regression\n")
    }else{
      Mx <- paste0(CF,collapse=" + ")
      cat(colnames(Tdata)[n],": significant independent variables in stepwise regression: ",CF,"\n")
    }
    Myx <- paste(vecT[n]," ~ ",Mx,sep="")
    indvar <- rownames(pvalue)
    #i=indvar[2]
    for(i in indvar){
      if(i %in% CF){
        fMyx <- Myx
        Mnf <- CF
        Mnr <- Mnf[!Mnf %in% i]
        if(length(Mnr)==0){
          rMx <- 1
        }else{
          rMx <- paste0(Mnr,collapse=" + ")
        }
        rMyx <- paste(vecT[n]," ~ ",rMx,sep="")
      }else{
        rMyx <- Myx
        fMyx <- paste(rMyx," + ",i,sep="")
      }
      lm.res <- lm(as.formula(rMyx),data=as.data.frame(TMdata))
      lm.ful <- lm(as.formula(fMyx),data=as.data.frame(TMdata))
      anovar <- anova(lm.res,lm.ful)
      pvalue[i,n*2-1] <- anovar$`Pr(>F)`[2]
    }
    
    if(correct == 'P'){
      pvalue[pvalue[,2*n-1] < vecThr[3],n*2] <- "* P"
      BaseP <- pvalue[pvalue[,n*2] %in% "* P",(n*2-1):(n*2)]
      if(nrow(BaseP) != 0){
        result[1+n] <- list(BaseP)
      }else{
        result[1+n] <- list("No significant variable at P levels")
      }
    }else if(correct == 'Bonferroni'){
      pvalue[pvalue[,2*n-1] < vecThr[3]/nbase,n*2] <- "* Bonferroni"
      BaseBonferroni <- pvalue[pvalue[,n*2] %in% "* Bonferroni",(n*2-1):(n*2)]
      if(nrow(BaseBonferroni) != 0){
        result[1+n] <- list(BaseBonferroni)
      }else{
        result[1+n] <- list("No significant variable at Bonferroni levels")
      }
    }else{
      stop("correct must be p or Bonferroni")
    }
    pvalue[pvalue[,n*2-1] < vecThr[3],n*2] <- "* P"
    pvalue[pvalue[,n*2-1] < vecThr[3]/nbase,n*2] <- "* Bonferroni"
  }
  result[1] <- list(pvalue)
  return(result)
}