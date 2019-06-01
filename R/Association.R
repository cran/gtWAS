Association <- function (Tdata,alldata,independent="B(E)",Elevels=c(0.05,0.95),selection='stepwise',select="SL",Choose="NULL",SL=c(0.05,0.15,0.15),correct="Bonferroni"){
  nT <- ncol(Tdata)
  if(independent == "E(B)" | independent == "B(E)") {
    nbase <- ncol(alldata)/2
  }else if (independent == "E" | independent == "B") {
    nbase <- ncol(alldata)
  }else {
    stop("independent must be B or E or B(E) or E(B)")
  }
  if (any(SL <= 0 | SL >= 1)) {
    stop("significant levels must be 0~1")
  }
  pvalue <- matrix(NA, nbase, nT * 2)
  pvalue <- as.data.frame(pvalue)
  colnames(pvalue) <- sort(c(paste(colnames(Tdata), "_Pvalue", sep = ""), paste(colnames(Tdata), "_Significant", sep = "")))
  vecT <- colnames(Tdata)
  if (independent == "E") {
    rownames(pvalue) <- paste("exp", c(1:nbase), sep = "")
  } else if (independent == "B") {
    rownames(pvalue) <- paste("base", c(1:nbase), sep = "")
  } else if (independent == "E(B)") {
    rownames(pvalue) <- paste("exp",c(1:nbase),":base",c(1:nbase),sep = "")
  } else if (independent == "B(E)") {
    rownames(pvalue) <- paste("base",c(1:nbase),":exp",c(1:nbase),sep = "")
  }
  if(independent != "B" && !is.null(Elevels)){
    if(independent == "E"){
      En <- 0
    }else{
      En <- nbase
    }
    Edata <- alldata[,En+1:nbase]
    nEL <- length(Elevels)
    Elevels <- sort(Elevels)
    obs <- nrow(alldata)
    qua <- sapply(Edata, quantile, probs=Elevels, na.rm=TRUE)
    
    nlist <- list()
    Lfactor <- factor(paste("L",1:(nEL+1),sep = ""))
    index <- nEL > 1
    for(n in 1:nbase){
      nlist[1] <- list(which(alldata[,En+n] <= qua[1,n]))
      nlist[nEL+1] <- list(which(alldata[,En+n] > qua[nEL,n]))
      if(index){
        for(i in 1:(nEL-1)){
          nlist[i+1] <- list(which(alldata[,En+n] <= qua[i+1,n] & alldata[,En+n] > qua[i,n]))
        }
      }
      for(k in 1:(nEL+1)){
        if(length(nlist[[k]])>0){
          alldata[nlist[[k]],En+n] <- Lfactor[k]
        }
      }
    }
  }
  #numeric base and expression data
  alldata[,1:nbase] <- sapply(alldata[,1:nbase], as.numeric)

  result <- list()
  for (n in 1:nT) {
    TMdata <- cbind(Tdata[colnames(Tdata)[n]], alldata)
    CF_RCpp <- stp(TMdata, independent, selection, select,SL[2],SL[3],Choose)
    CF <- CF_RCpp$variate[-1]
    if (length(CF) == 0) {
      Mx <- 1
      cat("there is no cofactors in stepwise regression for",colnames(Tdata)[n],"\n")
    } else {
      Mx <- paste0(CF, collapse = " + ")
      cat("significant independent variables in stepwise regression for",colnames(Tdata)[n],"\n",CF,"\n")
    }
    Myx <- paste(vecT[n], " ~ ", Mx, sep = "")
    indvar <- rownames(pvalue)
    for (i in indvar) {
      if (i %in% CF) {
        fMyx <- Myx
        Mnf <- CF
        Mnr <- Mnf[!Mnf %in% i]
        if (length(Mnr) == 0) {
          rMx <- 1
        }
        else {
          rMx <- paste0(Mnr, collapse = " + ")
        }
        rMyx <- paste(vecT[n], " ~ ", rMx, sep = "")
      }
      else {
        rMyx <- Myx
        fMyx <- paste(rMyx, " + ", i, sep = "")
      }
      lm.res <- lm(as.formula(rMyx), data = as.data.frame(TMdata))
      lm.ful <- lm(as.formula(fMyx), data = as.data.frame(TMdata))
      anovar <- anova(lm.res, lm.ful)
      pvalue[i, n * 2 - 1] <- anovar$`Pr(>F)`[2]
    }
    if (correct == "P") {
      pvalue[pvalue[, 2 * n - 1] < SL[1], n * 2] <- "* P"
      BaseP <- pvalue[pvalue[, n * 2] %in% "* P", (n * 2 - 1):(n * 2)]
      if (nrow(BaseP) != 0) {
        result[1 + n] <- list(BaseP)
      } else {
        result[1 + n] <- list("No significant variable at P levels")
      }
    }else if (correct == "Bonferroni") {
      pvalue[pvalue[, 2 * n - 1] < SL[1]/nbase, n * 2] <- "* Bonferroni"
      BaseBonferroni <- pvalue[pvalue[, n * 2] %in% "* Bonferroni",(n * 2 - 1):(n * 2)]
      if (nrow(BaseBonferroni) != 0) {
        result[1 + n] <- list(BaseBonferroni)
      } else {
        result[1 + n] <- list("No significant variable at Bonferroni levels")
      }
    } else {
      stop("correct must be p or Bonferroni")
    }
    pvalue[pvalue[, n * 2 - 1] < SL[1], n * 2] <- "* P"
    pvalue[pvalue[, n * 2 - 1] < SL[1]/nbase, n * 2] <- "* Bonferroni"
  }
  result[1] <- list(pvalue)
  names(result) <- c("Pvalue",colnames(Tdata))
  return(result)
}