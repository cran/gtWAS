stepwise <- function (AllData, independent, selection = "stepwise", select = "SL", sle = 0.05, sls = 0.05, Choose = "SBC"){
	if(independent == "E(B)"){
		nX <- (ncol(AllData)-1)/2
		Basedat <- AllData[,c(2:(nX+1))]
		#xf_ncl <- sum(sapply(Basedat,nlevels))
		vecName <- paste("exp",c(1:nX),":","base",c(1:nX),sep="")
		Mxf <- paste0(paste("exp",c(1:nX),":base",c(1:nX),sep=""),collapse=" + ")
	}else if(independent == "B(E)"){
		nX <- (ncol(AllData)-1)/2
		Basedat <- AllData[,c(2:(nX+1))]
		#xf_ncl <- sum(sapply(Basedat,nlevels))
		vecName <- paste("base",c(1:nX),":","exp",c(1:nX),sep="")
		Mxf <- paste0(paste("exp",c(1:nX),":base",c(1:nX),sep=""),collapse=" + ")
	}else if(independent == "E"){
		nX <- ncol(AllData)-1
		vecName <- paste("exp",c(1:nX),sep="")
		Mxf <- paste0(paste("exp",c(1:nX),sep=""),collapse=" + ")
	}else if(independent == "B"){
		nX <- ncol(AllData)-1
		vecName <- paste("base",c(1:nX),sep="")
		Mxf <- paste0(paste("base",c(1:nX),sep=""),collapse=" + ")
	}else{
		stop("independent must be B or E or B(E) or E(B)")
	}
	nObs <- nrow(AllData)
	vecT <- colnames(AllData)[1]
	Myxf <- paste(vecT," ~ ",Mxf,sep="")
  if (select == "BIC" | Choose == "BIC") {
    if (nX > nObs) {
        sigma <- 0
    }else {
        lmf <- lm(as.formula(Myxf),data=AllData)
        sigma <- sum(deviance(lmf)/df.residual(lmf))
    }
  }else {
    sigma <- 0
  }
  slcOpt <- list(serial = 'numeric', bestValue = 'numeric', bestPoint = 'numeric', enOrRe = 'logical', nVarIn = 'numeric')
  class(slcOpt) <- select
  chsOpt <- list(serial = 'numeric', bestValue = 'numeric', bestPoint = 'numeric', enOrRe = 'logical', nVarIn = 'numeric')
  class(chsOpt) <- Choose
  slcOpt$serial <- 1
  slcOpt$bestPoint <- 0
  slcOpt$enOrRe <- TRUE
  slcOpt$nVarIn <- 1
	
	Myb <- paste(vecT," ~ 1",sep="")
	lmb <- lm(as.formula(Myb),data=AllData)
    if (class(slcOpt) == "SL") {
        slcOpt$bestValue <- 1
    }else {
        slcOpt$bestValue <- InfoCriteria(class(slcOpt),lmb,nObs,sigma)
    }
    chsOpt$bestValue <- InfoCriteria(class(chsOpt),lmb,nObs,sigma)
	
    addVar <- TRUE
    varIn <- rep(0, nX)
    while (TRUE) {
        findIn <- if (addVar == TRUE) FALSE else TRUE
        pointer <- if (addVar == TRUE) 1 else -1
        p <- slcOpt$nVarIn[slcOpt$serial]
        stepvalue <- StepOne(findIn, independent, class(slcOpt), varIn, AllData, sigma)
        if (stepvalue$rank0 == stepvalue$bestlm$rank && findIn == FALSE) {
            break
        }else {
			if (class(slcOpt) == 'SL') {
				if (findIn == TRUE) {
					indicator <- stepvalue$PvorCr > sls
				} else {
					indicator <- stepvalue$PvorCr < sle
				}
			} else{
				indicator <- round(stepvalue$PvorCr,digits=7) <= round(slcOpt$bestValue[slcOpt$serial],digits=7)
			}
            if(indicator == TRUE){			
                if (addVar == TRUE) {
                  Order <- stepvalue$pointer
                }else {
                  stay <- which(varIn %in% 1)
                  Order <- stay[stepvalue$pointer]
                }
                slcOpt$serial <- slcOpt$serial + 1
                slcOpt$bestPoint[slcOpt$serial] <- Order
                slcOpt$bestValue[slcOpt$serial] <- stepvalue$PvorCr
                slcOpt$enOrRe[slcOpt$serial] <- addVar
                slcOpt$nVarIn[slcOpt$serial] <- if (addVar == TRUE) p + 1 else p - 1
                chsOpt$bestValue[slcOpt$serial] <- InfoCriteria(class(chsOpt),stepvalue$bestlm,nObs,sigma)
                varIn[Order] <- varIn[Order] + pointer
                if (selection == "forward") {
                  next
                }else if (selection == "stepwise") {
                  if (addVar == FALSE) {
                    next
                  }else if (addVar == TRUE) {
                    addVar <- FALSE
                    next
                  }
               }
            }else {
                if (selection == "stepwise" && addVar == FALSE) {
                  addVar <- TRUE
                  next
                }else {
                  break
                }
            }
        }
    }
    varName <- array(FALSE, slcOpt$serial)
    varName[1:slcOpt$serial] <- c("intercept", vecName[slcOpt$bestPoint[-1]])
    results <- data.frame(step = 0:(slcOpt$serial - 1), name = varName, 
        enter = slcOpt$enOrRe, best = slcOpt$bestPoint, nX = slcOpt$nVarIn, 
        slcv = slcOpt$bestValue, chsv = chsOpt$bestValue)
    return(results)
}
