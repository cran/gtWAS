StepOne <- function(findIn,independent,criteria,varIn,TMdata,sigma){
	nObs <- nrow(TMdata)
	vecT <- colnames(TMdata)[1]
	stay <- which(varIn %in% 1)
	leave <- which(varIn %in% 0)
	pointer <- 1
	if(length(stay) == 0){
		Mx0 <- 1
	}else{
		if(independent=="E"){	# just for expression data
			Mx0 <- paste0(paste("exp",stay,sep=""),collapse=" + ")
		}else if(independent=="B"){ # just for base data
			Mx0 <- paste0(paste("base",stay,sep=""),collapse=" + ")
			}else if(independent=="E(B)"){ # for expression nested in base data
			Mx0 <- paste0(paste("exp",stay,":base",stay,sep=""),collapse=" + ")
		}else if(independent=="B(E)"){	# for base nested in expression data
		Mx0 <- paste0(paste("base",stay,":exp",stay,sep=""),collapse=" + ")
		}else{
			stop("independent must be B or E or B(E) or E(B)")
		}
	}
	Myx0 <- paste(vecT," ~ ",Mx0,sep="")
	lm0 <- lm(as.formula(Myx0),data=TMdata)
	bestlm <- lm0
	rank0 <- lm0$rank
	rank_return <- rank0
	if(criteria == "SL"){
		PvorCr0 <- 1 
	}else{
		PvorCr0 <- 1e+10
	}
	if(findIn == FALSE){
		if(length(leave) > 0){
			for(i in leave){
				if(independent=="E"){	# just for expression data
					Mx <- paste0(paste("exp",c(stay,i),sep=""),collapse=" + ")
				}else if(independent=="B"){ # just for base data
					Mx <- paste0(paste("base",c(stay,i),sep=""),collapse=" + ")
				}else if(independent=="E(B)"){ # just for expression and base data
					Mx <- paste0(paste("exp",c(stay,i),":base",c(stay,i),sep=""),collapse=" + ")
				}else if(independent=="B(E)"){
				Mx <- paste0(paste("base",c(stay,i),":exp",c(stay,i),sep=""),collapse=" + ")
				}else{
					stop("independent must be B or E or B(E) or E(B)")
				}
		Myx <- paste(vecT," ~ ",Mx,sep="")
		lmf <- lm(as.formula(Myx),data=TMdata)
		rank1 <- lmf$rank
		if(rank0 - rank1 != 0){					 ##otherwise: rank0=rank1, variables entered is invalid
			if(criteria == "SL"){
						PvorCr <- anova(lm0,lmf)$'Pr(>F)'[2]
						if(is.finite(PvorCr)){
							if(PvorCr < PvorCr0){
								PvorCr0 <- PvorCr
								pointer <- i
								bestlm <- lmf
								rank_return <- rank1
							}
						}#else{
							#print("Selection stopped because all candidate effects for entry are linearly dependent on effects in the model.\n")
						#}
					}else{
						PvorCr <- ModelFit(criteria,lmf,nObs,sigma)
						if(PvorCr < PvorCr0){
							PvorCr0 <- PvorCr
							pointer <- i
							bestlm <- lmf
							rank_return <- rank1
						}
					}
				}
			}#i
	}
	}else{
		if(length(stay) > 0){
			for(i in stay){
				if(length(stay) == 1){
					Mx <- 1
				}else{
					if(independent=="E"){	# just for expression data
							Mx <- paste0(paste("exp",stay[!stay %in% i],sep=""),collapse=" + ")
						}else if(independent=="B"){ # just for base data
							Mx <- paste0(paste("base",stay[!stay %in% i],sep=""),collapse=" + ")
						}else if(independent=="E(B)"){ # just for expression nested in base data
							Mx <- paste0(paste("exp",stay[!stay %in% i],":base",stay[!stay %in% i],sep=""),collapse=" + ")
						}else if(independent=="B(E)"){ # just for base nested in expression data
						Mx <- paste0(paste("base",stay[!stay %in% i],":exp",stay[!stay %in% i],sep=""),collapse=" + ")
					}else{
						stop("independent must be B or E or B(E) or E(B)")
					}
				}
				Myx <- paste(vecT," ~ ",Mx,sep="")
				lmr <- lm(as.formula(Myx),data=TMdata)
				rank1 <- lmr$rank
				if(rank0 - rank1 != 0){
					if(criteria=="SL"){
						PvorCr <- anova(lmr,lm0)$'Pr(>F)'[2]
						if(is.finite(PvorCr)){
							if(PvorCr < PvorCr0){
								PvorCr0 <- PvorCr
								pointer <- i
								bestlm <- lmr
								rank_return <- rank1
							}
						}#else{
								#print("Selection stopped because all candidate effects for entry are linearly dependent on effects in the model.\n")
							#}
					}else{
						PvorCr <- ModelFit(criteria,lmr,nObs,sigma)
						if(PvorCr < PvorCr0){
							PvorCr0 <- PvorCr
							pointer <- i
							bestlm <- lmr
							rank_return <- rank1
						}
					}
				}
			}
		}
	}#findIn
	revalue <- list()
	revalue$PvorCr <- PvorCr0
	revalue$pointer <- pointer
	revalue$bestlm <- bestlm
	revalue$rank0 <- rank0
	#revalue$rank <- rank_return
	return(revalue)
}
