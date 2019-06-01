\name{ModelFit}
\alias{ModelFit}
\title{
Compute model fit statistics 
}
\description{
Compute model fit statistics based on a given criteria for linear model function
}
\usage{
ModelFit(criteria, lmresult, nObs, sigma_sqr)
}
\arguments{
  \item{criteria}{
The class of criteria including Akaike information criterion(AIC), the corrected form of Akaike information criterion(AICc), Bayesian information criterion(BIC), Schwarz criterion(SBC) and significant levels(SL)
}
  \item{lmresult}{
Result of linear model function
}
  \item{nObs}{
Number of observation
}
  \item{sigma_sqr}{
The estimation of pure error variance for the full model in regression
}
}
\value{
A numeric of model fit statistics 
}
\references{
Hurvich, C. M., & Tsai, C. (1989). Regression and time series model selection in small samples. Biometrika, 76(2), 297-307.

Judge, & GeorgeG. (1985). The Theory and practice of econometrics /-2nd ed. The Theory and practice of econometrics /. Wiley.

Mcquarrie, A. D. R., & Tsai, C. L. (1998). Regression and Time Series Model Selection. Regression and time series model selection /. World Scientific.

R.S. Sparks, W. Zucchini, & D. Coutsourides. (1985). On variable selection in multivariate regression. Communication in Statistics- Theory and Methods, 14(7), 1569-1587.

Sawa, T. (1978). Information criteria for discriminating among alternative regression models. Econometrica, 46(6), 1273-1291.

Schwarz, G. (1978). Estimating the dimension of a model. Annals of Statistics, 6(2), pags. 15-18.
}
\author{
JunhuiLi
}

\examples{
set.seed(4)
YX <- matrix(rnorm(200,20,4),20,10)
YX <- as.data.frame(YX)
colnames(YX) <- c("Y1","Y2",paste("X",c(1:8),sep=""))
lm_formula <- as.formula("Y1~X1+X2+X3+X4+X5")
lmresult <- lm(lm_formula,data=YX)
ModelFit("SBC", lmresult, nrow(YX), 0)
}
\keyword{Model Selection}

