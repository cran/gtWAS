\name{stepwise}
\alias{stepwise}
\title{
stepwise regression
}
\description{
Stepwise regression for model selection using linear model
}
\usage{
stepwise(AllData, independent, selection = "stepwise", select = "SL",
sle = 0.05, sls = 0.05, Choose = "SBC")
}
\arguments{
  \item{AllData}{
Data about dependent and independent variable data
}
  \item{independent}{
Indicator of independent variable to be used in linear model. 'B' is molecular marker effect, 'E' is expression data, 'B(E)' is expression effect nested within molecular marker effect and 'E(B)' is molecular marker effect nested within expression effect
}
  \item{selection}{
Model selection method including "forward" and "stepwise",forward selection starts with no effects in the model and adds effects, while stepwise regression is similar to the forward method except that effects already in the model do not necessarily stay there
}
  \item{select}{
Specifies the criterion that uses to determine the order in which effects enter and/or leave at each step of the specified selection method including Akaike Information Criterion(AIC), the Corrected form of Akaike Information Criterion(AICc),Bayesian Information Criterion(BIC),Schwarz criterion(SBC),Hannan and Quinn Information Criterion(HQ), Significant Levels(SL) and so on
}
  \item{sle}{
Specifies the significance level for entry
}
  \item{sls}{
Specifies the significance level for staying in the model
}
  \item{Choose}{
Chooses from the list of models at the steps of the selection process the model that yields the best value of the specified criterion. If the optimal value of the specified criterion occurs for models at more than one step, then the model with the smallest number of parameters is chosen. If you do not specify the Choose option, then the model selected is the model at the final step in the selection process
}
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
data(Tdata)
data(alldata)
independent <- "B"
nbase <- 100
AllData <- cbind(Tdata[colnames(Tdata)[1]],alldata[,1:nbase])
AllData <- sapply(AllData, as.numeric)
AllData <- as.data.frame(AllData)
stepwise(AllData,independent,selection="stepwise",select="SBC",sle=0.05,sls=0.05,Choose="SBC")
}
\keyword{ Model Selection }

