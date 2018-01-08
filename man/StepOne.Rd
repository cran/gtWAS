\name{StepOne}
\alias{StepOne}
\title{
Compute minimum p value and information criteria statistics in one step
}
\description{
Compute minimum p value and information criteria statistics in one step by adding or removing a variable
}
\usage{
StepOne(findIn, independent, criteria = "SL", varIn, TMdata, sigma = 0)
}
\arguments{
  \item{findIn}{
Logical value for adding or removing independent variables in regression model, the parameter is ture for removing a variable otherwise adding a variable
}
  \item{independent}{
Indicator of independent variable to be used in linear model. 'B' is molecular marker effect, 'E' is expression effect, 'B(E)' is expression effect nested within molecular marker effect and 'E(B)' is molecular marker effect nested within expression effect
}
  \item{criteria}{
Specifies the criterion that uses to determine the order in which effects enter and/or leave at each step of the specified selection method including Akaike Information Criterion(AIC), the Corrected form of Akaike Information Criterion(AICc),Bayesian Information Criterion(BIC),Schwarz criterion(SBC),Hannan and Quinn Information Criterion(HQ), Significant Levels(SL) and so on
}
  \item{varIn}{
Sequence of vector for every independent variables, 1 indicates this independent variable stays in the regression model, and 0 is not in the model
}
  \item{TMdata}{
Phenotype data
}
  \item{sigma}{
The estimation of pure error variance from the full model in regression
}
}
\value{
A list of minimum p value or information criteria statistics, sequence id of independent variable staying in the model, linear model regression and rank of last step linear model
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
TMdata <- cbind(Tdata,alldata[,1:100])
findIn = FALSE
independent = "B"
varIn <- rep(0,100)
StepOne(findIn,independent,criteria="SBC",varIn,TMdata,sigma=0)
}
\keyword{ Model Selection }

