\name{Association}
\alias{Association}
\title{
Genome and Transcriptome Wide Association
}
\description{
Reveal the association relationship between phenotype and molecular marker, expression effect, expression effect nested within molecular marker and molecular marker effect nested within expression effect
}
\usage{
Association(Tdata,alldata,independent="B(E)",Elevels=c(0.05,0.95),selection="stepwise",
      select="SL",Choose=NULL,SL=c(0.05,0.15,0.15),correct="Bonferroni")
}
\arguments{
  \item{Tdata}{
Phenotye data
}
  \item{alldata}{
Independent variables including molecular marker or corresponding expression effect related to marker on transcriptome level
}
  \item{independent}{
Indicator of independent variable to be used in linear model. 'B' is molecular marker effect, 'E' is expression effect, 'B(E)' is molecular marker nesting expression effect and 'E(B)' is expression effect nesting molecular marker effect
}
  \item{Elevels}{
Percentage of threshold value for different expression levels
}
\item{selection}{
Model selection method including "forward" and "stepwise",forward selection starts with no effects in the model and adds effects, while stepwise regression is similar to the forward method except that effects already in the model do not necessarily stay there
}
  \item{select}{
Specifies the criterion that uses to determine the order in which effects enter and/or leave at each step of the specified selection method including Akaike Information Criterion(AIC), the Corrected form of Akaike Information Criterion(AICc),Bayesian Information Criterion(BIC),Schwarz criterion(SBC), Significant Levels(SL) and so on
}
  \item{Choose}{
Chooses from the list of models at the steps of the selection process the model that yields the best value of the specified criterion. If the optimal value of the specified criterion occurs for models at more than one step, then the model with the smallest number of parameters is chosen. If you do not specify the Choose option, then the model selected is the model at the final step in the selection process
}
  \item{SL}{
Thresholds for significant levels of association and stepwise regression
}
  \item{correct}{
Bonferroni correct or the p value method for significant levels, default is bonferroni
}
}

\value{
p value of all effect and significant ones
}
\references{
Junhui Li, Haixiao Hu, Yujie Meng, Kun Cheng, Guoliang Li, Wenxin Liu, and Shaojiang Chen.(2016)Pleiotropic QTL detection for stalk traits in maize and related R package programming. Journal of China Agricultural University. DOI 10.11841/j.issn.1007-4333.2016.06.00

Gusev, A., Ko, A., Shi, H., Bhatia, G., Chung, W., & Penninx, B. W., et al. (2016). Integrative approaches for large-scale transcriptome-wide association studies. Nature Genetics, 48(3), 245.

Hurvich, C. M., & Tsai, C. (1989). Regression and time series model selection in small samples. Biometrika, 76(2), 297-307.

Judge, & GeorgeG. (1985). The Theory and practice of econometrics /-2nd ed. The Theory and practice of econometrics /. Wiley.

Mcquarrie, A. D. R., & Tsai, C. L. (1998). Regression and Time Series Model Selection. Regression and time series model selection /. World Scientific.

Sawa, T. (1978). Information criteria for discriminating among alternative regression models. Econometrica, 46(6), 1273-1291.

Schwarz, G. (1978). Estimating the dimension of a model. Annals of Statistics, 6(2), pags. 15-18.
}
\author{
JunhuiLi
}
\examples{
data(Tdata)
data(alldata)
Edata <- alldata[,1:100+100]
Bdata <- alldata[,1:100]
BE <- "B(E)"
EB <- "E(B)"
B <- "B"
E <- "E"

#for "B(E)"
#Association(Tdata,alldata,BE,Elevels=c(0.05,0.95),selection='stepwise',
#select="SL",Choose=NULL,SL=c(0.05,0.15,0.15),correct = "Bonferroni")

#for "E(B)" with Elevels = null
#Association(Tdata,alldata,EB,Elevels=NULL,selection='stepwise',
#select="SL",Choose=NULL,SL=c(0.05,0.15,0.15),correct = "Bonferroni")

#for "E" with Elevels = null
#Association(Tdata,Edata,E,Elevels=NULL,selection='stepwise',
#select="SL",Choose=NULL,SL=c(0.05,0.15,0.15),correct = "Bonferroni")

#for "B"
#Association(Tdata,Bdata,B,Elevels=NULL,selection='stepwise',
#select="SL",Choose=NULL,SL=c(0.05,0.15,0.15),correct = "Bonferroni")
}
\keyword{Association}

