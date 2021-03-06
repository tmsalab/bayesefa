% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{IFA_Mode_Jumper}
\alias{IFA_Mode_Jumper}
\title{Exploratory Item Factor Analysis for Ordinal Response Data}
\usage{
IFA_Mode_Jumper(Y, M, gamma, Ms, sdMH, burnin, chain_length = 10000L)
}
\arguments{
\item{Y}{A N by J matrix of item responses.}

\item{M}{The number of factors.}

\item{gamma}{Mode-jumping tuning parameter. Man & Culpepper used 0.5.}

\item{Ms}{A vector of the number of score categories. Note 2 = two parameter normal ogive, >2 ordinal normal ogive.}

\item{sdMH}{A vector of tuning parameters of length J for implementing the Cowles (1996) Metropolis-Hastings threshold sampler.}

\item{burnin}{Number of burn-in iterations to discard.}

\item{chain_length}{The total number of iterations (burn-in + post-burn-in).}
}
\value{
A list that contains nsamples = chain_length - burnin array draws from the posterior distribution:
\itemize{
\item \code{LAMBDA}: A J by M by nsamples array of sampled loading matrices on the standardized metric.
\item \code{PSIs}: A J by nsamples matrix of vector of variable uniquenesses on the standardized metric.
\item \code{ROW_OUT}: A matrix of sampled row indices of founding variables for mode-jumping algorithm.
\item \code{THRESHOLDS}: An array of sampled thresholds.
\item \code{INTERCEPTS}: Sampled variable thresholds on the standardized metric.
\item \code{ACCEPTED}: Acceptance rates for mode-jumping Metropolis-Hastings (MH) steps.
\item \code{MHACCEPT}: A J vector of acceptance rates for item threshold parameters. Note that binary items have an acceptance rate of zero, because MH steps are never performed.
}
}
\description{
Implement the Man and Culpepper (2020) mode-jumping algorithm to factor analyze ordinal response data. Missing values should be specified as a non-numeric value such as NA.
}
\examples{

#Load the psych package to apply the algorithm with the bfi dataset
library(psych)
data(bfi)

#subtract 1 from all items so scores start at 0
ydat<-as.matrix(bfi[,1:25])-1
J<-ncol(ydat)
N<-nrow(ydat)

#compute the number of item categories
Ms<-apply(ydat,2,max)+1

#specify burn-in and chain length
#in full application set these to 10000 and 20000
burnin = 10 
chain_length=20

out5<-IFA_Mode_Jumper(ydat,M=5,gamma=.5,Ms,sdMH=rep(.025,J),burnin,chain_length)

#check the acceptance rates for the threshold tuning parameter fall between .2 and .5
mh_acceptance_rates<-t(out5$MHACCEPT)

#compute mean thresholds
mthresholds<-apply(out5$THRESHOLDS,c(1,2),mean)

#compute mean intercepts
mintercepts<-apply(out5$INTERCEPTS,1,mean)

#rotate loadings to PLT
my_lambda_sample = out5$LAMBDA
my_M<-5
for (j in 1:dim(my_lambda_sample)[3]) {
  my_rotate = proposal2(1:my_M,my_lambda_sample[,,j],matrix(0,nrow=N,ncol = my_M))
  my_lambda_sample[,,j] = my_rotate$lambda
}

#compute mean of PLT loadings
mLambda<-apply(my_lambda_sample,c(1,2),mean)

#find promax rotation of posterior mean
rotatedLambda<-promax(mLambda)$loadings[1:J,1:my_M]

#save promax rotation matrix
promaxrotation<-promax(mLambda)$rotmat

#compute the factor correlation matrix
phi<-solve(t(promaxrotation)\%*\%promaxrotation)

}
\references{
Cowles, M. K. (1996), Accelerating Monte Carlo Markov chain convergence for cumulative link generalized linear models," Statistics and Computing, 6, 101-111.

Man, A. & Culpepper, S. A. (2020). A mode-jumping algorithm for Bayesian factor analysis. Journal of the American Statistical Association, doi:10.1080/01621459.2020.1773833.
}
\author{
Albert X Man, Steven Andrew Culpepper
}
