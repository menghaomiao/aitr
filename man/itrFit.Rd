\name{itrFit}
\alias{itrFit}
\title{Individualized Treatment Rule (ITR) Fitting}
\description{
Fit the ITR with reject and refine options under bent svm or dwd loss.
}
\usage{
itrFit(x, a, y, p=1/k, method='svm', s=1.1, kernel='linear', epsilon=0.01, d=2, cv=FALSE, lambda=5^(-5:5))
}
\arguments{
\item{x}{A data frame or a matrix.}
\item{a}{A factor or a numeric vector indicating the treatment assigned to patients.}
\item{y}{A numeric vector indicating the outcome after taking a treatment.}
\item{p}{A numeric vector with each element between 0 to 1, indicating the probability of the treatment been assigned to the patient.}
\item{method}{The method to be used. Currently only \code{method='svm'} and \code{method='dwd'} are supported.}
\item{s}{A numeric value greater than 1, indicating the bent slope of the loss.}
\item{kernel}{Kernel type to be used. Currently only \code{linear}, \code{gaussian}, and \code{polynomial} are supported.}
\item{epsilon}{A positive numeric value used only when \code{kernel='gaussian'}.}
\item{d}{A positive integer used only when \code{kernel='polynomial'}.}
\item{cv}{Logical.If \code{TRUE} then a 5-fold cross validation is performed.}
\item{lambda}{A numeric vector indicating the penalty for tunning.}
}
\details{
This function could fit the Individualized Treatment Rule under bent svm or bent dwd loss. The bent slope is controled by parameter s. When \code{s=1}, it will degenerate to the usual ITR fit. The current version allows \code{s=1} for bent svm loss, but not for bent dwd loss. For the purpose of saving computation time, it suggested to set \code{cv=TRUE} for bent svm loss and \code{cv=FALSE} for bent dwd loss.
}
\value{
\code{itrFit} returns an object of class \code{c('itrfit.svm', 'itrfit')} or \code{c('itrfit.dwd', 'itrfit')}. An object of class \code{'itrfit'} is a list containing at least the following components:
\item{optlambda}{the optimal lambda.}
\item{s}{the bent slope.}
\item{level}{the treatment level.}
\item{x}{the data used for fitting.}
\item{a}{the treatment used for fitting.}
\item{y}{the outcome used for fitting.}
\item{kernel}{the kernel matrix computed.}

In addition, for \code{'itrfit.svm'} object it will return fitted value \code{theta_s_gamma}, and for \code{'itrfit.dwd'} it will return the fitted coefficient \code{coef}.
}
\author{
\packageAuthor{itrrnr}
}
\examples{
data(sim1)
train = sim1[1:1000, ]
test = sim1[3001:5000, ]
res = itrFit(train[, 1:2], train[, 3], train[, 4], cv = T) #using bent svm loss
res = itrFit(train[, 1:2], train[, 3], train[, 4], method = 'dwd', lambda = 10^(0:2)) #using bent dwd loss
}