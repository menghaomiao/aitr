\name{itrFit}
\alias{itrFit}
\title{Individualized Treatment Rule (ITR) Fitting}
\description{
Fit the ITR with Reject and Refine Options (ITRRNR) under bent svm or dwd loss.
}
\usage{
itrFit(x, a, y, p = 1/k, s = 1.2, method = "svm", kernel = "linear", epsilon = 1/2/median(d)^2, d = 2, lambda = 5^(-5:5), cv = FALSE, tunning = TRUE)
}
\arguments{
\item{x}{a data frame or a matrix.}
\item{a}{a factor or an integer vector indicating the treatment assigned to patients.}
\item{y}{a numeric vector indicating the outcome after taking the treatment.}
\item{p}{a numeric vector with each element between 0 to 1, indicating the probability of the treatment been assigned to the patient.}
\item{s}{a numeric value greater than 1, indicating the bent slope of the loss.}
\item{method}{the method to be used. Currently only \code{method = "svm"} and \code{method = "dwd"} are supported.}
\item{kernel}{kernel type to be used. Currently only \code{linear}, \code{gaussian}, and \code{polynomial} are supported.}
\item{epsilon}{a positive numeric value used only when \code{kernel = "gaussian"}.}
\item{d}{a positive integer used only when \code{kernel = "polynomial"}.}
\item{lambda}{a numeric vector indicating the penalty for tunning.}
\item{cv}{logical.If \code{TRUE} then a 5-fold cross validation is performed.}
\item{tunning}{logical.If \code{TRUE} then a tunning procedure for \code{delta} and \code{fence} (used for predicting ITRRNR) is performed.}
}
\details{
This function could fit the Individualized Treatment Rule under bent svm or bent dwd loss. The bent slope is controled by parameter s. When \code{s = 1}, it will degenerate to the usual ITR fit. The current version allows \code{s = 1} for bent svm loss, but not for bent dwd loss. For the purpose of saving computation time, it suggested to set \code{cv = TRUE} for bent svm loss and \code{cv = FALSE} for bent dwd loss. It also suggested to set \code{tunning = TRUE} unless the reject and refine options are not needed.
}
\value{
\code{itrFit} returns an object of class \code{c("itrfit.svm", "itrfit")} or \code{c("itrfit.dwd", "itrfit")}. An object of class \code{"itrfit"} is a list containing at least the following components:

\item{obj_value}{the expected outcome under the original ITR. if \code{tunning=TRUE}, the weighted expected outcome under ITRRNR objective function is also returned.}
\item{optlambda}{the optimal lambda.}
\item{s}{the bent slope.}
\item{level}{the treatment levels.}
\item{x}{the input covariate \code{x} used for fitting.}
\item{a}{the input treatment \code{a} used for fitting.}
\item{y}{the input outcome \code{y} used for fitting.}
\item{kernel}{the kernel matrix computed.}
\item{predict}{an object of class \code{"ITR"}. If \code{tunning=TRUE}, then the fitted rule with reject and refine options is returned, otherwise the original ITR will be returned.}
\item{refine_par}{if \code{tunning=TRUE}, a tunning parameter used in predicting ITRRNR.}

In addition, for \code{"itrfit.svm"} object it will return fitted value \code{theta_s_gamma}, and for \code{"itrfit.dwd"} it will return the fitted coefficient \code{coef}.
}
\author{
\packageAuthor{itrrnr}
}
\examples{
data(sim1)
data = sim1[1:1000, ]

res = itrFit(data[, 1:2], data[, 3], data[, 4], cv = TRUE) #using bent svm loss
res = itrFit(data[, 1:2], data[, 3], data[, 4], method = "dwd", lambda = 10^(-1:2)) #using bent dwd loss
}