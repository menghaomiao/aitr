\name{itrFit}
\alias{itrFit}
\title{Fitting Individualized Treatment Recommendation (ITR)}
\description{
Fit alternative individualized treatment recommendation (A-ITR) or ITR given a loss function and a kernel.
}
\usage{
itrFit(x, a, y, p = 1/k, c = 1.2, method = "svm", kernel = "linear", epsilon = 1/2/median(d)^2, d = 2, lambda = 5^(-5:5), cv = FALSE, tunning = TRUE)
}
\arguments{
\item{x}{a data frame or a matrix.}
\item{a}{a factor or an integer vector indicating the treatment assigned to patients.}
\item{y}{a numeric vector indicating the outcome after taking the treatment.}
\item{p}{a numeric vector with each element between 0 to 1, indicating the probability of the treatment been assigned to the patient.}
\item{c}{a numeric value no less than 1, indicating the bent slope of the loss.}
\item{method}{the method to be used. Must be one of \code{"svm"}, \code{"dwd"}, \code{"exp"}, \code{"log"} and \code{"sq"}.}
\item{kernel}{kernel type to be used. Currently only \code{linear}, \code{gaussian}, and \code{polynomial} are supported.}
\item{epsilon}{a positive numeric value used only when \code{kernel = "gaussian"}.}
\item{d}{a positive integer used only when \code{kernel = "polynomial"}.}
\item{lambda}{a numeric vector indicating the penalty for tunning.}
\item{cv}{logical. If \code{TRUE} then a 5-fold cross validation is performed.}
\item{tunning}{logical. If \code{TRUE} then a tunning procedure for \code{delta} and \code{fence} (used for predicting A-ITR) is performed.}
}
\details{
This function fits an individualized treatment recommendation under a bent loss. The bent slope is controled by the parameter c. When \code{c = 1}, it degenerates to the traditional ITR fit. For the purpose of saving computation time, it is suggested to set \code{cv = TRUE} for bent svm loss and \code{cv = FALSE} for other bent losses. It is also suggested to set \code{tunning = TRUE} for fitting A-ITR.
}
\value{
\code{itrFit} returns an object of class \code{c("itrfit.svm", "itrfit")} or \code{c("itrfit.dif", "itrfit")}. An object of class \code{"itrfit"} is a list containing at least the following components:

\item{obj_value}{the expected outcome under the traditional ITR. If \code{tunning = TRUE}, the weighted expected outcome (for evaluating an A-ITR) is also returned.}
\item{optlambda}{the optimal lambda.}
\item{c}{the bent slope.}
\item{level}{the treatment levels.}
\item{x}{the input covariate \code{x} used for fitting.}
\item{a}{the input treatment \code{a} used for fitting.}
\item{y}{the input outcome \code{y} used for fitting.}
\item{method}{the method used.}
\item{kernel}{the kernel matrix computed.}
\item{predict}{an object of class \code{"ITR"}. If \code{tunning = TRUE}, return A-ITR, otherwise return ITR.}
\item{refine_par}{if \code{tunning = TRUE}, the tunning parameter used for predicting the A-ITR.}

In addition, the fitted coefficients are returned depending on the object class. Specifically, \code{theta_s_gamma} for \code{"itrfit.svm"} object, and \code{coef} for \code{"itrfit.dif"} object.
}
\author{
\packageAuthor{aitr}
}
\examples{
data(sim1)
data = sim1[1:1000, ]

res = itrFit(data[, 1:2], data[, 3], data[, 4], cv = TRUE) #using bent svm loss
res = itrFit(data[, 1:2], data[, 3], data[, 4], kernel = 'polynomial', d = 4, lambda = 5^(4:11)) #using bent svm loss with polynomial kernel
res = itrFit(data[, 1:2], data[, 3], data[, 4], method = 'square', kernel = 'polynomial', d = 2, lambda = 5^(0:6)) #using bent square loss
}