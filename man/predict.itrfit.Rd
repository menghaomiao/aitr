\name{predict.itrfit}
\alias{predict.itrfit}
\title{Predict Method for ITR Fit}
\description{
Predict the options based on \code{"itrfit"} object.
}
\usage{
## S3 method for class "itrfit"
predict(obj, newdata = NULL, option = "refine", delta = NULL, fence = NULL, c = obj$c)
}
\arguments{
\item{obj}{object of class inheriting from \code{"itrfit"}.}
\item{newdata}{an optional data frame or a data matrix in which to look for variables with which to predict. If omitted, the fitted options are used.}
\item{option}{the option used for prediction. \code{option = "normal"} will produce ITR while \code{option = "refine"} will produce A-ITR.}
\item{delta}{a numeric value greater than -1 used for predicting A-ITR, can be omitted if \code{option = "normal"} or \code{tunning = TRUE} when \code{itrFit} is called.}
\item{fence}{a numeric value greater than 0 used for predicting A-ITR, can be omitted if \code{option = "normal"} or \code{tunning = TRUE} when \code{itrFit} is called.}
\item{c}{a numeric value no less than 1, ignored unless \code{obj} is of class \code{"itrfit.dif"} and \code{obj$c = 1}.}
}
\details{
There are two prediction methods for \code{"itrfit"} object when \code{option = "refine"}. One is based on bent loss, where a one-step estimation is applied. Another is based on differentiable loss, which is a two-step estimation procedure. We recommand the first method in practice.
}
\value{
\code{predict.itrfit} preduces a 0-1 matrix of class \code{"ITR"}, where the column names are the treatment levels and ith row indicates which treatments are suggested to ith patient. Specifically, \code{1} stands for "suggested" and \code{0} stands for "not suggested".
}
\author{
\packageAuthor{aitr}
}
\examples{
data(sim1)
train = sim1[1:1000, ]
test = sim1[3001:5000, ]

res = itrFit(train[, 1:2], train[, 3], train[, 4], kernel = "polynomial", d = 4, lambda = 5^(4:11)) #bent hinge loss fit
rule = predict(res, test[, 1:2], option = "normal") #prediction for ITR
plot(test[, 1], test[, 2], col = pred_to_num(rule) + 1)

rule = predict(res, test[, 1:2]) #prediction for A-ITR using the one-step estimation
plot(test[, 1], test[, 2], col = pred_to_num(rule) + 1)

res = itrFit(train[, 1:2], train[, 3], train[, 4], c = 1, method = "log", kernel = "polynomial", d = 2, lambda = 5^(0:8)) #common logistic loss (unbent) fit
rule = predict(res, test[, 1:2], c = 1.2) #prediction for A-ITR using the two-step estimation
plot(test[, 1], test[, 2], col = pred_to_num(rule) + 1)
}