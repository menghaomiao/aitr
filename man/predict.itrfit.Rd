\name{predict.itrfit}
\alias{predict.itrfit}
\title{Predict method for Individualized Treatment Rule (ITR)}
\description{
Predict the options based on itrfit object.
}
\usage{
## S3 method for class \code{"itrfit"}
predict(obj, newdata = NULL, option = "refine", delta = NULL, fence = NULL)
}
\arguments{
\item{obj}{object of class inheriting from \code{"itrrnr"}.}
\item{newdata}{an optional data frame or data matrix in which to look for variables with which to predict. If omitted, the fitted options are used.}
\item{option}{the option used for prediction. \code{option = "normal"} will produce ITR while \code{option = "refine"} will produce ITR with reject and refine option.}
\item{delta}{a numeric value greater than -1 used for predicting ITRRNR, can be omitted if \code{option = "normal"} or \code{tunning = TRUE} when \code{itrFit} is called.}
\item{fence}{a numeric value greater than 0 used for predicting ITRRNR, can be omitted if \code{option = "normal"} or \code{tunning = TRUE} when \code{itrFit} is called.}
}
\value{
\code{predict.itrfit} preduces a 0-1 matrix, where the column names are the treatment levels and ith row indicates which treatments are suggested to ith patient. Specifically, 1 stands for "suggested" and 0 stands for "not suggested".
}
\author{
\packageAuthor{itrrnr}
}
\examples{
data(sim1)
train = sim1[1:1000, ]
test = sim1[3001:5000, ]
res = itrFit(train[, 1:2], train[, 3], train[, 4], kernel = "polynomial", d = 4, lambda = 5^(4:11))

rule = predict(res, test[, 1:2], option = "normal") #prediction for ITR
plot(test[, 1], test[, 2], col = pred_to_num(rule) + 1)

rule = predict(res, test[, 1:2]) #prediction for ITR with Reject and Rifine options
plot(test[, 1], test[, 2], col = pred_to_num(rule) + 1)
}