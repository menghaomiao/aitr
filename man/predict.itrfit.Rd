\name{predict.itrfit}
\alias{predict.itrfit}
\title{Predict method for Individualized Treatment Rule (ITR)}
\description{
Predict the options based on itrfit object.
}
\usage{
## S3 method for class 'itrfit'
itrFit(obj, newdata = NULL, option = "refine", delta = 0.01)
}
\arguments{
\item{obj}
Object of class inheriting from \code{'itrrnr'}.
\item{newdata}
An optional data frame or data matrix in which to look for variables with which to predict. If omitted, the fitted options are used.
\item{option}
The option used for prediction. \code{option='normal'} will produce ITR while \code{option='refine'} will produce ITR with reject and refine option.
\item{delta}
A positive numeric value for tunning refine options. Used only when \code{option='refine'}.
}
\value{
\code{predict.itrfit} preduces a matrix, where the column namea are the treatment levels and ith row indicates which treatments the ith patient should take. Specifically, 1 stands for "belong to" and 0 stands for "not belong to".
}
\author{
\packageAuthor{itrrnr}
}
\examples{
data(sim1)
train = sim1[1:1000, ]
test = sim1[3001:5000, ]
res = itrFit(train[, 1:2], train[, 3], train[, 4], kernel = 'polynomial', d = 4, lambda = 5^(4:11))

rule = predict(res, test[, 1:2], option = 'normal') #prediction for ITR
plot(test[, 1], test[, 2], col = pred_to_num(rule) + 1)

rule = predict(res, test[, 1:2], delta = 0.01) #prediction for ITR with Reject and Rifine options
plot(test[, 1], test[, 2], col = pred_to_num(rule) + 1)
}