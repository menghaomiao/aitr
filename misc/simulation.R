#######################################################################
##Simulation Setting: (Y is outcome, X is 2-dim predictor, A is action)
##X1, X2, X3, X4 from Unif(-1, 1) and error from N(0, 1/2)
##Case I:
##E(Y|X, A=1)=1+3*X1^2+3*X2^2
##E(Y|X, A=2)=3-X1^2/2-X2^2/2
##E(Y|X, A=3)=3+X1+X2
##Case II:
##E(Y|X, A=1)=2-cos(pi*(X1-X2)/2)
##E(Y|X, A=2)=2-cos(pi*(X1+X2)/2)
##E(Y|X, A=3)=2+cos(pi*(X1-X2)/2)
##E(Y|X, A=4)=2+cos(pi*(X1+X2)/2)
##Case III:
##E(Y|X, A=1)=max(2.5, 2.3+X1^2+X2^2)
##E(Y|X, A=2)=2.7-2*X1+exp(X3^2)-X4^3
##E(Y|X, A=3)=min(3, 3.2-X1^2-X2^2)
#######################################################################
library(aitr)
sim=function(n, p=5, setting=1) {
 x=matrix(2*runif(n*p)-1, n)
 y=sim.y(x, setting)
 k=ncol(y)
 y=y+matrix(rnorm(n*k)/4, n)
 a=sample(1:k, n, T)
 data=cbind(x, a, y)
 colnames(data)=c(paste0('x', 1:p), 'a', paste0('y', 1:k))
 return(data)
}
sim.y=function(x, setting=1) {
 x1=x[, 1];x2=x[, 2]
 if (setting==1) {
  y1=1+3*x1^2+3*x2^2
  y2=3-x1^2/2-x2^2/2
  y3=3+x1+x2
  truey=cbind(y1, y2, y3)
 } else if (setting==2) {
  y1=2-cos(pi*(x1-x2)/2)
  y2=2-cos(pi*(x1+x2)/2)
  y3=2+cos(pi*(x1-x2)/2)
  y4=2+cos(pi*(x1+x2)/2)
  truey=cbind(y1, y2, y3, y4)
 } else if (setting==3) {
  x3=x[, 3];x4=x[, 4]
  y1=pmax(2.5, 2.3+x1^2+x2^2)
  y2=2.7-2*x1+exp(x3^2)-x4^3
  y3=pmin(3, 3.2-x1^2-x2^2)
  truey=cbind(y1, y2, y3)
 } else stop('setting must be 1, 2, or 3!')
 return(truey)
}
getTR=function(truey, setting, c=1.2) {
 if (!missing(setting)) truey=sim.y(truey, setting)
 colnames(truey)=1:ncol(truey)
 truerule=truey/sapply(1:nrow(truey), function(i) min(truey[i, ]))
 temp=which(t(truerule)==1, arr.ind=T)
 temp=temp[, 2:1]
 truerule[truerule<=c]=1
 truerule[truerule>c]=0
 class(truerule)=c('ITR', 'matrix')
 attr(truerule, 'outcome.ratio')=c
 return(truerule)
}
getR=function(truerule) {
 region=rowSums(truerule)
 region[region==ncol(truerule)]=0
 region[region>1]=2
 return(region)
}
evalR=function(rule, truerule, testy, FUN=min) {
 f=match.fun(FUN)
 storage.mode(rule)='logical'
 y=sapply(1:nrow(testy), function(i) f(testy[i, ][rule[i, ]]))
 region=getR(truerule)
 tapply(y, region, mean)
}
regfit=function(x, a, y, testx) {
 k=max(a)
 sapply(1:k, function(i) randomForest::randomForest(x, y, testx, subset=a==i)$test$predicted)
}
linfit=function(x, a, y, testx) {
 k=max(a)
 predy=lapply(1:k, function(i) predict(lm(y~., data.frame(x), a==i), data.frame(testx)))
 do.call('cbind', predy)
}
comp=function(setting=1, p=5, n=1000) {
 compR=function(rule1, rule2, truerule, testy, testa) {
  r1min=evalR(rule1, truerule, testy)
  r2min=evalR(rule2, truerule, testy)
  r2max=evalR(rule2, truerule, testy, max)
  out=rbind(r1min, r2min, r2max)
  testy=testy[cbind(1:length(testa), testa)]
  loss1=evalITR(rule1, testa, testy)
  loss2=evalITR(rule2, testa, testy)
  cbind(out, loss=c(loss1, loss2, 0))
 }
 p=max(p, 2*(setting>2)+2)
 train=sim(n, p, setting)
 x=train[, 1:p];a=train[, p+1];y=train[, -(1:(p+1))][cbind(1:n, a)]
 test=sim(1000, p, setting)
 testx=test[, 1:p];testa=test[, p+1];testy=test[, -(1:(p+1))]
 bestrule=getTR(testx, setting, 1)
 truerule=getTR(testx, setting)
 out0=compR(bestrule, truerule, truerule, testy, testa)
 res=itrFit(x, a, y, kernel='polynomial', lambda=2^(-(2:9)), cv=T)
 rule1=predict(res, testx, 'normal')
 rule2=predict(res, testx)
 out1=compR(rule1, rule2, truerule, testy, testa)
#res=itrFit(x, a, y, c=1, method='sq', kernel='gaussian', lambda=4^(-(3:10)))
 if (setting==2) {
  res=itrFit(x, a, y, c=1, method='sq', kernel='polynomial', lambda=4^(-7:0))
 } else {
  res=itrFit(x, a, y, c=1, method='sq', kernel='polynomial', d=4, lambda=4^(-5:2))
 }
 rule1=predict(res, testx, 'normal')
 rule2=predict(res, testx, c=1.2)
 out2=compR(rule1, rule2, truerule, testy, testa)
#predy=linfit(x, a, y, testx)
 predy=regfit(x, a, y, testx)
 rule1=getTR(predy, c=1)
 rule2=getTR(predy)
 out3=compR(rule1, rule2, truerule, testy, testa)
 rbind(out0, out1, out2, out3)
}
comp()
comp(2)
comp(3)