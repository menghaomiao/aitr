ddel=function(x, method) {
 switch(method, dwd=ddwd(x), exponential=exp(x), logistic=1-1/(1+exp(x)))
}

ddwd=function(x) {
 ind=x>=-0.5
 ifelse(ind, 1, 0)+ifelse(ind, 0, 0.25/x^2)
}