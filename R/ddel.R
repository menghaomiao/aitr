ddel=function(x, method) {
 switch(method, dwd=ddwd(x), exponential=exp(x/5)/5, logistic=exp(x/2)/(1+exp(x/2))/2)
}

ddwd=function(x) {
 ind=x>=-0.5
 ifelse(ind, 1, 0)+ifelse(ind, 0, 0.25/x^2)
}