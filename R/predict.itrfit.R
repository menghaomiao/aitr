predict.itrfit=function(obj, newdata=NULL, option='refine', delta=0.01) {
 level=obj$level
 k=length(level)
 method=class(obj)[1]
 x=obj$x
 K=obj$kernel
 if (!is.null(newdata)) {
  newx=as.matrix(newdata)
  kernel=attr(K, 'type')
  if (kernel=='linear') {
   K=newx%*%t(x)
  } else if (kernel=='gaussian') {
   epsilon=attr(K, 'epsilon')
   K=exp(-epsilon*as.matrix(pdist(newx, x))^2)
  } else if (kernel=='polynomial') {
   d=attr(K, 'degree')
   K=(1+newx%*%t(x))^d
  }
 }
 if (method=='itrfit.svm') {
  a=obj$a
  n=length(a)
  W_aW=matrix(-1/(k-1), n, k)
  W_aW[cbind(1:n, a)]=1
  theta_s_gamma=obj$theta_s_gamma
  inner=get_inner(as.matrix(theta_s_gamma), K+1, W_aW)
 } else if (method=='itrfit.dwd') {
  A=obj$coef
  inner=cbind(1, K)%*%A%*%W.gen(k)
 }
 if (option=='refine') {
  rule=pred_refine(inner, delta)
 } else {
  rule=pred(inner)
 }
 colnames(rule)=level
 return(rule)
}