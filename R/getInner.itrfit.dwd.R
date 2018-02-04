getInner.itrfit.dwd=function(obj, newdata=NULL) {
 A=obj$coef
 k=length(obj$level)
 K=obj$kernel
 if (!is.null(newdata)) {
  x=obj$x
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
 inner=cbind(1, K)%*%A%*%W.gen(k)
 return(inner)
}