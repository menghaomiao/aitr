getInner.itrfit.dwd=function(obj, newdata=NULL) {
 level=obj$level
 k=length(level)
 A=obj$coef
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
 colnames(inner)=level
 return(inner)
}