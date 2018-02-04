getInner.itrfit.svm=function(obj, newdata=NULL) {
 a=obj$a
 n=length(a)
 level=obj$level
 k=length(level)
 W_aW=matrix(-1/(k-1), n, k)
 W_aW[cbind(1:n, a)]=1
 theta_s_gamma=as.matrix(obj$theta_s_gamma)
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
 inner=get_inner(theta_s_gamma, K+1, W_aW)
 colnames(inner)=level
 return(inner)
}