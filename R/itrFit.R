itrFit=function(x, a, y, p=1/k, method='svm', s=1.1, kernel='linear', epsilon=1/2/median(d)^2, d=2, cv=F, lambda=5^(-5:5)) {
 n=length(a)
 a=as.factor(a)
 level.name=levels(a)
 a=as.numeric(a)
 k=max(a)
 lambda=sort(lambda, decreasing=T)
 m=length(lambda)
 ind=as.matrix(dist(a, 'manhattan')) #ind=ind_{ai=al}-1/(k-1)*ind_{ai\neq al}
 ind[which(ind!=0)]=-1/(k-1)
 ind[which(ind==0)]=1
 w=y/p
 sminus=s-1
 x=as.matrix(x)
 if (kernel=='linear') {
  K=x%*%t(x)
 } else if (kernel=='gaussian') {
  d=as.matrix(dist(x))
  K=exp(-epsilon*d^2)
 } else if (kernel=='polynomial') {
  K=(1+x%*%t(x))^d
 } else stop('kernel must be one of the followings: linear, gaussian, polynomial.')
 if (method=='svm') {
  K=K+1
  WWK=ind*K
  diagK=diag(K)
  W_aW=matrix(-1/(k-1), n, k) #W_aW=<W_ai, W_j>
  W_aW[cbind(1:n, a)]=1
  if (cv) {
   inner=matrix(0, n, k*m)
   folds=sample.int(n)%%5+1
   for (i in 1:5) {
    id=folds==i
    theta_s_gamma=svmfit_C(WWK[!id, !id], diagK[!id], w[!id], sminus, lambda)
    inner[id, ]=get_inner(theta_s_gamma, K[id, !id], W_aW[!id, ])
   }
  } else {
   theta_s_gamma=svmfit_C(WWK, diagK, w, sminus, lambda)
   inner=get_inner(theta_s_gamma, K, W_aW)
  }
  dim(inner)=c(n, k, m)
 } else if (method=='dwd') {
  WWK=ind*(K+1)
  W=t(W.gen(k)[, a])
  inner=array(0, c(n, k, m))
  if (cv) {
   folds=sample.int(n)%%5+1
   for (i in 1:5) {
    id=folds==i
    A=dwdfit_C(WWK[!id, !id], K[!id, !id], W[!id, ], w[!id], sminus, lambda)
    for (j in 1:m) {
     inner[id, , j]=cbind(1, K[id, !id])%*%A[, , j]%*%W.gen(k)
    }
   }
  } else {
   A=dwdfit_C(WWK, K, W, w, sminus, lambda)
   for (j in 1:m) {
    inner[, , j]=cbind(1, K)%*%A[,  ,j]%*%W.gen(k)
   }
  }
 } else stop('current version only supports svm and dwd loss.')
 error=numeric(m)
 for (i in 1:m) {
  rule=sapply(1:n, function(j) which.max(inner[j, , i]))
  error[i]=sum(w[a==rule])/sum((a==rule)/p)
 }
 opt=which.min(error)
 if (opt==1 | opt==m) warning('The lambda may not be optimal. Another range may be considered.')
 optlambda=lambda[opt]
 res=list(optlambda=optlambda, s=s, level=level.name, x=x, a=a, y=y)
 if (method=='svm') {
  if (cv) theta_s_gamma=as.numeric(svmfit_C(WWK, diagK, w, sminus, optlambda))
  else theta_s_gamma=theta_s_gamma[, opt]
  res$theta_s_gamma=theta_s_gamma
  res$kernel=K-1
  class(res)=c('itrfit.svm', 'itrfit')
 } else {
  if (cv) A=dwdfit_C(WWK, K, W, w, sminus, lambda[1:opt])
  else A=A[, , opt]
  res$coef=A
  res$kernel=K
  class(res)=c('itrfit.dwd', 'itrfit')
 }
 attr(res$kernel, 'type')=kernel
 if (kernel=='gaussian') attr(res$kernel, 'epsilon')=epsilon
 if (kernel=='polynomial') attr(res$kernel, 'degree')=d
 return(res)
}