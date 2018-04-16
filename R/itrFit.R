itrFit=function(x, a, y, p=1/k, s=1.2, method='svm', kernel='linear', epsilon=1/2/median(d)^2, d=2, lambda=5^(-5:5), cv=F, tunning=T) {
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
 if (any(p<0 | p>1)) stop('p must be between 0 to 1!')
 if (length(s)>1 | s<1) stop('s must be a number greater than 1!')
 w=y/p
 sminus=s-1
 x=as.matrix(x)
 if (kernel=='linear') {
  ker=x%*%t(x)
 } else if (kernel=='gaussian') {
  d=as.matrix(dist(x))
  ker=exp(-epsilon*d^2)
 } else if (kernel=='polynomial') {
  ker=(1+x%*%t(x))^d
 } else stop('kernel must be one of the followings: linear, gaussian, polynomial.')
 K=ker+diag(n)*0.01*mean(diag(ker))
 WWK=ind*(K+1)
 if (method=='svm') {
  diagK=diag(K)+1
  W_aW=matrix(-1/(k-1), n, k) #W_aW=<W_ai, W_j>
  W_aW[cbind(1:n, a)]=1
  if (cv) {
   inner=matrix(0, n, k*m)
   folds=sample.int(n)%%5+1
   for (i in 1:5) {
    id=folds==i
    theta_s_gamma=svmfit_C(WWK[!id, !id], diagK[!id], w[!id], sminus, lambda)
    inner[id, ]=get_inner(theta_s_gamma, ker[id, !id], W_aW[!id, ])
   }
  } else {
   theta_s_gamma=svmfit_C(WWK, diagK, w, sminus, lambda)
   inner=get_inner(theta_s_gamma, ker, W_aW)
  }
  dim(inner)=c(n, k, m)
 } else if (method=='dwd') {
  Wbasis=W.gen(k)
  W=t(Wbasis[, a])
  inner=array(0, c(n, k, m))
  if (cv) {
   folds=sample.int(n)%%5+1
   for (i in 1:5) {
    id=folds==i
    A=dwdfit_C(WWK[!id, !id], K[!id, !id], W[!id, ], w[!id], sminus, lambda)
    for (j in 1:m) {
     inner[id, , j]=cbind(1, ker[id, !id])%*%A[, , j]%*%Wbasis
    }
   }
  } else {
   A=dwdfit_C(WWK, K, W, w, sminus, lambda)
   for (j in 1:m) {
    inner[, , j]=cbind(1, ker)%*%A[,  ,j]%*%Wbasis
   }
  }
 } else stop('current version only supports svm and dwd loss.')
 obj=numeric(m)
 for (i in 1:m) {
  rule=sapply(1:n, function(j) which.max(inner[j, , i]))
  obj[i]=sum(w[a==rule])/sum((a==rule)/p)
 }
 opt=which.min(obj)
 if (opt==1 | opt==m) warning('The lambda may not be optimal. Another range may be considered.')
 optlambda=lambda[opt]
 optinner=inner[, , opt]
 res=list(obj_value=c(itr=obj[opt]), optlambda=optlambda, s=s, level=level.name, x=x, a=a, y=y, kernel=ker)
 if (method=='svm') {
  if (cv) {
   theta_s_gamma=svmfit_C(WWK, diagK, w, sminus, optlambda)
   inner=get_inner(theta_s_gamma, ker, W_aW)
  } else {
   theta_s_gamma=theta_s_gamma[, opt]
   inner=optinner
  }
  res$theta_s_gamma=as.numeric(theta_s_gamma)
  class(res)=c('itrfit.svm', 'itrfit')
 } else {
  if (cv) {
   A=dwdfit_C(WWK, K, W, w, sminus, lambda[1:opt])
   inner=cbind(1, ker)%*%A%*%Wbasis
  } else {
   A=A[, , opt]
   inner=optinner
  }
  res$coef=A
  class(res)=c('itrfit.dwd', 'itrfit')
 }
 attr(res$kernel, 'type')=kernel
 if (kernel=='gaussian') attr(res$kernel, 'epsilon')=epsilon
 if (kernel=='polynomial') attr(res$kernel, 'degree')=d
 if (s>1 & tunning) {
  tunning=tune(optinner, a, y, p, s)
  obj=tunning$obj_value
  res$obj_value=c(res$obj_value, itrrnr=obj)
  if (cv) tunning=tune(inner, a, y, p, s)
  res$predict=tunning$rule
  colnames(res$predict)=level.name
  rownames(res$predict)=rownames(x)
  attr(res$predict, 'outcome_ratio')=s
  class(res$predict)=c('ITR', 'matrix')
  res$refine_par=tunning$refine_par
 } else {
  colnames(inner)=level.name
  rownames(inner)=rownames(x)
  res$predict=pred(inner)
 }
 return(res)
}