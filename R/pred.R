pred=function(inner) {
 n=nrow(inner)
 res=matrix(0, n, ncol(inner))
 rule=sapply(1:n, function(i) which.max(inner[i, ]))
 res[cbind(1:n, rule)]=1
 colnames(res)=colnames(inner)
 return(res)
}