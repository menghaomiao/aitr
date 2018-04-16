pred=function(inner) {
 d=dim(inner)
 res=matrix(0, d[1], d[2])
 rule=sapply(1:d[1], function(i) which.max(inner[i, ]))
 res[cbind(1:d[1], rule)]=1
 colnames(res)=colnames(inner)
 rownames(res)=rownames(inner)
 class(res)=c('ITR', 'matrix')
 attr(res, 'outcome_ratio')=1
 return(res)
}