predict.itrfit.del=function(obj, newdata=NULL, option='refine', s=obj$s, delta=NULL, fence=NULL) {
 if (obj$s==1 & option=='refine') {
  inner=getInner(obj, newdata)
  d=dim(inner)
  inner=getInner(obj, newdata)
  rule=matrix(0, d[1], d[2])
  deri=ddel(inner, obj$method)
  maxd=sapply(1:d[1], function(i) max(deri[i, ]))
  rule[maxd/deri<=s]=1
  colnames(rule)=colnames(inner)
  rownames(rule)=rownames(inner)
  class(rule)=c('ITR', 'matrix')
  attr(rule, 'outcome.ratio')=s
  return(rule)
 } else {
  if (s!=obj$s & option=='refine') {
   message(paste0('no refine method for s=', s, '. Use s=', obj$s, ' instead.'))
  }
  NextMethod()
 }
}