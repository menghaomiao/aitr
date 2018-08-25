predict.itrfit.del=function(obj, newdata=NULL, option='refine', c=obj$c, delta=NULL, fence=NULL) {
 if (obj$c==1 & option=='refine') {
  inner=getInner(obj, newdata)
  d=dim(inner)
  inner=getInner(obj, newdata)
  rule=matrix(0, d[1], d[2])
  deri=ddel(inner, obj$method)
  maxd=sapply(1:d[1], function(i) max(deri[i, ]))
  rule[maxd/deri<=c]=1
  colnames(rule)=colnames(inner)
  rownames(rule)=rownames(inner)
  class(rule)=c('ITR', 'matrix')
  attr(rule, 'outcome.ratio')=c
  return(rule)
 } else {
  if (c!=obj$c & option=='refine') {
   message(paste0('no refine method for c=', c, '. Use c=', obj$c, ' instead.'))
  }
  NextMethod()
 }
}