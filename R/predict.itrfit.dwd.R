predict.itrfit.dwd=function(obj, newdata=NULL, option='refine', s=obj$s, delta=NULL, fence=NULL) {
 if (obj$s==1 & option=='refine') {
  inner=getInner(obj, newdata)
  rule=matrix(0, nrow(inner), ncol(inner))
  d=ddwd(inner)
  rule[1/d<=s]=1
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