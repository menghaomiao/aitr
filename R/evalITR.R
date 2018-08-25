evalITR=function(rule, a, y, p=1/d[2]) {
 if (!is(rule, 'ITR')) stop('rule must be of class "ITR"!')
 d=dim(rule)
 c=attr(rule, 'outcome.ratio')
 if (!is.matrix(a)) {
  a=as.character(a)
  level=colnames(rule)
  if (all(unique(a)%in%level)) {
   a=match(a, level)
  } else {
   a=as.integer(a)
   maxa=max(a)
   if (is.na(maxa) | maxa>d[2]) stop('invalid second argument!')
  }
  if (length(a)!=d[1] | length(y)!=d[1]) stop('dimensions mismatch. Be sure to use the correct "a" and "y".')
  rulea=matrix(0, d[1], d[2])
  rulea[cbind(1:d[1], a)]=1
 } else rulea=a
 rulesize=rowSums(rule)
 inrule=rowSums(rule*rulea)==1
 sum((y/p/(1+(rulesize-1)*c))[inrule])/sum(inrule/p/rulesize)
}