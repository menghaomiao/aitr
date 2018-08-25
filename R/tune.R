tune=function(inner, a, y, p, c) {
 evalfit=function(rule) {
  rulesize=rowSums(rule)
  inrule=rowSums(rule*rulea)==1
  sum((y/p/(1+(rulesize-1)*c))[inrule])/sum(inrule/p/rulesize)
 }
 n=nrow(inner)
 k=ncol(inner)
 rulea=matrix(0, n, k)
 rulea[cbind(1:n, a)]=1
 stdinner=-inner/sapply(1:n, function(i) min(inner[i, ]))
 delta=seq(-1/2, 1/2, len=100)
 obj=numeric(100)
 for (i in 1:100) {
  rule=refine(stdinner, delta[i])
  obj[i]=evalfit(rule)
 }
 delta=delta[which.min(obj)]
 rule=refine(stdinner, delta)
 cutq=900:999/1000
 fence=-quantile(inner[inner<0], cutq)
 for (i in 1:100) {
  rej=reject(inner, fence[i])
  temprule=rule
  temprule[rej, ]=1
  obj[i]=evalfit(temprule)
 }
 fence=fence[which.min(obj)]
 rule[reject(inner, fence), ]=1
 list(refine_par=c(delta=delta, fence=fence), rule=rule, obj_value=min(obj))
}