pred_refine=function(inner, delta=0.01) {
 rule=matrix(0, nrow(inner), ncol(inner))
 inner[inner>=-delta & inner<=delta]=0
 rule[inner>=0]=1
 return(rule)
}