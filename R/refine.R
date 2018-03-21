refine=function(inner, delta) {
 rule=matrix(0, nrow(inner), ncol(inner))
 rule[inner>=delta]=1
 ind=rowSums(rule)==0
 if (sum(ind)>0) rule[ind, ][inner[ind, ]>=-delta]=1
 return(rule)
}