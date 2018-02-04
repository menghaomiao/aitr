pred_refine=function(inner, delta=0.01) {
 res=matrix(0, nrow(inner), ncol(inner))
 inner[inner>=-delta & inner<=delta]=0
 res[inner>=0]=1
 colnames(res)=colnames(inner)
 return(res)
}