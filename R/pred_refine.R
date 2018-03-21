pred_refine=function(inner, delta, fence) {
 rej=reject(inner, fence)
 inner=-inner/sapply(1:nrow(inner), function(i) min(inner[i, ]))
 res=refine(inner, delta)
 res[rej, ]=1
 colnames(res)=colnames(inner)
 return(res)
}