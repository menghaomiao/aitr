pred_to_num=function(prediction) {
 d=dim(prediction)
 res=numeric(d[1])
 temp=rowSums(prediction)
 ind=which(temp==1)
 if (length(ind)>0) res[ind]=sapply(ind, function(i) which(prediction[i, ]==1))
 res[temp>=2]=d[2]+1
 res[temp==d[2]]=0
 return(res)
}