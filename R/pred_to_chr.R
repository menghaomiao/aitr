pred_to_chr=function(prediction) {
 n=nrow(prediction)
 rule=character(n)
 level=colnames(prediction)
 storage.mode(prediction)='logical'
 for (i in 1:n) rule[i]=paste(level[prediction[i, ]], collapse='_')
 return(rule)
}