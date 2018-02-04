predict.itrfit=function(obj, newdata=NULL, option='refine', delta=0.01) {
 inner=getInner(obj, newdata)
 if (option=='refine') {
  rule=pred_refine(inner, delta)
 } else {
  rule=pred(inner)
 }
 return(rule)
}