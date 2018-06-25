ddwd=function(x) {
 ind=x>=-0.5
 ifelse(ind, 1, 0)+ifelse(ind, 0, 0.25/x^2)
}