reject=function(inner, fence) {
 rej=inner>=-fence & inner<=fence
 rowSums(rej)==ncol(inner)
}