evalFit=function(rule, testrule, testw, s) {
 mean(testw*(rowSums(rule*testrule)==1)/(1+(rowSums(rule)-1)*s))
}