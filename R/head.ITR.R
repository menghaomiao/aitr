head.ITR=function(rule) {
 c=attr(rule, 'outcome.ratio')
 cat(paste('Treatments with outcome ratio less than', c, 'are refined.\n'))
 head.matrix(rule)
}