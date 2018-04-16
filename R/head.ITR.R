head.ITR=function(rule) {
 s=attr(rule, 'outcome_ratio')
 cat(paste('Treatments with outcome ratio less than', s, 'are refined.\n'))
 head.matrix(rule)
}