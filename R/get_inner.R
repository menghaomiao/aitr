get_inner=function(alpha_s_gamma, XX, W_aW) {
## inner=<W_j, f(x_i)>, k*1 vector for each observation
## n*k matrix for each lambda, column bind for all lambdas
 alpha_s_gamma=as.matrix(alpha_s_gamma)
 m=ncol(alpha_s_gamma)
 inner=lapply(1:m, function(j) -XX%*%(alpha_s_gamma[, j]*W_aW))
 do.call('cbind', inner)
}