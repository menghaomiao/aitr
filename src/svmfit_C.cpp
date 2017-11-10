#include<Rcpp.h>
using namespace Rcpp;

double box(double v, double b) {
 if (v<0) {
  v=0;
 } else if (v>b) {
  v=b;
 }
 return v;
}

// [[Rcpp::export]]
NumericMatrix svmfit_C(NumericMatrix WWK, NumericVector diagK, NumericVector w, double sminus, NumericVector lambda) {

/**************************************************************************
 This function computes (theta+gamma*(s-1))/n/lambda for different lambdas
 It return a matrix with n rows and m columns
 lambda: needs to be sorted from the largest to the smallest
 n: number of observations
 m: length of lambda
**************************************************************************/

 int n=WWK.ncol(), m=lambda.size();
 NumericMatrix theta_s_gamma(n, m);
 NumericVector theta(n), diff(n);
 double theta_new, nlambda, temp, epsilon=0.0000001*n;

 if (sminus>0) { /* svm fit wit bent hinge loss */

  NumericVector s_gamma(n), ws=w*sminus;
  double s_gamma_new;

  /* iterate for lambda */
  for (int i=0; i<m; i++) {

   nlambda=n*lambda[i];

   /* coordinate decent algorithm */
   for (int iter=0; iter<100; iter++) {

    /* update theta and gamma*/
    for (int l=0; l<n; l++) {

     temp=sum((theta+s_gamma)*WWK(_, l));
     s_gamma_new=box(s_gamma[l]-temp/diagK[l], ws[l]);
     diff[l]=s_gamma_new-s_gamma[l];
     s_gamma[l]=s_gamma_new;

     temp+=diff[l]*WWK(l, l)-nlambda;
     theta_new=box(theta[l]-temp/diagK[l], w[l]);
     diff[l]+=theta_new-theta[l];
     theta[l]=theta_new;

    } /* coordinate decent alorighm for one iteration ends */

    /* check convergence */
    if (sum(pow(diff, 2))<epsilon) break;

   } /* computation for one lambda ends */

   theta_s_gamma(_, i)=(theta+s_gamma)/nlambda;
  } /* iteration for all lambda ends */

 } else { /* svm fit with common hinge loss */

  /* iterate for lambda */
  for (int i=0; i<m; i++) {

   nlambda=n*lambda[i];

   for (int iter=0; iter<100; iter++) {

    /* update theta */
    for (int l=0; l<n; l++) {

     temp=sum(theta*WWK(_, l))-nlambda;
     theta_new=box(theta[l]-temp/diagK[l], w[l]);
     diff[l]=theta_new-theta[l];
     theta[l]=theta_new;

    } /* coordinate decent alorighm for one iteration ends */

    /* check convergence */
    if (sum(pow(diff, 2))<epsilon) break;

   } /* computation for one lambda ends */

   theta_s_gamma(_, i)=theta/nlambda;
  } /* iteration for all lambda ends */
 }

 return theta_s_gamma;
}