// [[Rcpp::depends(RcppArmadillo)]]

#include<RcppArmadillo.h>
using namespace arma;

vec firdwd(vec inner, int n) {
 vec d=ones(n);
 double u;
 for (int i=0; i<n; i++) {
  u=inner[i];
  if (u<-0.5) d[i]=1/(4*u*u);
 }
 return d;
}

vec secdwd(vec inner, int n) {
 vec d=zeros(n);
 double u;
 for (int i=0; i<n; i++) {
  u=inner[i];
  if (u<-0.5) d[i]=-1/(2*u*u*u);
 }
 return d;
}

mat update_alpha(mat A, mat K, mat P, mat W, mat ZB, mat wW, mat wWW, double nlrho) {

 int n=K.n_rows, p=K.n_cols, kminus=W.n_cols;
 double partial, secpartial, temp, diff, epsilon=0.0000001*p*kminus;
 mat oldA(p, kminus);
 vec u=sum(W%(K*A), 1);

 for (int iter=0; iter<100; iter++) {

  oldA=A;

  /* update alpha0 */
  for (int j=0; j<kminus; j++) {

   /* update alpha_0j */
   for (int i=0; i<10; i++) {
    partial=sum(firdwd(u, n)%wW.col(j))+nlrho*A(0, j)+ZB(0, j);
    if (fabs(partial)<0.0000001) break;
    secpartial=sum(secdwd(u, n)%wWW.col(j))+nlrho;
    if (fabs(secpartial)<0.001) secpartial=0.001;
    temp=partial/secpartial;
    A(0, j)-=temp;
    u=u-temp*W.col(j);
   } /* newton raphson for alpha_0j */

  }

  /* update alpha */
  for (int j=0; j<kminus; j++) {
   for (int q=1; q<p; q++) {

    /* update alpha_qj */
    for (int i=0; i<10; i++) {
     partial=sum(firdwd(u, n)%wW.col(j)%K.col(q))+sum(P.col(q)%(nlrho*A.col(j)+ZB.col(j)));
     if (fabs(partial)<0.0000001) break;
     secpartial=sum(secdwd(u, n)%wWW.col(j)%K.col(q)%K.col(q))+nlrho*P(q, q);
     if (fabs(secpartial)<0.001) secpartial=0.001;
     temp=partial/secpartial;
     A(q, j)-=temp;
     u=u-temp*K.col(q)%W.col(j);
    } /* newton raphson for alpha_qj */

   }
  }

 diff=sum(sum(pow(A-oldA, 2)));
 if (diff<epsilon) break;
 }

 return A;
}

vec update_gamma(vec gamma, mat K, mat WWK, vec ub, vec w) {

 int n=K.n_rows;
 vec oldgamma(n);
 double gamma_new, diff, epsilon=0.0000001*n;;

 for (int iter=0; iter<100; iter++) {

  oldgamma=gamma;

  /* update gamma */
  for (int l=0; l<n; l++) {
   gamma_new=gamma[l]-(sum(gamma%WWK.col(l))-ub[l])/(K(l, l+1)+1);
   if (gamma_new<0) {
    gamma[l]=0;
   } else if (gamma_new>w[l]) {
    gamma[l]=w[l];
   } else {
    gamma[l]=gamma_new;
   }
  }

 diff=sum(pow(gamma-oldgamma, 2));
 if (diff<epsilon) break;
 }

 return gamma;
}

double update_rho(double rho, double r1, double r2) {
 double r=r1/r2;
 if (r>100) {
  rho=2*rho;
 } else if (r<0.01) {
  rho=rho/2;
 }
 return rho;
}

// [[Rcpp::export]]
cube dwdfit_C(mat WWK, mat K, mat W, vec w, double sminus, vec lambda, double maxiter=100) {

 int n=K.n_rows, p=n+1, kminus=W.n_cols, m=lambda.size();
 mat A(p, kminus), ZB(p, kminus), wW(n, kminus), wWW(n, kminus), P(p, p);
 vec u(n);
 double nlambda;
 cube coef(p, kminus, m);

 for (int j=0; j<kminus; j++) {
  wW.col(j)=w%W.col(j);
  wWW.col(j)=wW.col(j)%W.col(j);
 }
 P.fill(0), A.fill(0), ZB.fill(0), u.fill(1);
 P(0, 0)=1;
 P.submat(1, 1, p-1, p-1)=K;
 K=join_rows(u, K);

 if (sminus>0) { /* dwd fit for bent loss */

  mat B(p, kminus), Z(p, kminus), oldB(p, kminus), temp(p, kminus);
  vec b(n), ub(n), gamma(n);
  double nlrho, r1, r2, rho=1, epsilon=p*kminus*0.0000001;
  Z.fill(0), oldB.fill(0);

  /* iterate for lambda */
  for (int i=0; i<m; i++) {
 
   nlambda=n*lambda[i];

   /* ADMM algorithm start */
   for (int iter=0; iter<maxiter; iter++) {

    nlrho=nlambda+rho;

    /* update alpha */
    A=update_alpha(A, K, P, W, ZB, wW, wWW, nlrho);
    u=sum(W%(K*A), 1);
    b=sum(W%(K*Z), 1);
    ub=(rho*u+b)/sminus;

    /* update beta */
    gamma=update_gamma(gamma, K, WWK, ub, w);
    temp.fill(0);
    for (int j=0; j<kminus; j++) {
     temp.submat(1, j, n, j)=gamma%W.col(j);
    }
    temp.row(0)=sum(temp);
    B=A-(sminus*temp-Z)/rho;

    /* check convergence */
    r1=0, r2=0;
    temp=A-B;
    for (int j=0; j<kminus; j++) {
     r1+=as_scalar(temp.col(j).t()*P*temp.col(j));
    }
    temp=B-oldB;
    for (int j=0; j<kminus; j++) {
     r2+=as_scalar(temp.col(j).t()*P*temp.col(j));
    }
    r2=rho*rho*r2;
    if ((r1<epsilon) & (r2<epsilon)) break;
    oldB=B;

    /* update z */
    ZB=Z-rho*B;
    Z=ZB+rho*A;

    /* update rho */
    rho=update_rho(rho, r1, r2);
   }

   coef.slice(i)=A;
  }

 } else { /* dwd fit for common loss */

  /* iterate for lambda */
  for (int i=0; i<m; i++) {

   nlambda=n*lambda[i];
   A=update_alpha(A, K, P, W, ZB, wW, wWW, nlambda);
   coef.slice(i)=A;

  }
 }

 return coef;
}