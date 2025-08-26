#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector remove_first(NumericVector x) {
  int n = x.length();
  IntegerVector idx = seq(1, n-1);
  return x[idx];
}

// [[Rcpp::export]]
NumericVector subset(NumericVector x, int min, int max) {
  IntegerVector idx = seq(min, max);
  return x[idx];
}



// [[Rcpp::export]]
NumericVector ll_Poisson_Gamma(NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas){
  
  int r = beta1.length() -1;
  int s = gamma1.length();
  int p = beta2.length() -1;
  int k = gamma2.length()-1;
  
  NumericVector u_aux(4);
  
  u_aux[0] =r;
  u_aux[1] =s;
  u_aux[2] =p;
  u_aux[3] = k;
  
  int u = max(u_aux);
  int N= Y.nrow();
  
  NumericVector Y1= Y(_,0);
  NumericVector Y2= Y(_,1);
  
  NumericVector ll1(N); NumericVector ll2(N);
  NumericVector mu1(N); NumericVector mu2(N);
  
  NumericVector beta1_aux = remove_first(beta1);
  NumericVector beta2_aux = remove_first(beta2);
  NumericVector delta1 = delta_gammas[0];
  NumericVector delta2 = delta_gammas[1];
  
 
 for(int t = u; t < N; t++) {
 
    
    NumericVector sub_y1 = log(rev(subset(Y1, t-r, t-1)));
    double m1 = sum(sub_y1*beta1_aux);
    NumericVector sub_y2 = rev(subset(Y2, t-s, t-1))+1;
    double m2 = sum(delta1*gamma1*log(sub_y2));
    
    mu1[t] = exp(beta1[0] + m1 +m2);
    ll1[t]= R::dgamma(Y1[t],1/phi,(mu1[t]*phi),true);
    
    NumericVector sub3 = rev(subset(Y2, t-p, t-1)) +1;
    NumericVector sub4 = rev(subset(Y1, t-k, t));
   double m3 =  sum(log(sub3)*beta2_aux);
   double m4 = sum(delta2*gamma2*log(sub4));
    
    mu2[t] = exp(beta2[0] + m3 +m4);
    
    ll2[t] = R::dpois(Y2[t], mu2[t], true);   
       
      
    
  }  
 
 NumericVector output(2);
 output[0] = sum(ll1);
 output[1] = sum(ll2);
 
  return output;
  
    
}


// [[Rcpp::export]]
NumericVector ll_Gamma_Poisson(NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas){
  
  int r = beta1.length() -1;
  int s = gamma1.length();
  int p = beta2.length() -1;
  int k = gamma2.length()-1;
  
  NumericVector u_aux(4);
  
  u_aux[0] =r;
  u_aux[1] =s;
  u_aux[2] =p;
  u_aux[3] = k;
  
  int u = max(u_aux);
  int N= Y.nrow();
  
  NumericVector Y1= Y(_,0);
  NumericVector Y2= Y(_,1);
  
  NumericVector ll1(N); NumericVector ll2(N);
  NumericVector mu1(N); NumericVector mu2(N);
  
  NumericVector beta1_aux = remove_first(beta1);
  NumericVector beta2_aux = remove_first(beta2);
  NumericVector delta1 = delta_gammas[0];
  NumericVector delta2 = delta_gammas[1];
  
  
  for(int t = u; t < N; t++) {
    
    
    NumericVector sub_y1 = log(rev(subset(Y1, t-r, t-1))+1);
    double m1 = sum(sub_y1*beta1_aux);
    NumericVector sub_y2 = rev(subset(Y2, t-s, t-1));
    double m2 = sum(delta1*gamma1*log(sub_y2));
    
    mu1[t] = exp(beta1[0] + m1 +m2);
    ll1[t]= R::dpois(Y1[t], mu1[t],true);
    
    NumericVector sub3 = rev(subset(Y2, t-p, t-1));
    NumericVector sub4 = rev(subset(Y1, t-k, t))+1;
    double m3 =  sum(log(sub3)*beta2_aux);
    double m4 = sum(delta2*gamma2*log(sub4));
    
    mu2[t] = exp(beta2[0] + m3 +m4);
    
    ll2[t] = R::dgamma(Y2[t], 1/phi,(mu2[t]*phi),true);
    
    
    
  }  
  
  NumericVector output(2);
  output[0] = sum(ll1);
  output[1] = sum(ll2);
  
  return output;
  
  
}


// [[Rcpp::export]]
NumericVector ll_Geo_Gamma(NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas){
  
  int r = beta1.length() -1;
  int s = gamma1.length();
  int p = beta2.length() -1;
  int k = gamma2.length()-1;
  
  NumericVector u_aux(4);
  
  u_aux[0] =r;
  u_aux[1] =s;
  u_aux[2] =p;
  u_aux[3] = k;
  
  int u = max(u_aux);
  int N= Y.nrow();
  
  NumericVector Y1= Y(_,0);
  NumericVector Y2= Y(_,1);
  
  NumericVector ll1(N); NumericVector ll2(N);
  NumericVector mu1(N); NumericVector mu2(N);
  
  NumericVector beta1_aux = remove_first(beta1);
  NumericVector beta2_aux = remove_first(beta2);
  NumericVector delta1 = delta_gammas[0];
  NumericVector delta2 = delta_gammas[1];
  
  
  for(int t = u; t < N; t++) {
    
    
    NumericVector sub_y1 = log(rev(subset(Y1, t-r, t-1)));
    double m1 = sum(sub_y1*beta1_aux);
    NumericVector sub_y2 = rev(subset(Y2, t-s, t-1))+1;
    double m2 = sum(delta1*gamma1*log(sub_y2));
    
    mu1[t] = exp(beta1[0] + m1 +m2);
    ll1[t]= R::dgamma(Y1[t],1/phi,(mu1[t]*phi),true);
    
    NumericVector sub3 = rev(subset(Y2, t-p, t-1)) +1;
    NumericVector sub4 = rev(subset(Y1, t-k, t));
    double m3 =  sum(log(sub3)*beta2_aux);
    double m4 = sum(delta2*gamma2*log(sub4));
    
    mu2[t] = beta2[0] + m3 +m4;
    
    ll2[t] =  R::dgeom(Y2[t], exp(mu2[t])/(1+exp(mu2[t])), true);
    
    
    
  }  
  
  NumericVector output(2);
  output[0] = sum(ll1);
  output[1] = sum(ll2);
  
  return output;
  
  
}


// [[Rcpp::export]]
NumericVector ll_Geo_LN(NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas){
  
  int r = beta1.length() -1;
  int s = gamma1.length();
  int p = beta2.length() -1;
  int k = gamma2.length()-1;
  
  NumericVector u_aux(4);
  
  u_aux[0] =r;
  u_aux[1] =s;
  u_aux[2] =p;
  u_aux[3] = k;
  
  int u = max(u_aux);
  int N= Y.nrow();
  
  NumericVector Y1= Y(_,0);
  NumericVector Y2= Y(_,1);
  
  NumericVector ll1(N); NumericVector ll2(N);
  NumericVector mu1(N); NumericVector mu2(N);
  
  NumericVector beta1_aux = remove_first(beta1);
  NumericVector beta2_aux = remove_first(beta2);
  NumericVector delta1 = delta_gammas[0];
  NumericVector delta2 = delta_gammas[1];
  
  
  for(int t = u; t < N; t++) {
    
    
    NumericVector sub_y1 = log(rev(subset(Y1, t-r, t-1)));
    double m1 = sum(sub_y1*beta1_aux);
    NumericVector sub_y2 = rev(subset(Y2, t-s, t-1))+1;
    double m2 = sum(delta1*gamma1*log(sub_y2));
    
    mu1[t] = exp(beta1[0] + m1 +m2);
    ll1[t]=  R::dlnorm(Y1[t], log(mu1[t]), phi, true);
      
    NumericVector sub3 = rev(subset(Y2, t-p, t-1)) +1;
    NumericVector sub4 = rev(subset(Y1, t-k, t));
    double m3 =  sum(log(sub3)*beta2_aux);
    double m4 = sum(delta2*gamma2*log(sub4));
    
    mu2[t] = beta2[0] + m3 +m4;
   
    ll2[t] =  R::dgeom(Y2[t], exp(mu2[t])/(1+exp(mu2[t])), true);
    
    
    
  }  
  
  NumericVector output(2);
  output[0] = sum(ll1);
  output[1] = sum(ll2);
  
  return output;
  
  
}


// [[Rcpp::export]]
NumericVector ll_Geo_Geo(NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas){
  
  int r = beta1.length() -1;
  int s = gamma1.length();
  int p = beta2.length() -1;
  int k = gamma2.length()-1;
  
  NumericVector u_aux(4);
  
  u_aux[0] =r;
  u_aux[1] =s;
  u_aux[2] =p;
  u_aux[3] = k;
  
  int u = max(u_aux);
  int N= Y.nrow();
  
  NumericVector Y1= Y(_,0);
  NumericVector Y2= Y(_,1);
  
  NumericVector ll1(N); NumericVector ll2(N);
  NumericVector mu1(N); NumericVector mu2(N);
  
  NumericVector beta1_aux = remove_first(beta1);
  NumericVector beta2_aux = remove_first(beta2);
  NumericVector delta1 = delta_gammas[0];
  NumericVector delta2 = delta_gammas[1];
  
  
  for(int t = u; t < N; t++) {
    
    NumericVector sub_y1 = log(rev(subset(Y1, t-r, t-1)) +1);
    double m1 = sum(sub_y1*beta1_aux);
    NumericVector sub_y2 = rev(subset(Y2, t-s, t-1))+1;
    double m2 = sum(delta1*gamma1*log(sub_y2));
    
    mu1[t] = beta1[0] + m1 +m2;
    ll1[t]=  R::dgeom(Y1[t], exp(mu1[t])/(1+exp(mu1[t])), true);
    
    NumericVector sub3 = rev(subset(Y2, t-p, t-1)) +1;
    NumericVector sub4 = rev(subset(Y1, t-k, t));
    double m3 =  sum(log(sub3)*beta2_aux);
    double m4 = sum(delta2*gamma2*log(sub4+1));
    
    mu2[t] = beta2[0] + m3 +m4;
    
    ll2[t] =  R::dgeom(Y2[t], exp(mu2[t])/(1+exp(mu2[t])), true);
    
  }  
  
  NumericVector output(2);
  output[0] = sum(ll1);
  output[1] = sum(ll2);
  
  return output;
  
  
}



// [[Rcpp::export]]
NumericVector ll_Pois_Pois(NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas){
  
  int r = beta1.length() -1;
  int s = gamma1.length();
  int p = beta2.length() -1;
  int k = gamma2.length()-1;
  
  NumericVector u_aux(4);
  
  u_aux[0] =r;
  u_aux[1] =s;
  u_aux[2] =p;
  u_aux[3] = k;
  
  int u = max(u_aux);
  int N= Y.nrow();
  
  NumericVector Y1= Y(_,0);
  NumericVector Y2= Y(_,1);
  
  NumericVector ll1(N); NumericVector ll2(N);
  NumericVector mu1(N); NumericVector mu2(N);
  
  NumericVector beta1_aux = remove_first(beta1);
  NumericVector beta2_aux = remove_first(beta2);
  NumericVector delta1 = delta_gammas[0];
  NumericVector delta2 = delta_gammas[1];
  
  
  for(int t = u; t < N; t++) {
    
    NumericVector sub_y1 = log(rev(subset(Y1, t-r, t-1)) +1);
    double m1 = sum(sub_y1*beta1_aux);
    NumericVector sub_y2 = rev(subset(Y2, t-s, t-1))+1;
    double m2 = sum(delta1*gamma1*log(sub_y2));
    
    mu1[t] = exp(beta1[0] + m1 +m2);
    ll1[t]=  R::dpois(Y1[t], mu1[t], true);
    
    NumericVector sub3 = rev(subset(Y2, t-p, t-1)) +1;
    NumericVector sub4 = rev(subset(Y1, t-k, t));
    double m3 =  sum(log(sub3)*beta2_aux);
    double m4 = sum(delta2*gamma2*log(sub4 + 1));
    
    mu2[t] = exp(beta2[0] + m3 +m4);
    
    ll2[t] =  R::dpois(Y2[t], mu2[t], true);
    
  }  
  
  NumericVector output(2);
  output[0] = sum(ll1);
  output[1] = sum(ll2);
  
  return output;
  
  
}


// [[Rcpp::export]]
List mh_gamma1(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_g1, NumericVector prior_sd_g1){

  
  NumericVector gamma1_prime = clone(gamma1);
  
  gamma1_prime[g_index] = R::rnorm(gamma1[g_index], sigma_g1[g_index]);
  
  double ll_diff = ll_Poisson_Gamma( beta1,  beta2,  gamma1_prime,  gamma2,  phi,  Y,  delta_gammas)[0]-
                   ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0];

  double prior_diff = R::dnorm(gamma1_prime[g_index],0, prior_sd_g1[g_index],true)-
                         R::dnorm(gamma1[g_index],0, prior_sd_g1[g_index],true);

  double alpha = ll_diff+prior_diff;

  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;

if( compare[0]==true ){
  output = gamma1_prime;
  accept=1;
} else {
  output = gamma1;
  accept = 0;
}
return List::create(_["param"] = output, 
                    _["accepted"] = accept);

}



// [[Rcpp::export]]
List mh_gamma1_ga_pois(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_g1, NumericVector prior_sd_g1){
  
  
  NumericVector gamma1_prime = clone(gamma1);
  
  gamma1_prime[g_index] = R::rnorm(gamma1[g_index], sigma_g1[g_index]);
  
  double ll_diff = ll_Gamma_Poisson( beta1,  beta2,  gamma1_prime,  gamma2,  phi,  Y,  delta_gammas)[0]-
    ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0];
  
  double prior_diff = R::dnorm(gamma1_prime[g_index],0, prior_sd_g1[g_index],true)-
    R::dnorm(gamma1[g_index],0, prior_sd_g1[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = gamma1_prime;
    accept=1;
  } else {
    output = gamma1;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}



// [[Rcpp::export]]
List mh_gamma1_geoln(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_g1, NumericVector prior_sd_g1){
  
  
  NumericVector gamma1_prime = clone(gamma1);
  
  gamma1_prime[g_index] = R::rnorm(gamma1[g_index], sigma_g1[g_index]);
  
  double ll_diff = ll_Geo_LN( beta1,  beta2,  gamma1_prime,  gamma2,  phi,  Y,  delta_gammas)[0]-
    ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0];
  
  double prior_diff = R::dnorm(gamma1_prime[g_index],0, prior_sd_g1[g_index],true)-
    R::dnorm(gamma1[g_index],0, prior_sd_g1[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = gamma1_prime;
    accept=1;
  } else {
    output = gamma1;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}


// [[Rcpp::export]]
List mh_gamma1_geogamma(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_g1, NumericVector prior_sd_g1){
  
  
  NumericVector gamma1_prime = clone(gamma1);
  
  gamma1_prime[g_index] = R::rnorm(gamma1[g_index], sigma_g1[g_index]);
  
  double ll_diff = ll_Geo_Gamma( beta1,  beta2,  gamma1_prime,  gamma2,  phi,  Y,  delta_gammas)[0]-
    ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0];
  
  double prior_diff = R::dnorm(gamma1_prime[g_index],0, prior_sd_g1[g_index],true)-
    R::dnorm(gamma1[g_index],0, prior_sd_g1[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = gamma1_prime;
    accept=1;
  } else {
    output = gamma1;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}

// [[Rcpp::export]]
List mh_gamma2(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_g2, NumericVector prior_sd_g2){
  
  
  NumericVector gamma2_prime = clone(gamma2);
  
  gamma2_prime[g_index] = R::rnorm(gamma2[g_index], sigma_g2[g_index]);
  
  double ll_diff = ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2_prime,  phi,  Y,  delta_gammas)[1]-
    ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1];
  
  double prior_diff = R::dnorm(gamma2_prime[g_index],0, prior_sd_g2[g_index],true)-
    R::dnorm(gamma2[g_index],0, prior_sd_g2[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = gamma2_prime;
    accept=1;
  } else {
    output = gamma2;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}


// [[Rcpp::export]]
List mh_gamma2_ga_pois(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_g2, NumericVector prior_sd_g2){
  
  
  NumericVector gamma2_prime = clone(gamma2);
  
  gamma2_prime[g_index] = R::rnorm(gamma2[g_index], sigma_g2[g_index]);
  
  double ll_diff = ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2_prime,  phi,  Y,  delta_gammas)[1]-
    ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1];
  
  double prior_diff = R::dnorm(gamma2_prime[g_index],0, prior_sd_g2[g_index],true)-
    R::dnorm(gamma2[g_index],0, prior_sd_g2[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = gamma2_prime;
    accept=1;
  } else {
    output = gamma2;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}


// [[Rcpp::export]]
List mh_gamma2_geoln(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_g2, NumericVector prior_sd_g2){
  
  
  NumericVector gamma2_prime = clone(gamma2);
  
  gamma2_prime[g_index] = R::rnorm(gamma2[g_index], sigma_g2[g_index]);
  
  double ll_diff = ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2_prime,  phi,  Y,  delta_gammas)[1]-
    ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1];
  
  double prior_diff = R::dnorm(gamma2_prime[g_index],0, prior_sd_g2[g_index],true)-
    R::dnorm(gamma2[g_index],0, prior_sd_g2[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = gamma2_prime;
    accept=1;
  } else {
    output = gamma2;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}


// [[Rcpp::export]]
List mh_gamma2_geogamma(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_g2, NumericVector prior_sd_g2){
  
  
  NumericVector gamma2_prime = clone(gamma2);
  
  gamma2_prime[g_index] = R::rnorm(gamma2[g_index], sigma_g2[g_index]);
  
  double ll_diff = ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2_prime,  phi,  Y,  delta_gammas)[1]-
    ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1];
  
  double prior_diff = R::dnorm(gamma2_prime[g_index],0, prior_sd_g2[g_index],true)-
    R::dnorm(gamma2[g_index],0, prior_sd_g2[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = gamma2_prime;
    accept=1;
  } else {
    output = gamma2;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}

// [[Rcpp::export]]
List mh_beta2(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_b2, NumericVector prior_sd_b2){
  
  
  NumericVector beta2_prime = clone(beta2);
  
  beta2_prime[g_index] = R::rnorm(beta2[g_index], sigma_b2[g_index]);
  
  double ll_diff = ll_Poisson_Gamma( beta1,  beta2_prime,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1]-
    ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1];
  
  double prior_diff = R::dnorm(beta2_prime[g_index],0, prior_sd_b2[g_index],true)-
    R::dnorm(beta2[g_index],0, prior_sd_b2[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = beta2_prime;
    accept=1;
  } else {
    output = beta2;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}


// [[Rcpp::export]]
List mh_beta2_ga_pois(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_b2, NumericVector prior_sd_b2){
  
  
  NumericVector beta2_prime = clone(beta2);
  
  beta2_prime[g_index] = R::rnorm(beta2[g_index], sigma_b2[g_index]);
  
  double ll_diff = ll_Gamma_Poisson( beta1,  beta2_prime,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1]-
    ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1];
  
  double prior_diff = R::dnorm(beta2_prime[g_index],0, prior_sd_b2[g_index],true)-
    R::dnorm(beta2[g_index],0, prior_sd_b2[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = beta2_prime;
    accept=1;
  } else {
    output = beta2;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}

// [[Rcpp::export]]
List mh_beta2_geoln(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_b2, NumericVector prior_sd_b2){
  
  
  NumericVector beta2_prime = clone(beta2);
  
  beta2_prime[g_index] = R::rnorm(beta2[g_index], sigma_b2[g_index]);
  
  double ll_diff = ll_Geo_LN( beta1,  beta2_prime,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1]-
    ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1];
  
  double prior_diff = R::dnorm(beta2_prime[g_index],0, prior_sd_b2[g_index],true)-
    R::dnorm(beta2[g_index],0, prior_sd_b2[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = beta2_prime;
    accept=1;
  } else {
    output = beta2;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}


// [[Rcpp::export]]
List mh_beta2_geogamma(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_b2, NumericVector prior_sd_b2){
  
  
  NumericVector beta2_prime = clone(beta2);
  
  beta2_prime[g_index] = R::rnorm(beta2[g_index], sigma_b2[g_index]);
  
  double ll_diff = ll_Geo_Gamma( beta1,  beta2_prime,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1]-
    ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1];
  
  double prior_diff = R::dnorm(beta2_prime[g_index],0, prior_sd_b2[g_index],true)-
    R::dnorm(beta2[g_index],0, prior_sd_b2[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = beta2_prime;
    accept=1;
  } else {
    output = beta2;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}



// [[Rcpp::export]]
List mh_beta1(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_b1, NumericVector prior_sd_b1){
  
  
  NumericVector beta1_prime = clone(beta1);
  
  beta1_prime[g_index] = R::rnorm(beta1[g_index], sigma_b1[g_index]);
  
  double ll_diff = ll_Poisson_Gamma( beta1_prime,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0]-
    ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0];
  
  double prior_diff = R::dnorm(beta1_prime[g_index],0, prior_sd_b1[g_index],true)-
    R::dnorm(beta1[g_index],0, prior_sd_b1[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = beta1_prime;
    accept=1;
  } else {
    output = beta1;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}

// [[Rcpp::export]]
List mh_beta1_ga_pois(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_b1, NumericVector prior_sd_b1){
  
  
  NumericVector beta1_prime = clone(beta1);
  
  beta1_prime[g_index] = R::rnorm(beta1[g_index], sigma_b1[g_index]);
  
  double ll_diff = ll_Gamma_Poisson( beta1_prime,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0]-
    ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0];
  
  double prior_diff = R::dnorm(beta1_prime[g_index],0, prior_sd_b1[g_index],true)-
    R::dnorm(beta1[g_index],0, prior_sd_b1[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = beta1_prime;
    accept=1;
  } else {
    output = beta1;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}


// [[Rcpp::export]]
List mh_beta1_geoln(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_b1, NumericVector prior_sd_b1){
  
  
  NumericVector beta1_prime = clone(beta1);
  
  beta1_prime[g_index] = R::rnorm(beta1[g_index], sigma_b1[g_index]);
  
  double ll_diff = ll_Geo_LN( beta1_prime,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0]-
    ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0];
  
  double prior_diff = R::dnorm(beta1_prime[g_index],0, prior_sd_b1[g_index],true)-
    R::dnorm(beta1[g_index],0, prior_sd_b1[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = beta1_prime;
    accept=1;
  } else {
    output = beta1;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}

// [[Rcpp::export]]
List mh_beta1_geogamma(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, NumericVector sigma_b1, NumericVector prior_sd_b1){
  
  
  NumericVector beta1_prime = clone(beta1);
  
  beta1_prime[g_index] = R::rnorm(beta1[g_index], sigma_b1[g_index]);
  
  double ll_diff = ll_Geo_Gamma( beta1_prime,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0]-
    ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0];
  
  double prior_diff = R::dnorm(beta1_prime[g_index],0, prior_sd_b1[g_index],true)-
    R::dnorm(beta1[g_index],0, prior_sd_b1[g_index],true);
  
  double alpha = ll_diff+prior_diff;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = beta1_prime;
    accept=1;
  } else {
    output = beta1;
    accept = 0;
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept);
  
}

// [[Rcpp::export]]
int delta_g1_cpp(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double omega_g1){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[0];
  NumericVector delta_0_aux = delta_0[0];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  
  double ll_diff = ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_0)[0]-
    ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_1)[0];
  
  
  
  double log_prob_dj = log((1-omega_g1)) - log(omega_g1) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
    
  int draw_j =R::rbinom(1, prob_dj);
    
    return draw_j;
    

}


// [[Rcpp::export]]
int delta_g1_cpp_ga_pois(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double omega_g1){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[0];
  NumericVector delta_0_aux = delta_0[0];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  
  double ll_diff = ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_0)[0]-
    ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_1)[0];
  
  
  
  double log_prob_dj = log((1-omega_g1)) - log(omega_g1) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
  
  int draw_j =R::rbinom(1, prob_dj);
  
  return draw_j;
  
  
}


// [[Rcpp::export]]
int delta_g1_cpp_geoln(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double omega_g1){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[0];
  NumericVector delta_0_aux = delta_0[0];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  
  double ll_diff = ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_0)[0]-
    ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_1)[0];
  
  
  
  double log_prob_dj = log((1-omega_g1)) - log(omega_g1) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
  
  int draw_j =R::rbinom(1, prob_dj);
  
  return draw_j;
  
  
}


// [[Rcpp::export]]
int delta_g1_cpp_geogeo(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, double omega_g1){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[0];
  NumericVector delta_0_aux = delta_0[0];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  
  double ll_diff = ll_Geo_Geo( beta1,  beta2,  gamma1,  gamma2,  Y,  delta_0)[0]-
    ll_Geo_Geo( beta1,  beta2,  gamma1,  gamma2,  Y,  delta_1)[0];
  
  
  
  double log_prob_dj = log((1-omega_g1)) - log(omega_g1) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
  
  int draw_j =R::rbinom(1, prob_dj);
  
  return draw_j;
  
  
}

// [[Rcpp::export]]
int delta_g1_cpp_poispois(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, double omega_g1){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[0];
  NumericVector delta_0_aux = delta_0[0];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  
  double ll_diff = ll_Pois_Pois( beta1,  beta2,  gamma1,  gamma2,  Y,  delta_0)[0]-
    ll_Pois_Pois( beta1,  beta2,  gamma1,  gamma2,  Y,  delta_1)[0];
  
  double log_prob_dj = log((1-omega_g1)) - log(omega_g1) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
  
  int draw_j =R::rbinom(1, prob_dj);
  
  return draw_j;
  
  
}


// [[Rcpp::export]]
int delta_g1_cpp_geogamma(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double omega_g1){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[0];
  NumericVector delta_0_aux = delta_0[0];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  
  double ll_diff = ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_0)[0]-
    ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_1)[0];
  
  
  
  double log_prob_dj = log((1-omega_g1)) - log(omega_g1) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
  
  int draw_j =R::rbinom(1, prob_dj);
  
  return draw_j;
  
  
}

// [[Rcpp::export]]
int delta_g2_cpp(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double omega_g2){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[1];
  NumericVector delta_0_aux = delta_0[1];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  double ll_diff = ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_0)[1]-
    ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_1)[1];
  
  double log_prob_dj = log((1-omega_g2)) - log(omega_g2) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
  
  int draw_j =R::rbinom(1, prob_dj);
  
  return draw_j;
  
}


// [[Rcpp::export]]
int delta_g2_cpp_ga_pois(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double omega_g2){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[1];
  NumericVector delta_0_aux = delta_0[1];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  double ll_diff = ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_0)[1]-
    ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_1)[1];
  
  double log_prob_dj = log((1-omega_g2)) - log(omega_g2) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
  
  int draw_j =R::rbinom(1, prob_dj);
  
  return draw_j;
  
}


// [[Rcpp::export]]
int delta_g2_cpp_geoln(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double omega_g2){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[1];
  NumericVector delta_0_aux = delta_0[1];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  double ll_diff = ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_0)[1]-
    ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_1)[1];
  
  double log_prob_dj = log((1-omega_g2)) - log(omega_g2) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
  
  int draw_j =R::rbinom(1, prob_dj);
  
  return draw_j;
  
}


// [[Rcpp::export]]
int delta_g2_cpp_geogeo(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, double omega_g2){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[1];
  NumericVector delta_0_aux = delta_0[1];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  double ll_diff = ll_Geo_Geo( beta1,  beta2,  gamma1,  gamma2,  Y,  delta_0)[1]-
    ll_Geo_Geo( beta1,  beta2,  gamma1,  gamma2,  Y,  delta_1)[1];
  
  double log_prob_dj = log((1-omega_g2)) - log(omega_g2) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
  
  int draw_j =R::rbinom(1, prob_dj);
  
  return draw_j;
  
}


// [[Rcpp::export]]
int delta_g2_cpp_poispois(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, double omega_g2){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[1];
  NumericVector delta_0_aux = delta_0[1];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  double ll_diff = ll_Pois_Pois( beta1,  beta2,  gamma1,  gamma2,  Y,  delta_0)[1]-
    ll_Pois_Pois( beta1,  beta2,  gamma1,  gamma2,  Y,  delta_1)[1];
  
  double log_prob_dj = log((1-omega_g2)) - log(omega_g2) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
  
  int draw_j =R::rbinom(1, prob_dj);
  
  return draw_j;
  
}


// [[Rcpp::export]]
int delta_g2_cpp_geogamma(int index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double omega_g2){
  
  
  List delta_1 = clone(delta_gammas);
  List delta_0 = clone(delta_gammas);
  
  NumericVector delta_1_aux= delta_1[1];
  NumericVector delta_0_aux = delta_0[1];
  delta_1_aux[index]=1;
  delta_0_aux[index]=0;
  
  double ll_diff = ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_0)[1]-
    ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_1)[1];
  
  double log_prob_dj = log((1-omega_g2)) - log(omega_g2) + ll_diff;
  
  double prob_dj = 1/(1 + exp(log_prob_dj));
  
  int draw_j =R::rbinom(1, prob_dj);
  
  return draw_j;
  
}


NumericVector proposal_ln(double param, double sigma){
  NumericVector param_prime = rlnorm(1, log(param), sigma);
  return(param_prime);
}


// [[Rcpp::export]]
List mh_phi(NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double sigma_phi, double prior_rate_phi){
  
  double phi_prime = proposal_ln(phi, sigma_phi)[0];
  
  double ll_diff = ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi_prime,  Y,  delta_gammas)[0]-
    ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0];
  
  double prior_diff = R::dexp(phi_prime, prior_rate_phi,true)-R::dexp(phi, prior_rate_phi,true);

  double qratio =  log(phi) - log(phi_prime) ;
  
  double alpha = ll_diff+prior_diff+qratio;
  double ll;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = phi_prime;
    accept=1;
    ll = sum(ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi_prime,  Y,  delta_gammas));
  } else {
    output = phi;
    accept = 0;
    ll = sum(ll_Poisson_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas));
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept,
                      _["ll"] = ll);
  
}

// [[Rcpp::export]]
List mh_phi_ga_pois(NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double sigma_phi, double prior_rate_phi){
  
  double phi_prime = proposal_ln(phi, sigma_phi)[0];
  
  double ll_diff = ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi_prime,  Y,  delta_gammas)[1]-
    ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[1];
  
  double prior_diff = R::dexp(phi_prime, prior_rate_phi,true)-R::dexp(phi, prior_rate_phi,true);
  
  double qratio =  log(phi) - log(phi_prime) ;
  
  double alpha = ll_diff+prior_diff+qratio;
  double ll;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = phi_prime;
    accept=1;
    ll = sum(ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi_prime,  Y,  delta_gammas));
  } else {
    output = phi;
    accept = 0;
    ll = sum(ll_Gamma_Poisson( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas));
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept,
                      _["ll"] = ll);
  
}

// [[Rcpp::export]]
List mh_phi_geoln(NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double sigma_phi, double prior_rate_phi){
  
  double phi_prime = proposal_ln(phi, sigma_phi)[0];
  
  double ll_diff = ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi_prime,  Y,  delta_gammas)[0]-
    ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0];
  
  double prior_diff = R::dexp(phi_prime, prior_rate_phi,true)-R::dexp(phi, prior_rate_phi,true);
  
  double qratio =  log(phi) - log(phi_prime) ;
  
  double alpha = ll_diff+prior_diff+qratio;
  double ll;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = phi_prime;
    ll = sum(ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi_prime,  Y,  delta_gammas));
    accept=1;
  } else {
    output = phi;
    accept = 0;
    ll =sum(ll_Geo_LN( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas));
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept,
                      _["ll"] = ll);
  
}


// [[Rcpp::export]]
List mh_phi_geogamma(NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, double phi, NumericMatrix Y, List delta_gammas, double sigma_phi, double prior_rate_phi){
  
  double phi_prime = proposal_ln(phi, sigma_phi)[0];
  
  double ll_diff = ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi_prime,  Y,  delta_gammas)[0]-
    ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas)[0];
  
  double prior_diff = R::dexp(phi_prime, prior_rate_phi,true)-R::dexp(phi, prior_rate_phi,true);
  
  double qratio =  log(phi) - log(phi_prime) ;
  
  double alpha = ll_diff+prior_diff+qratio;
  double ll;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector output;
  
  if( compare[0]==true ){
    output = phi_prime;
    ll = sum(ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi_prime,  Y,  delta_gammas));
    accept=1;
  } else {
    output = phi;
    accept = 0;
    ll =sum(ll_Geo_Gamma( beta1,  beta2,  gamma1,  gamma2,  phi,  Y,  delta_gammas));
  }
  return List::create(_["param"] = output, 
                      _["accepted"] = accept,
                      _["ll"] = ll);
  
}


// [[Rcpp::export]]
List mh_gamma1_geo_geo(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, NumericVector sigma_g1, NumericVector prior_sd_g1){
  NumericVector gamma1_prime = clone(gamma1);
  gamma1_prime[g_index] = R::rnorm(gamma1[g_index], sigma_g1[g_index]);
  
  double ll_diff = ll_Geo_Geo(beta1, beta2, gamma1_prime, gamma2, Y, delta_gammas)[0] -
    ll_Geo_Geo(beta1, beta2, gamma1, gamma2, Y, delta_gammas)[0];
  
  double prior_diff = R::dnorm(gamma1_prime[g_index], 0, prior_sd_g1[g_index], true) -
    R::dnorm(gamma1[g_index], 0, prior_sd_g1[g_index], true);
  
  double alpha = ll_diff + prior_diff;
  LogicalVector compare = log(runif(1)) < alpha;
  
  return List::create(
    _["param"] = compare[0] ? gamma1_prime : gamma1,
    _["accepted"] = compare[0] ? 1 : 0
  );
}

// [[Rcpp::export]]
List mh_gamma2_geo_geo(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, NumericVector sigma_g2, NumericVector prior_sd_g2){
  NumericVector gamma2_prime = clone(gamma2);
  gamma2_prime[g_index] = R::rnorm(gamma2[g_index], sigma_g2[g_index]);
  
  double ll_diff = ll_Geo_Geo(beta1, beta2, gamma1, gamma2_prime, Y, delta_gammas)[1] -
    ll_Geo_Geo(beta1, beta2, gamma1, gamma2, Y, delta_gammas)[1];
  
  double prior_diff = R::dnorm(gamma2_prime[g_index], 0, prior_sd_g2[g_index], true) -
    R::dnorm(gamma2[g_index], 0, prior_sd_g2[g_index], true);
  
  double alpha = ll_diff + prior_diff;
  LogicalVector compare = log(runif(1)) < alpha;
  
  return List::create(
    _["param"] = compare[0] ? gamma2_prime : gamma2,
    _["accepted"] = compare[0] ? 1 : 0
  );
}

// [[Rcpp::export]]
List mh_beta1_geo_geo(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, NumericVector sigma_b1, NumericVector prior_sd_b1){
  NumericVector beta1_prime = clone(beta1);
  beta1_prime[g_index] = R::rnorm(beta1[g_index], sigma_b1[g_index]);
  
  double ll_diff = ll_Geo_Geo(beta1_prime, beta2, gamma1, gamma2, Y, delta_gammas)[0] -
    ll_Geo_Geo(beta1, beta2, gamma1, gamma2, Y, delta_gammas)[0];
  
  double prior_diff = R::dnorm(beta1_prime[g_index], 0, prior_sd_b1[g_index], true) -
    R::dnorm(beta1[g_index], 0, prior_sd_b1[g_index], true);
  
  double alpha = ll_diff + prior_diff;
  LogicalVector compare = log(runif(1)) < alpha;
  
  return List::create(
    _["param"] = compare[0] ? beta1_prime : beta1,
    _["accepted"] = compare[0] ? 1 : 0
  );
}

// [[Rcpp::export]]
List mh_beta2_geo_geo(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, NumericVector sigma_b2, NumericVector prior_sd_b2){
  NumericVector beta2_prime = clone(beta2);
  beta2_prime[g_index] = R::rnorm(beta2[g_index], sigma_b2[g_index]);
  
  double ll_diff = ll_Geo_Geo(beta1, beta2_prime, gamma1, gamma2, Y, delta_gammas)[1] -
    ll_Geo_Geo(beta1, beta2, gamma1, gamma2, Y, delta_gammas)[1];
  
  double prior_diff = R::dnorm(beta2_prime[g_index], 0, prior_sd_b2[g_index], true) -
    R::dnorm(beta2[g_index], 0, prior_sd_b2[g_index], true);
  
  double alpha = ll_diff + prior_diff;
  LogicalVector compare = log(runif(1)) < alpha;
  
  return List::create(
    _["param"] = compare[0] ? beta2_prime : beta2,
    _["accepted"] = compare[0] ? 1 : 0
  );
}


// [[Rcpp::export]]
List mh_gamma1_pois_pois(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, NumericVector sigma_g1, NumericVector prior_sd_g1){
  NumericVector gamma1_prime = clone(gamma1);
  gamma1_prime[g_index] = R::rnorm(gamma1[g_index], sigma_g1[g_index]);
  
  double ll_diff = ll_Pois_Pois(beta1, beta2, gamma1_prime, gamma2, Y, delta_gammas)[0] -
    ll_Pois_Pois(beta1, beta2, gamma1, gamma2, Y, delta_gammas)[0];
  
  double prior_diff = R::dnorm(gamma1_prime[g_index], 0, prior_sd_g1[g_index], true) -
    R::dnorm(gamma1[g_index], 0, prior_sd_g1[g_index], true);
  
  double alpha = ll_diff + prior_diff;
  LogicalVector compare = log(runif(1)) < alpha;
  
  return List::create(
    _["param"] = compare[0] ? gamma1_prime : gamma1,
    _["accepted"] = compare[0] ? 1 : 0
  );
}

// [[Rcpp::export]]
List mh_gamma2_pois_pois(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, NumericVector sigma_g2, NumericVector prior_sd_g2){
  NumericVector gamma2_prime = clone(gamma2);
  gamma2_prime[g_index] = R::rnorm(gamma2[g_index], sigma_g2[g_index]);
  
  double ll_diff = ll_Pois_Pois(beta1, beta2, gamma1, gamma2_prime, Y, delta_gammas)[1] -
    ll_Pois_Pois(beta1, beta2, gamma1, gamma2, Y, delta_gammas)[1];
  
  double prior_diff = R::dnorm(gamma2_prime[g_index], 0, prior_sd_g2[g_index], true) -
    R::dnorm(gamma2[g_index], 0, prior_sd_g2[g_index], true);
  
  double alpha = ll_diff + prior_diff;
  LogicalVector compare = log(runif(1)) < alpha;
  
  return List::create(
    _["param"] = compare[0] ? gamma2_prime : gamma2,
    _["accepted"] = compare[0] ? 1 : 0
  );
}

// [[Rcpp::export]]
List mh_beta1_pois_pois(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, NumericVector sigma_b1, NumericVector prior_sd_b1){
  NumericVector beta1_prime = clone(beta1);
  beta1_prime[g_index] = R::rnorm(beta1[g_index], sigma_b1[g_index]);
  
  double ll_diff = ll_Pois_Pois(beta1_prime, beta2, gamma1, gamma2, Y, delta_gammas)[0] -
    ll_Pois_Pois(beta1, beta2, gamma1, gamma2, Y, delta_gammas)[0];
  
  double prior_diff = R::dnorm(beta1_prime[g_index], 0, prior_sd_b1[g_index], true) -
    R::dnorm(beta1[g_index], 0, prior_sd_b1[g_index], true);
  
  double alpha = ll_diff + prior_diff;
  LogicalVector compare = log(runif(1)) < alpha;
  
  return List::create(
    _["param"] = compare[0] ? beta1_prime : beta1,
    _["accepted"] = compare[0] ? 1 : 0
  );
}

// [[Rcpp::export]]
List mh_beta2_pois_pois(int g_index, NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, NumericMatrix Y, List delta_gammas, NumericVector sigma_b2, NumericVector prior_sd_b2){
  NumericVector beta2_prime = clone(beta2);
  beta2_prime[g_index] = R::rnorm(beta2[g_index], sigma_b2[g_index]);
  
  double ll_diff = ll_Pois_Pois(beta1, beta2_prime, gamma1, gamma2, Y, delta_gammas)[1] -
    ll_Pois_Pois(beta1, beta2, gamma1, gamma2, Y, delta_gammas)[1];
  
  double prior_diff = R::dnorm(beta2_prime[g_index], 0, prior_sd_b2[g_index], true) -
    R::dnorm(beta2[g_index], 0, prior_sd_b2[g_index], true);
  
  double alpha = ll_diff + prior_diff;
  LogicalVector compare = log(runif(1)) < alpha;
  
  return List::create(
    _["param"] = compare[0] ? beta2_prime : beta2,
    _["accepted"] = compare[0] ? 1 : 0
  );
}
