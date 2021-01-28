#include <RcppArmadillo.h>
#include <math.h>
#include <time.h>
#include <mvnorm.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <bits/stdc++.h> 
  
using namespace Rcpp;

// log-likelihood
/*
 * y = (y_1, ... , y_n)  vec length = n
 * cc = (c_0, ..., c_n-1, c_n)  vec length = n + 1
 * AA = (A_0 = 0, A_1, ... , A_n)
 */
double loglik(const arma::vec& y, const arma::vec& cc, const arma::vec& AA, 
              double & b, double & gamma, double & sigma2, double & tau2)
{
  int n = y.n_elem;
  arma::vec llik(n);
  
  for(int k = 0; k < n; k++) { 
    llik(k) = R::dnorm(y(k), b + gamma * cc(k) + AA(k+1), std::sqrt(sigma2 + tau2), true) ;
  }
  return arma::accu(llik); 
}


// prior on gamma (AR parameter)
/*
 * gamma ~ Beta(hyp_gamma1, hyp_gamma2)
 */
double logprior_gamma(double & gamma, double & hyp_gamma1, double & hyp_gamma2) 
{
  double out;
  out = R::dbeta( gamma, hyp_gamma1, hyp_gamma2, true);
  return(out);
}

// log-posterior (MH step on gamma)
double logpost_gamma(const arma::vec& y, const arma::vec& cc, const arma::vec& AA, 
                     double & b, double & gamma, double & sigma2, double & tau2, 
                     double & hyp_gamma1, double & hyp_gamma2)
{
  double out;
  out = loglik(y, cc, AA, b, gamma, sigma2, tau2) + logprior_gamma(gamma, hyp_gamma1, hyp_gamma2) ;
  return(out);
}


// prior on A for non-zero values
/*
 * A ~ G0 = Gamma(hyp_A1, hyp_A2) 
 * the mean is hyp_A1/hyp_A2
 */
double logprior_A(double & A, double & hyp_A1, double & hyp_A2) 
{
  double out;
  out = R::dgamma( A, hyp_A1, 1/hyp_A2, true );
  return(out);
}

// log-posterior on A
// the function takes as input only those y_i s.t. A_i = A_k for some k = 1,..., K (fixed).
double logpost_A(const arma::vec& y, const arma::vec& cc, double & A, 
                 double & b, double & gamma, double & sigma2, double & tau2, 
                 double & hyp_A1, double & hyp_A2)
{
  double out ;
  arma::vec AAA(y.n_elem + 1) ; AAA.fill(A) ;
  out = loglik(y, cc, AAA, b, gamma, sigma2, tau2) + logprior_A(A, hyp_A1, hyp_A2) ;
  return(out) ;
}

// sampling from prior on A (spike and slab)
double sample_mix(double & p, double & hyp_A1, double & hyp_A2)
{
  double out = 0 ;
  int k = Rcpp::rbinom( 1, 1, p )[0] ;
  if( k == 1 ) {
    out = R::rgamma(hyp_A1, 1/hyp_A2) ;
  }
  return(out) ;
}



/* 
 * functions for sampling from gMFM
 */

// sample from dirichlet distribution
arma::vec rdirichlet(arma::vec alpha_m) 
{
  int distribution_size = alpha_m.n_elem ;
  arma::vec distribution = arma::zeros(distribution_size) ;
  
  double sum_term = 0 ;
  // draw Gamma variables
  for (int j = 0; j < distribution_size; j++) 
  {
    double cur = R::rgamma(alpha_m[j], 1.0) ;
    distribution(j) = cur ;
    sum_term += cur ;
  }
  // normalize
  for (int j = 0; j < distribution_size; j++) 
  {
    distribution(j) = distribution(j)/sum_term ;
  }
  return(distribution) ;
}


// product of gamma functions in the likelihood (step 3 in Algorithm 1 of Fruehwirth-Schnatter et al. 2020)
double log_prod_frac(int maxx,
                 const arma::vec cluster_size,
                 double beta,
                 int maxlab)
{
  double out ;
  arma::vec tmp(maxlab) ;
  for(int h = 0; h < maxlab ; h ++)
  {
    tmp(h) = lgamma( cluster_size(h) + beta/maxx ) - lgamma( 1 + beta/maxx ) ;
  }
  out = arma::accu(tmp) ;
  return(out) ;
}


// prior on Dirichlet hyperparameter Dir(beta/K)
double logprior_beta(double beta, 
                     double hyp_beta1, double hyp_beta2)
{
  double out ;
  out = R::dgamma(beta, hyp_beta1, 1/hyp_beta2, true) ;
  //out = R::df(beta, hyp_beta1, hyp_beta2, true) ;
  return(out) ;
}

// log-posterior for MH step on beta
double logpost_beta(double beta, 
                    double hyp_beta1, double hyp_beta2,
                    int T,
                    const arma::vec clusterO_size,
                    int maxlabel,
                    int maxL)
{
  double out ;
  out = logprior_beta(beta, hyp_beta1, hyp_beta2) + maxlabel * log(beta) + lgamma(beta) - 
    lgamma(T + beta) + log_prod_frac(maxL, clusterO_size, beta, maxlabel) ;
  return(out) ;
}


// prior on maxL
/*
 * max - 1 ~ BNB(hyp1, hyp2, hyp3)
 */
double lbeta(double a, double b)
{
  double out = 0;
  out = lgamma(a) + lgamma(b) - lgamma(a+b) ;
  return(out) ;
}

double logprior_maxx(int maxx, double hyp_max1,
                      double hyp_max2, double hyp_max3)
{
  double out ;
  out = lgamma( hyp_max1 + maxx - 1 ) + lbeta( hyp_max1 + hyp_max2, maxx - 1 + hyp_max3) -
    lgamma( hyp_max1 ) - lgamma( maxx ) - lbeta( hyp_max2, hyp_max3 ) ;
  return(out) ;
}

// log-posterior for MH step on maxL
double logpost_maxx(int maxx, 
                    double hyp_max1, 
                    double hyp_max2, 
                    double hyp_max3, 
                    double d_par,
                    const arma::vec cluster_size,
                    int max_lab)
{
  double out ;
  out = logprior_maxx(maxx, hyp_max1, hyp_max2, hyp_max3) + max_lab * log(d_par) + lgamma(maxx + 1) - max_lab * log(maxx) - 
    lgamma(maxx - max_lab + 1) + log_prod_frac(maxx, cluster_size, d_par, max_lab) ;
  return(out) ;
}


// Sampling of generalized mixtures of finite mixtures
/*
 * nested telescoping sampler (see Appendix)
 */
Rcpp::List nested_gMFM(const arma::vec& y, const arma::vec& g, 
                          arma::vec clusterD, arma::vec clusterO, 
                          const arma::vec& cc,
                          arma::vec A, double b, double gamma, 
                          double p,
                          double sigma2, double tau2,
                          double & alpha, double & beta, 
                          arma::vec pi_k,
                          arma::mat omega_lk,
                          int maxK, int maxL,
                          double & hyp_alpha1, double & hyp_alpha2, double & eps_alpha, 
                          double & hyp_beta1, double & hyp_beta2, double & eps_beta, 
                          double & hyp_A1, double & hyp_A2,  
                          double & hyp_maxK1, double & hyp_maxK2, double & hyp_maxK3, 
                          double & hyp_maxL1, double & hyp_maxL2, double & hyp_maxL3, 
                          double & eps_A,
                          int & eps_maxK, int & eps_maxL,
                          int & check)
{
  /*
   * clusterD (length = J) cluster allocation of the distributions, with values in {1,2,3,...}
   * clusterO (length = T) cluster allocation of the observations, with values in {0,1,2,3,...}, where the 0 is the absence of a spike
   */
 
  int T = y.n_elem ;
  arma::vec clusterD_long(T) ; 
  for(int t = 0; t < T; t++) { clusterD_long(t) = clusterD( g(t)-1 ) ; }
  
  arma::vec unique_g = arma::unique(g) ;
  int J = unique_g.n_elem ;
  
  double oldA ;
  double newA ;
  double oldalpha ;
  double newalpha ;
  double oldbeta ;
  double newbeta ;
  int oldmaxK ;
  int newmaxK ;
  int oldmaxL ;
  int newmaxL ;
  double ratio ;
 

  // step 1: sample the weights on the distributions
  arma::vec dir_paramD(maxK) ;
  for(int k = 1; k < maxK + 1; k++)
  {
    arma::uvec ind_k = find(clusterD == k) ;
    dir_paramD(k-1) = alpha/maxK + ind_k.n_elem ;
  }
  pi_k.fill(0) ;
  arma::vec sample_dir_d = rdirichlet(dir_paramD) ;
  for(int k = 0; k < maxK; k++)
  {
    pi_k(k) = sample_dir_d(k) ;
  }
  
 
  // step 2: sample the weights on the observations
  // for each k in {1,...,maxK} we have maxL weights: matrix(maxL, maxK)
  for(int k = 1; k < J*3 + 1; k++)
  {
    omega_lk.col(k-1).fill(0) ;
    arma::uvec ind_clusterD_k = find(clusterD_long == k) ; // indices of the y_t s.t. clusterD_t = k
    arma::vec subcluster = clusterO.elem( ind_clusterD_k ) ; // cluster labels for those y ~ G*k
    arma::vec dir_param(maxL) ;
    for(int l = 0; l < maxL ; l++)
    {
      arma::uvec subcluster_l = find( subcluster == l ) ;
      dir_param(l) = beta/(maxL+1) + subcluster_l.n_elem ;
    }
    arma::vec sample_dir = rdirichlet(dir_param) ;
    for(int l = 0; l < maxL ; l++)
    {
      omega_lk(l, k-1) = sample_dir(l) ;
    }
  }
  
  // step 3: sample the distributional cluster indicator
  NumericVector probK(maxK) ;
  IntegerVector clusterD_id =  Rcpp::seq(1, maxK);
  for(int j = 0; j < J; j++)
  {
    arma::uvec ind_t = find(g == (j+1)) ;
    arma::vec mixdens(ind_t.n_elem) ;
    
    for(int k = 0; k < maxK; k++)
    {
      probK[k] = -9999 ;
      arma::vec omega_col = omega_lk.col(k) ;
      for(int t = 0; t < ind_t.n_elem; t++)
      {
        int cl = clusterO( ind_t(t) ) ;
        mixdens(t) = log( omega_col(cl) )  + 
          R::dnorm(y(ind_t(t)) - b - gamma * cc(ind_t(t)) - A(cl), 0, std::sqrt(sigma2 + tau2), true) ;
      }
      
      probK[k] =  log( pi_k(k) ) + arma::accu(mixdens) ;
    }
    probK =  probK - max(probK)  ;
    probK =  exp(probK)  ;
    clusterD(j) = Rcpp::sample(clusterD_id, 1, false, probK )[0] ;
  } 
  for(int t = 0; t < T; t++) { clusterD_long(t) = clusterD( g(t) - 1 ) ; }


  // step 3b: relabel the clusters so that the first k = maxdist are non-empty
  // determine the clusters size and the maximum occupied cluster (maxdist)
  arma::vec clusterD_size(maxK) ;
  for(int k = 0; k < maxK; k++)
  {
    arma::uvec ind = find( clusterD == k + 1 ) ;   
    clusterD_size(k) = ind.n_elem ;
  }

  int cumsumD = 0 ;
  for(int k = 0; k < maxK ; k++)
  {
    int kk = k ;
    cumsumD = cumsumD + clusterD_size(k) ;
    if( (clusterD_size(k) == 0) & ( cumsumD < J ) )
    {
      do{ kk++ ; } while (clusterD_size(kk) == 0);
      arma::uvec ind = find( clusterD == kk + 1 ) ;
      clusterD.elem(ind).fill(k + 1) ;
      
      clusterD_size(k) = clusterD_size(kk) ;
      clusterD_size(kk) = 0 ;
      
      cumsumD = cumsumD + clusterD_size(k) ;
    }
    if(cumsumD == J) { k = maxK - 1 ; }
  }  
  int maxdistr = max(clusterD) ; 

  // step 4: sample the observational cluster indicator
  NumericVector probL(maxL) ;
  IntegerVector clusterO_id = Rcpp::seq(0, maxL-1) ;
  for(int t = 0; t < T; t++)
  {
    for(int l = 0; l < maxL; l++)
    {
      probL[l] = log( omega_lk( l, clusterD_long(t) - 1 ) ) + 
        R::dnorm(y(t) - b - gamma * cc(t), A(l), std::sqrt(sigma2 + tau2), true) ;
    }
    probL =  probL - max(probL) ;
    probL = exp(probL) ;
    clusterO(t) = Rcpp::sample( clusterO_id, 1, false, probL )[0] ;
    if( A(clusterO(t)) == 0 ) { clusterO(t) = 0 ; }
  } 
  
  // step 4b: relabel the clusters so that the first l=maxlabel are non-empty
  // determine the clusters size and the maximum occupied cluster (maxlabel)
  arma::vec clusterO_size(maxL) ;
  for(int l = 0; l < maxL; l++)
  {
    arma::uvec ind = find( clusterO == l ) ;
    clusterO_size(l) = ind.n_elem ;
  }
  
  int cumsum = 0 ;
  for(int l = 0; l < maxL ; l++)
  {
    int ll = l ;
    cumsum = cumsum + clusterO_size(l) ;
    if( (clusterO_size(l) == 0) & ( cumsum < T ) )
    {
      do{ ll++ ; } while (clusterO_size(ll) == 0);
      arma::uvec ind = find( clusterO == ll ) ;
      clusterO.elem(ind).fill(l) ;
      
      clusterO_size(l) = clusterO_size(ll) ;
      clusterO_size(ll) = 0 ;
      
      A(l) = A(ll) ;
      A(ll) = 0 ;
      
      cumsum = cumsum + clusterO_size(l) ;
    }
    if(cumsum == T) { l = maxL - 1 ; }
  }
  int maxlabel = max(clusterO) ; 

  // step 5: sample the cluster parameters for the non-empty clusters
  for(int l = 1; l < maxlabel + 1; l++)
  {
    arma::uvec ind_l = find( clusterO == l ) ;
    arma::vec sub_y = y(ind_l) ;
    // MH step on A(l)
    oldA = A(l) ;
    newA = oldA ;
    
    newA = oldA + R::runif(-eps_A, eps_A) ;
    ratio = exp( logpost_A(sub_y, cc(ind_l), 
                           newA, b, 
                           gamma, sigma2, tau2,
                           hyp_A1, hyp_A2) -
                  logpost_A(sub_y, cc(ind_l), 
                            oldA, b, 
                            gamma, sigma2, tau2,
                            hyp_A1, hyp_A2) ) ;
    if(R::runif(0, 1) < ratio) oldA = newA ;
    A(l) = oldA ; 
  }

  // step 6: sample the number of components on the distributions
  oldmaxK = maxK ;
  newmaxK = maxK ;
  NumericVector prob_maxK(eps_maxK) ;
  prob_maxK.fill(1.0/eps_maxK) ;
  IntegerVector maxK_id =  Rcpp::seq(maxdistr, maxdistr + eps_maxK - 1);
  newmaxK =  Rcpp::sample(maxK_id, 1, false, prob_maxK)[0] ;
  if(newmaxK > 3 * J) {newmaxK = 3 * J ;}
  ratio = exp( logpost_maxx(newmaxK, 
                             hyp_maxK1, hyp_maxK2, hyp_maxK3,
                             alpha,
                             clusterD_size,
                             maxdistr) -
                logpost_maxx(oldmaxK, 
                             hyp_maxK1, hyp_maxK2, hyp_maxK3,
                             alpha,
                             clusterD_size,
                             maxdistr) ) ;
  if(R::runif(0, 1) < ratio) oldmaxK = newmaxK ;
  if(oldmaxK == 1) oldmaxK = maxK ;
  maxK = oldmaxK ;  

  // step 7: sample the number of components on the observations
  oldmaxL = maxL ;
  newmaxL = maxL ;
  
  NumericVector prob_maxL(eps_maxL) ;
  prob_maxL.fill(1.0/eps_maxL) ;
  IntegerVector maxL_id =  Rcpp::seq(maxlabel + 1, maxlabel + eps_maxL );
  
  newmaxL =  Rcpp::sample(maxL_id, 1, false, prob_maxL)[0] ;
  ratio = exp( logpost_maxx(newmaxL + 1, 
                            hyp_maxL1, hyp_maxL2, hyp_maxL3,
                            beta,
                            clusterO_size,
                            maxlabel + 1) -
                logpost_maxx(oldmaxL + 1, 
                            hyp_maxL1, hyp_maxL2, hyp_maxL3,
                            beta,
                            clusterO_size,
                            maxlabel + 1) ) ;
  if(R::runif(0, 1) < ratio) oldmaxL = newmaxL ;
  int uu = A.n_elem ;
  maxL = std::min(oldmaxL, uu) ; 

  // step 8 update hyperparameter on Dirichlet distr on the distributions (MH step)
  oldalpha = alpha ;
  newalpha = alpha ;
  newalpha = exp( R::rnorm( log(oldalpha), eps_alpha ) ) ;
  ratio = exp( logpost_beta(newalpha, 
                            hyp_alpha1, hyp_alpha2,
                            T,
                            clusterD_size,
                            maxdistr,
                            maxK) -
                logpost_beta(oldalpha, 
                            hyp_alpha1, hyp_alpha2,
                            T,
                            clusterD_size,
                            maxdistr,
                            maxK) ) ;
  if(R::runif(0, 1) < ratio) oldalpha = newalpha ;
  alpha = oldalpha ;

  // step 9: update hyperparameter on Dirichlet distr on the observations (MH step)
  oldbeta = beta ;
  newbeta = beta ;
  newbeta = exp( R::rnorm( log(oldbeta), eps_beta ) ) ;
  ratio = exp( logpost_beta(newbeta, 
                            hyp_beta1, hyp_beta2,
                            T,
                            clusterO_size,
                            maxlabel + 1,
                            maxL + 1) -
                logpost_beta(oldbeta, 
                             hyp_beta1, hyp_beta2,
                             T,
                             clusterO_size,
                             maxlabel + 1,
                             maxL + 1) ) ;
  if(R::runif(0, 1) < ratio) oldbeta = newbeta ;
  beta = oldbeta ; 
  
  // step 10: sample a cluster parameter for the empty components
  int empty_comp = maxL - maxlabel ;
  if(empty_comp > 0)
  {
    for(int h = maxlabel + 1; h < maxL; h++)
    {
      A(h) = sample_mix(p, hyp_A1, hyp_A2) ;
    }
  }

  return Rcpp::List::create(Rcpp::Named("clusterO") = clusterO,
                            Rcpp::Named("clusterD") = clusterD,
                            Rcpp::Named("A") = A,
                            Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("pi_k") = pi_k,
                            Rcpp::Named("omega_lk") = omega_lk,
                            Rcpp::Named("maxK") = maxK,
                            Rcpp::Named("maxL") = maxL);
}



// Gibbs sampler
/*
 * this is the main function
 */
// [[Rcpp::export]]
Rcpp::List calcium_gibbs(int Nrep, 
                         arma::vec y, arma::vec g,
                         arma::vec cal,
                         arma::vec clO, arma::vec clD, 
                         arma::vec A_start,
                         double b_start, double gamma_start, 
                         double sigma2_start, 
                         double tau2_start, // calcium variance
                         double p_start,
                         double alpha_start, // Dirichlet param on distributions
                         double beta_start, // Dirichlet param on observations
                         int maxK_start, int maxL_start,
                         double c0, double varC0, // mean and variance of c0 (calcium at time 0)
                         double hyp_A1, double hyp_A2, // shape and rate parameters Gamma distr on A
                         double hyp_b1, double hyp_b2, // prior mean and variance Normal distr on b
                         double hyp_sigma21, double hyp_sigma22, // shape and rate parameters Gamma distr on 1/sigma2 (precision)
                         double hyp_tau21, double hyp_tau22, // shape and rate parameters Gamma distr on 1/tau2 (precision state equation)
                         double hyp_gamma1, double hyp_gamma2, // shape1 and shape2 parameters Beta distr on gamma (AR parameter)
                         double hyp_p1, double hyp_p2, // shape1 and shape2 parameters Beta distr on p (spike probability)
                         double hyp_alpha1, double hyp_alpha2,
                         double hyp_beta1, double hyp_beta2, 
                         double hyp_maxK1, double hyp_maxK2, double hyp_maxK3, 
                         double hyp_maxL1, double hyp_maxL2, double hyp_maxL3, 
                         double eps_alpha, double eps_beta, 
                         double eps_gamma, // MH step size
                         double eps_A,
                         int eps_maxK, int eps_maxL) 
{
  // allocate output matrices
  int n = y.n_elem ;
  arma::vec unique_g = arma::unique(g) ;
  int J = unique_g.n_elem ;
  
  arma::mat out_c(n+1, Nrep) ; // 0 1 ... n
  arma::mat out_A = arma::zeros(A_start.n_elem, Nrep) ;
  arma::vec out_b(Nrep) ;
  arma::vec out_gamma = arma::zeros(Nrep) ;
  arma::vec out_sigma2 = arma::zeros(Nrep) ;
  arma::vec out_tau2 = arma::zeros(Nrep) ;
  arma::mat clusterO = arma::zeros(n, Nrep)  ;
  arma::mat clusterD = arma::zeros(J, Nrep)  ;
  arma::vec out_p(Nrep) ;
  arma::vec out_alpha(Nrep) ;
  arma::vec out_beta(Nrep) ;
  arma::vec out_maxK(Nrep) ;
  arma::vec out_maxL(Nrep) ;
  
  // initialize the chains
  out_c(0,0) = c0 ;
  out_b(0) = b_start ;
  out_gamma(0) = gamma_start ;
  out_sigma2(0) = sigma2_start ;
  out_tau2(0) = tau2_start ;
  out_p(0) = p_start ;
  out_alpha(0) = alpha_start ;
  out_beta(0) = beta_start ;
  out_maxK(0) = maxK_start ;
  out_maxL(0) = maxL_start ;
  
  // init kalman filter quantities
  out_c.col(0) = cal;
  arma::vec filter_mean(n+1) ; // 0 1 ... n
  arma::vec filter_var(n+1) ; // 0 1 ... n
  arma::vec R(n+1) ; arma::vec a(n+1) ; // 0 1 ... n
  double back_mean; double back_var ;
  
  // MH quantities
  double oldgamma ; double newgamma ;
  double ratio ;
  
  // CAM quantities
  double n_clus ;
  arma::vec AA = arma::zeros(n+1) ;
  Rcpp::List out_slice ;
  
  arma::mat omega_lk(A_start.n_elem, J*3, arma::fill::zeros) ;
  for(int k = 0; k < J ; k ++)
  {
    for(int l = 0; l < maxL_start; l++)
    {
      omega_lk(l,k) = 1.0/maxL_start ; 
    }
  }
  arma::vec pi_k(J*3) ;
  pi_k.fill(1.0/(J*3)) ;
  
  clusterO.col(0) = clO ;
  clusterD.col(0) = clD ;
  out_A.col(0) = A_start; 
  
  
  ///////////////////// debug
  arma::mat pi_k_out = arma::zeros(J*3, Nrep)  ;
  //////////////////////////////
  
  bool display_progress = true ;
  Progress p(Nrep, display_progress) ;
  
  for(int i = 0; i < Nrep -1 ; i++)
  {
    if( Progress::check_abort() ) { return -1.0 ; }
    p.increment();
    int check = 0 ;
    
    AA(0) = 0 ;
    for(int j = 1; j < n+1; j++) { AA(j) =  out_A(clusterO(j-1,i), i) ; }
    
    /*
     * Sampling calcium level
     * Kalman filter + backward sampling
     */
    a(0) = 0 ; // a0
    R(0) = varC0 ; // R0
    filter_mean(0) = 0 ;
    filter_var(0) = varC0 ;
    
    for(int j = 1; j < n +1 ; j++)
    {
      a(j) = out_gamma(i) * filter_mean(j-1) + AA(j);
      R(j) = out_gamma(i) * out_gamma(i) * filter_var(j-1) + out_tau2(i) ;
      
      filter_mean(j) = a(j) + R(j) / (R(j) + out_sigma2(i) ) * (y(j-1) - out_b(i) - a(j)) ;
      filter_var(j) = out_sigma2(i) * R(j) / (R(j) + out_sigma2(i)) ;
    }
    out_c(n, i+1) = R::rnorm(filter_mean(n), filter_var(n)) ;
    
    for(int j = n-1; j > -1; j--)
    {
      back_mean = filter_mean(j) + out_gamma(i) * filter_var(j) / R(j+1) * (out_c(j+1, i+1) - a(j+1)) ;
      back_var = filter_var(j) - pow(out_gamma(i) * filter_var(j), 2) / R(j+1) ;
      
      out_c(j, i+1) = R::rnorm(back_mean, back_var) ;
    }
    // out_c.col(i+1) = out_c.col(i) ;
    
    
    /*
     * Sampling b
     */
    arma::vec z(n) ; 
    for(int j = 0; j < n; j++) { z(j) = y(j) - out_c(j+1, i+1) ; } 
    
    out_b(i+1) = R::rnorm( (hyp_b2 * hyp_b1 + 1/out_sigma2(i) * arma::accu(z)) / (hyp_b2 + n / out_sigma2(i)) , 
          std::sqrt( 1/ (hyp_b2 + n / out_sigma2(i)) ) ) ;
    // out_b(i+1) = out_b(i) ;
    
    
    /*
     * Sampling sigma2
     */
    arma::vec sq(n) ; 
    for(int j = 0; j < n; j++) { sq(j) = (z(j) - out_b(i+1)) * (z(j) - out_b(i+1)) ; }
    
    out_sigma2(i+1) = 1/ R::rgamma(hyp_sigma21 + n/2, 1/(hyp_sigma22 + 0.5 * arma::accu(sq)) ) ;
    // out_sigma2(i+1) = out_sigma2(i) ;
    
    
    /*
     * Sampling tau2
     */
    arma::vec sq2(n) ; 
    for(int j = 0; j < n ; j++) { sq2(j) = (out_c(j+1,i+1) - out_gamma(i) * out_c(j,i+1) - AA(j+1)) * 
       (out_c(j+1,i+1) - out_gamma(i) * out_c(j,i+1) - AA(j+1)) ; }
    out_tau2(i+1) = 1/ R::rgamma(hyp_tau21 + n/2, 1/(hyp_tau22 + 0.5 * arma::accu(sq2)) ) ;
    // out_tau2(i+1) = out_tau2(i) ;
    
    
    /*
     * Sampling gamma
     * MH step with uniform random walk
     */
    oldgamma = out_gamma(i) ;
    newgamma = oldgamma + R::runif(-eps_gamma, eps_gamma) ;
    ratio = exp( logpost_gamma(y, out_c.col(i+1), 
                         AA, out_b(i+1), 
                         newgamma, out_sigma2(i+1), out_tau2(i+1),
                         hyp_gamma1, hyp_gamma2) -
                  logpost_gamma(y, out_c.col(i+1), 
                          AA, out_b(i+1), 
                          oldgamma, out_sigma2(i+1), out_tau2(i+1),
                          hyp_gamma1, hyp_gamma2) ) ;
    if(R::runif(0, 1) < ratio) oldgamma = newgamma ;
    out_gamma(i+1) = oldgamma ;
    
    
    /*
     * Sampling of clusters and cluster parameters
     * Slice sampler
     */
    out_slice = nested_gMFM(y, g, 
                              clusterD.col(i), clusterO.col(i), 
                              out_c.col(i+1),
                              out_A.col(i), out_b(i+1), out_gamma(i+1), 
                              out_p(i),
                              out_sigma2(i+1), out_tau2(i+1),
                              out_alpha(i), out_beta(i), 
                              pi_k, omega_lk,
                              out_maxK(i),
                              out_maxL(i),
                              hyp_alpha1, hyp_alpha2, eps_alpha, 
                              hyp_beta1, hyp_beta2, eps_beta, 
                              hyp_A1, hyp_A2,  
                              hyp_maxK1, hyp_maxK2, hyp_maxK3,
                              hyp_maxL1, hyp_maxL2, hyp_maxL3,
                              eps_A,
                              eps_maxK,
                              eps_maxL,
                              check) ;

    arma::vec out_slice_clusO = out_slice["clusterO"] ;
    clusterO.col(i+1) = out_slice_clusO ;
    //clusterO.col(i+1) = clusterO.col(i) ;
    
    arma::vec out_slice_clusD = out_slice["clusterD"] ;
    clusterD.col(i+1) = out_slice_clusD ;
    //clusterD.col(i+1) = clusterD.col(i) ;
    
    arma::vec out_slice_A = out_slice["A"] ;
    out_A.col(i+1) = out_slice_A ;
   // out_A.col(i+1) = out_A.col(i) ;
    
    double out_slice_alpha = out_slice["alpha"] ;
    out_alpha(i+1) = out_slice_alpha ;
 
    double out_slice_beta = out_slice["beta"] ;
    out_beta(i+1) = out_slice_beta ;
    
    arma::mat omega_lk_tmp = out_slice["omega_lk"] ;
    omega_lk = omega_lk_tmp ;

    int out_slice_maxK = out_slice["maxK"] ;
    out_maxK(i+1) = out_slice_maxK ;
    
    int out_slice_maxL = out_slice["maxL"] ;
    out_maxL(i+1) = out_slice_maxL ;
    
    arma::vec pi_k_tmp = out_slice["pi_k"] ;
    pi_k = pi_k_tmp ;
    for(int k = 0; k < out_maxK(i+1); k++) pi_k_out(k,i+1) = pi_k(k) ;
    
    /*
     * Sampling p (spike and slab proportion)
     */
    arma::vec line(n) ; line = clusterO.col(i+1) ;
    double n0 = std::count(line.begin(), line.end(), 0) ;
    
    out_p(i+1) = R::rbeta(hyp_p1 + n - n0, hyp_p2 + n0) ;
    if(out_p(i+1) > 0.4) { check = 1 ; }
    //   out_p(i+1) = out_p(i) ;

    
    
    //// END Gibbs sampler ////
    if(check == 1) { 
      Rcout << "Stop at iter. " << i << "\n" ;
      i = Nrep - 2 ; 
    }
  }
  return Rcpp::List::create(Rcpp::Named("calcium") = out_c,
                            Rcpp::Named("A") = out_A,
                            Rcpp::Named("b") = out_b,
                            Rcpp::Named("gamma") = out_gamma,
                            Rcpp::Named("sigma2") = out_sigma2,
                            Rcpp::Named("clusterO") = clusterO,
                            Rcpp::Named("clusterD") = clusterD,
                            Rcpp::Named("p") = out_p,
                            Rcpp::Named("tau2") = out_tau2,
                            Rcpp::Named("alpha") = out_alpha,
                            Rcpp::Named("beta") = out_beta,
                            Rcpp::Named("maxK") = out_maxK,
                            Rcpp::Named("maxL") = out_maxL,
                            Rcpp::Named("pi_k_out") = pi_k_out) ;
}



