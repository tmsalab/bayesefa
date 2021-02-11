#include <RcppArmadillo.h>

#ifndef COMMONEFAFUNCTIONS_H
#define COMMONEFAFUNCTIONS_H

Rcpp::List proposal2(arma::uvec& new_r_idx,arma::mat& lambda_mean, arma::mat& f_mean);
  
Rcpp::List mode_jump(const arma::mat& X, arma::mat& lambda_mean, arma::mat& f_mean, 
                     arma::mat& invClam, const arma::vec sigma, arma::uvec r_idx, double my_gamma);

arma::mat rmvnorm(unsigned int n, const arma::vec& mu, const arma::mat& S);

double sim_gamma_type(double x,double alph,double g,double d);

void update_F_matrix(const arma::mat& Y,const arma::mat& I_K,arma::mat& F,arma::mat& Lambda,
                     arma::vec& psis_inv);

void update_invClam(const arma::mat& Lambda,arma::mat& invClam);

arma::uvec my_setdiff(arma::uvec& x, const arma::uvec& y);

void update_Lambda_loadings_hard_zero(const arma::mat& Y,arma::uvec& r_idx,arma::mat& F,
                                      arma::mat& Lambda,arma::vec& psis_inv,
                                      const arma::mat invClam,double my_gamma);

arma::mat kappa_initialize(const arma::vec& Ms);

#endif
