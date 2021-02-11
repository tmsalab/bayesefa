#include <RcppArmadillo.h>
#include "commonEFAfunctions.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


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

//' @title Rotate loadings and factor scores to a permuted positive lower triangle
//' @description Rotates loading matrix according to specified founding variable row indices. 
//' @param new_r_idx vector of row indices of length equal to the number of factors.
//' @param lambda Loading matrix.
//' @param factors A n x m matrix of factor scores that is rotated. 
//' @author Albert X Man
//' @return A list of rotated loadings and factors.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List proposal2(arma::uvec& new_r_idx, // Input in R format, aka counting starts from 1 not 0
                     arma::mat& lambda, 
                     arma::mat& factors) {
  arma::mat Q, R;
  //extract lambda_r
  arma::mat lambda_sub = lambda.rows(new_r_idx - 1);
  
  //compute QR decomposition
  qr(Q, R, lambda_sub.t());
  arma::vec sgn = sign(R.diag());
  //compute Q for rotating to new row index configuration
  Q = Q * diagmat(sgn);
  arma::mat lambda_new = lambda * Q;
  arma::mat f_new = factors * Q;
  
  return Rcpp::List::create(Rcpp::Named("lambda_new",lambda_new),
                            Rcpp::Named("f_new",f_new));
}


//' @title Propose Mode-Jumping Rotation
//' @description Function to propose a new rotation 
//' @param X A N by J matrix of responses.
//' @param lambda_mean A J by M matrix of loadings.
//' @param f_mean A N by M matrix of factor scores.
//' @param invClam The loading precision hyper-parameter.
//' @param sigma This does not appear to be used.
//' @param r_idx A M vector of permuted positive lower triangular (PPLT) row indices.
//' @param my_gamma The mode-jumping tuning parameter.
//' @author Albert X Man
//' @return A list containing:
//' 
//' - `lambda`: Rotated loading matrix. 
//' - `r_idx`: PPLT row indices.
//' - `f_mat`: Rotated factor scores.
//' - `accepted`: An indicator denoting whether a rotated candidate was accepted.
//' 
//' @noRd
// [[Rcpp::export]]
Rcpp::List mode_jump(const arma::mat& X, 
                     arma::mat& lambda_mean, 
                     arma::mat& f_mean, 
                     arma::mat& invClam, 
                     const arma::vec sigma, 
                     arma::uvec r_idx, // Input in Rcpp format; index starts at 0 not 1
                     double my_gamma) {
  int J = lambda_mean.n_rows;
  int M = lambda_mean.n_cols;
  bool accepted = false;
  
  //generate and sample proposed row indices for PLT configuration
  arma::uvec new_r_idx;
  arma::uvec one_to_J = arma::shuffle(arma::linspace<arma::uvec>(1, J, J));
  new_r_idx = sort(one_to_J.subvec(0, M-1));
  
  //rotate current lambda and F to proposed configuration
  Rcpp::List forward_prop = proposal2(new_r_idx, lambda_mean, f_mean);
  
  //computations for current state
  arma::mat lambda_sub = lambda_mean.rows(r_idx);
  double log_like_old = 0;
  // Lambda prior variance
  //compute the log determinant for current state
  log_like_old += accu(log(pow(lambda_sub.diag(), my_gamma)));
  
  //computations for proposed rotation
  //extract rotated lambda and F
  arma::mat lambda_proposed = forward_prop["lambda_new"];
  arma::mat F_proposed = forward_prop["f_new"];
  
  //compute logged determinant
  arma::mat lambda_sub_proposed = lambda_proposed.rows(new_r_idx-1);
  double log_like_prop = 0;
  log_like_prop += accu(log(pow(lambda_sub_proposed.diag(), my_gamma)));
  
  // double marginal_like = NULL;
  // double marginal_like_prop = NULL;
  
  //compute acceptance probability
  // double A = std::min(1.0, exp(log_like_prop - log_like_old));
  // if (!std::isfinite(A)) {
  //   std::cout << "NA acceptance probability\n";
  //   return Rcpp::List::create(Rcpp::Named("lambda",lambda_mean),
  //                             Rcpp::Named("r_idx",r_idx),
  //                             Rcpp::Named("f_mat",f_mean),
  //                             Rcpp::Named("accepted",accepted),
  //                             Rcpp::Named("marginal_like",log_like_old));
  // }
  
  double u = arma::randu<double>();
  
  if (log(u) < log_like_prop - log_like_old) {
    // if (u < A) {
    arma::mat temp = forward_prop["f_new"];
    f_mean = temp;
    lambda_mean = lambda_proposed;
    r_idx = new_r_idx - 1;
    accepted = true;
  }
  return Rcpp::List::create(Rcpp::Named("lambda",lambda_mean),
                            Rcpp::Named("r_idx",r_idx),
                            Rcpp::Named("f_mat",f_mean),
                            Rcpp::Named("accepted",accepted));
}

//' @title Generate Random Multivariate Normal Distribution
//' @description Creates a random Multivariate Normal when given number of obs, mean, and sigma. 
//' @param n An \code{int}, which gives the number of observations.  (> 0)
//' @param mu A \code{vector} length m that represents the means of the normals.
//' @param S A \code{matrix} with dimensions m x m that provides Sigma, the covariance matrix. 
//' @return A \code{matrix} that is a Multivariate Normal distribution
//' @author James J Balamuta
//' @noRd
//' 
//' @examples 
//' #Call with the following data:
//' rmvnorm(2, c(0,0), diag(2))
//' 
// [[Rcpp::export]]
arma::mat rmvnorm(unsigned int n, const arma::vec& mu, const arma::mat& S) {
  unsigned int ncols = S.n_cols;
  arma::mat Y(n, ncols);
  Y.imbue( norm_rand ) ;
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(S);
}

//' @title Simulate a Gamma-Type Random Variable
//' @description Simulates a gamma-type random variable with the slice sampler. 
//' @param x The value of the given parameter.
//' @param alph The gamma-type shape parameter.
//' @param g The value of the g parameter.
//' @param d The value of the d parameter
//' @author Albert X Man
//' @return A scalar that is a gamma-type distribution.
//' @noRd
//' 
// [[Rcpp::export]]
double sim_gamma_type(double x,double alph,double g,double d){
  //update u given x
  double v1 = R::runif(0.0,1.0);
  double u = v1*exp(-pow(x-g,2)/d);
  //update x given u
  double v2 = R::runif(0.0,1.0);
  arma::vec acand(2);
  acand(0) = 0.;
  acand(1) = g - sqrt(-1.*d*log(u));
  double a = arma::max(acand);
  double b = g + sqrt(-1.*d*log(u));
  double x1 = pow(v2*pow(b,alph)+(1-v2)*pow(a,alph),1./alph);
  return x1;
}

//update F
//this could be updated to sample F within an i loop, which may be advantageous for ordinal samples
//' @title Sample Factor Scores
//' @description Sample factor scores from the posterior distribution. 
//' @param Y A N by J matrix of responses or augmented data.
//' @param I_K A M by M identity matrix.
//' @param F A N by M matrix of factor scores.
//' @param Lambda A J by M matrix of loadings.
//' @param psis_inv A J vector of inverse uniquenessess.
//' @author Albert X Man, Steven A Culpepper
//' 
//' @noRd
// [[Rcpp::export]]
void update_F_matrix(const arma::mat& Y,
                     const arma::mat& I_K,
                     arma::mat& F,
                     arma::mat& Lambda,
                     arma::vec& psis_inv){
  unsigned int N = F.n_rows;
  unsigned int M = F.n_cols;
  arma::mat Psi_inv = arma::diagmat(psis_inv);
  arma::mat Psi_inv_Lam = Psi_inv*Lambda;
  arma::mat VF = inv(I_K + Lambda.t()*Psi_inv_Lam);
  arma::mat mF = Y*Psi_inv_Lam*VF;
  F = mF + arma::randn<arma::mat>(N,M)*arma::chol(VF); 
  
}


//update precision for elements of lambda
//note there is just single element for the prior
//' @title Sample the Loading Precision Hyperparameter.
//' @description Sample the loading precision parameter from the posterior distribution.
//' @param Lambda A J by M matrix of loadings.
//' @param invClam The loading precision hyper-parameter.
//' @author Albert X Max
//' 
//' @noRd
// [[Rcpp::export]]
void update_invClam(const arma::mat& Lambda,
                    arma::mat& invClam) {
  // No spike-and-slab, only column-wise variance
  unsigned int J = Lambda.n_rows;
  int M = Lambda.n_cols;
  double u = R::runif(0.0,1.0);
  double d = accu(pow(Lambda, 2));
  invClam.fill(R::qgamma(u,double(J*M-((M-1)*M/2) + 2)/2.0,2.0/(2.0 + d),1,0));
}


//' @title Set difference function
//' @description Find indices in x that are not included in y. 
//' @param x The first set.
//' @param y The second set
//' @return A vec of elements that are in x and not y.
//' @author Albert X Man
//' @noRd
//' 
// [[Rcpp::export]]
arma::uvec my_setdiff(arma::uvec& x, const arma::uvec& y){
  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uword q1 = arma::conv_to<arma::uword>::from(arma::find(x == y(j)));
    x.shed_row(q1);
  }
  return x;
}

//update loadings given Y, row indices, F, and uniquenesses
//' @title Sample the Loadings.
//' @description Sample the loadings from the posterior distribution.
//' @param Y A N by J matrix of responses or augmented data.
//' @param r_idx A M vector of permuted positive lower triangular (PPLT) row indices.
//' @param F A N by M matrix of factor scores.
//' @param Lambda A J by M matrix of loadings.
//' @param psis_inv A J vector of inverse uniquenessess.
//' @param invClam The loading precision hyper-parameter.
//' @param my_gamma The mode-jumping tuning parameter.
//' @author Albert X Man
//' 
//' @noRd
// [[Rcpp::export]]
void update_Lambda_loadings_hard_zero(const arma::mat& Y,
                                      arma::uvec& r_idx,
                                      arma::mat& F,
                                      arma::mat& Lambda,
                                      arma::vec& psis_inv,
                                      const arma::mat invClam,
                                      double my_gamma){
  unsigned int J = Lambda.n_rows;
  unsigned int M = F.n_cols;
  //update Lambda 
  // Lambda.fill(0);
  arma::mat FpF = F.t()*F;
  arma::mat FpY = F.t()*Y;
  arma::rowvec lambdaj(M);
  
  // Constrained case
  unsigned int j;
  for(unsigned int m=0;m<M;++m){
    j = r_idx(m);
    arma::mat FpF_1toj = FpF(arma::span(0,m),arma::span(0,m));
    arma::mat invClamj = arma::diagmat(invClam.submat(j,0,j,m));
    arma::vec FpYj = FpY(arma::span(0,m),j);
    arma::mat C_1toj = inv(psis_inv(j)*FpF_1toj+invClamj);
    arma::vec m_1toj = psis_inv(j)*C_1toj*FpYj;
    double lambdajj=sim_gamma_type(Lambda(j,m),my_gamma+1.0,m_1toj(m),2.*C_1toj(m,m));
    
    if (!arma::is_finite(lambdajj)) {
      lambdajj = 0;
    }
    // update lambdaj conditioned upon lambdajj
    if (m > 0) {
      arma::vec m1 = m_1toj(arma::span(0,m-1));
      arma::mat S1 = C_1toj(arma::span(0,m-1),arma::span(0,m-1));
      arma::vec S12= C_1toj(arma::span(0,m-1),m);
      arma::vec mb = m1+S12*(lambdajj-m_1toj(m))/C_1toj(m,m);
      arma::mat Sb = S1-S12*S12.t()/C_1toj(m,m);
      arma::rowvec lambda_1tjm1 = rmvnorm(1,mb,Sb);
      Lambda(j,arma::span(0,m-1))=lambda_1tjm1;
    }
    Lambda(j,m) = lambdajj;
    
    //fill remaining elements in PLT rows with zeros
    if ((m+1) < M) {
      arma::rowvec fill_zeros = arma::zeros<arma::rowvec>(M-(m+1));
      Lambda(j,arma::span(m+1,M-1)) = fill_zeros;
    }
  }
  
  // Unconstrained case
  arma::uvec one_to_J = arma::linspace<arma::uvec>(0, J-1, J);
  arma::uvec not_r_idx = my_setdiff(one_to_J, r_idx);//find indices not in r
  unsigned int I = not_r_idx.n_elem;
  for(unsigned int i=0; i<I; ++i) {
    j = not_r_idx(i);
    arma::mat Cj = inv(psis_inv(j)*FpF+arma::diagmat(invClam.row(j)));
    arma::vec mj = psis_inv(j)*Cj*FpY.col(j);
    lambdaj = rmvnorm(1,mj,Cj);
    Lambda.row(j) = lambdaj;
  }
}

//' @title Initialize Thresholds.
//' @description Initialize category threshold parameters.
//' @param Ms A J vector indicating the number of categories.
//' @author Steven A Culpepper
//' @noRd
// [[Rcpp::export]]
arma::mat kappa_initialize(const arma::vec& Ms){
  unsigned int J = Ms.n_elem;
  unsigned int M = max(Ms);
  arma::mat KAP0(M+1,J);
  (KAP0.row(0)).fill(-arma::datum::inf);
  KAP0.row(1).fill(0.0);
  for(unsigned int j=0;j<J;j++){
    for(unsigned int m=2;m<M+1;m++){
      if(m<Ms(j)){
        KAP0(m,j)=(m-1)*1;
      }
      if(m==Ms(j)){
        KAP0(m,j)=arma::datum::inf;
      }
      if(m>Ms(j)){
        KAP0(m,j)=arma::datum::nan;
      }
    }  
  }
  return KAP0;
}

