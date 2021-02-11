#include <RcppArmadillo.h>
#include <commonEFAfunctions.h>



//unique to EFA code
//sample unique factor variances
//' @title Sample Uniquenesses.
//' @description Sample the unique factor variances from the posterior distribution.
//' @param Y A N by J matrix of responses or augmented data.
//' @param F A N by M matrix of factor scores.
//' @param Lambda A J by M matrix of loadings.
//' @param psis_inv A J vector of inverse uniquenessess.
//' @author Albert X Man
//' 
//' @noRd
// [[Rcpp::export]]
void update_uniquenesses(const arma::mat& Y,
                         arma::mat& F,
                         arma::mat& Lambda,
                         arma::vec& psis_inv) {
  // One uniqueness variance per parameter
  // Shared between residuals and loadings
  
  // unsigned int M = Lambda.n_cols;
  int J = Lambda.n_rows;
  unsigned int N = F.n_rows;
  
  double b1 = 1;
  double b2 = 1;
  arma::mat Y_hat = F * Lambda.t();
  arma::mat Y_resids = Y - Y_hat;
  arma::colvec d2 = trans( sum( pow(Y_resids, 2), 0) );

  double u;
  int n_obs;
  for (int j = 0; j < J; j++) {
    u = R::runif(0.0,1.0);
    n_obs = N;
    // qgamma is parameterized as shape, scale rather than shape and rate
    psis_inv(j) = (R::qgamma(u,double(n_obs + 2.0*b1)/2.0,2.0/(2.0*b2 + d2[j]),1,0));
  }
}





//' @title Exploratory Factor Analysis of Continuous Response Data
//' @description Implement the Man and Culpepper (2020) mode-jumping algorithm to factor analyze continuous response data. 
//' @param Y A N by J matrix of mean-centered, continuous variables. 
//' @param M An integer specifying the number of factors. 
//' @param gamma The value of the mode-jumping tuning parameter. Man and Culpepper (2020) used gamma = 0.5. 
//' @param burnin Number of burn-in iterations to discard.
//' @param chain_length The total number of iterations (burn-in + post-burn-in).
//' 
//' @return A list that contains nsamples = chain_length - burnin array draws from the posterior distribution:
//' 
//' - `LAMBDA`: A J by M by nsamples array of sampled loading matrices.
//' - `PSIs`: A J by nsamples matrix of vector of variable uniquenesses.
//' - `ROW_OUT`: A matrix of sampled row indices of founding variables for mode-jumping algorithm.
//' - `F_OUT`: An array of sampled factor scores. 
//' - `ACCEPTED`: Acceptance rates for mode-jumping Metropolis-Hastings (MH) steps.
//' 
//' @author Steven Andrew Culpepper, Albert Man
//' 
//' @references
//' 
//' Man, A. & Culpepper, S. A. (2020). A mode-jumping algorithm for Bayesian factor analysis. Journal of the American Statistical Association, doi:10.1080/01621459.2020.1773833.
//' 
//' @export
//' 
//' @examples
//' 
//' data(exchangerate)
//' 
//' #Retain complete cases and drop Month_Year column
//' X<-exchangerate[complete.cases(exchangerate),-1]
//' X<-apply(X,2, diff)
//' X<-as.matrix(scale(X))
//' 
//' #Specify the number of factors
//' my_M<-2
//' 
//' #Run the mode-jumping EFA algorithm
//' burn_in<-150
//' chain_length<-300
//' out <- EFA_Mode_Jumper(X,my_M,gamma=0.5,burnin=burn_in,chain_length)
//'   
//' #Rotate all samples to permutation positive lower triangular (PPLT) 
//' #structure with USD and FRF as factor founding variables
//'   my_lambda_sample = out$LAMBDA
//'     for (j in 1:dim(my_lambda_sample)[3]) {
//'       my_rotate = proposal2(c(1,4),my_lambda_sample[,,j],out$F_OUT[,,j])
//'       my_lambda_sample[,,j] = my_rotate$lambda
//'     }
//'     
//' #compute posterior mean of PPLT loadings
//' mLambda<-apply(my_lambda_sample,c(1,2),mean)
//'       
// [[Rcpp::export]]
Rcpp::List EFA_Mode_Jumper(const arma::mat& Y,
                           unsigned int M,
                           double gamma,
                           unsigned int burnin,
                           unsigned int chain_length=10000){
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  double my_gamma = gamma;
  unsigned int chain_m_burn = chain_length-burnin;
  unsigned int tmburn;
  arma::mat I_K = arma::eye(M,M);
  
  //Saving output
  arma::cube LAMBDA(J,M,chain_m_burn);
  arma::cube F_OUT(N,M,chain_m_burn);
  arma::mat PSIs(J,chain_m_burn);
  arma::umat ROW_OUT(M, chain_m_burn);
  // arma::mat MH_PROBS(4, chain_m_burn);
  arma::vec ACCEPTED(chain_m_burn);
  arma::vec LIKELIHOOD(chain_m_burn);
  
  // Initialize r_idx as PLT 
  arma::uvec r_idx;
  // arma::uvec one_to_J = arma::shuffle(arma::linspace<arma::uvec>(0, J-1, J));
  arma::uvec one_to_J = arma::linspace<arma::uvec>(0, J-1, J);
  r_idx = sort(one_to_J.subvec(0, M-1));
  
  arma::mat F = arma::randn<arma::mat>(N,M);
  // arma::mat Cstart=inv(I_K+F.t()*F);
  
  arma::mat Lambda = arma::randn<arma::mat>(J,M);
  for(unsigned int j=0;j<M;++j){
    Lambda(j,j) = sim_gamma_type(1.,my_gamma+1.0,0.,2.);
    for(unsigned int k=j+1;k<M;++k){
      Lambda(j,k) = 0.;
    }
  }
  
  arma::vec psis_inv(J);
  // arma::rowvec varY=arma::var(Y);
  for(unsigned int j=0;j<J;++j){
    arma::rowvec lambdaj = Lambda.row(j);
    arma::vec Y_m_Flam = Y.col(j)-F*lambdaj.t();
    double u = R::runif(0.0,1.0);
    //qgamma is parameterized as shape, scale rather than shape and rate
    psis_inv(j) = R::qgamma(u,double(N+2.)/2.0,
             2.0/(2.+arma::dot(Y_m_Flam,Y_m_Flam)),1,0);
  }
  
  // arma::vec omega(M, arma::fill::ones);
  arma::mat invClam(J, M, arma::fill::ones);

  // Main Algorithm Loop
  double accepted;
  // arma::vec metropolis_probs;
  Rcpp::List L;
  // double marginal_like;
  
  // arma::rowvec Y_mean = arma::zeros<arma::rowvec>(J);
  // arma::vec marginal_like_vec;
  
  for(unsigned int t = 0; t < chain_length; ++t){
    
    accepted = 0; 
    update_Lambda_loadings_hard_zero(Y, r_idx, F, Lambda, psis_inv, invClam, my_gamma);
    update_F_matrix(Y,I_K, F, Lambda, psis_inv);
    update_invClam(Lambda, invClam);
    update_uniquenesses(Y, F, Lambda, psis_inv);
    //perform mode jump step
    L = mode_jump(Y, Lambda, F, invClam, pow(psis_inv, -0.5), r_idx, my_gamma);
    arma::mat temp = L["lambda"];
    Lambda = temp;
    arma::uvec temp2 = L["r_idx"];
    r_idx = temp2;
    arma::mat temp3 = L["f_mat"];
    F = temp3;
    double temp5 = L["accepted"];
    accepted = temp5;
    // double temp6 = L["marginal_like"];
    // marginal_like = temp6;
    
    // Save output after burn-in
    if(t>(burnin-1)){
      tmburn = t-burnin;
      ACCEPTED(tmburn) = accepted;
      LAMBDA.slice(tmburn) = Lambda;
      F_OUT.slice(tmburn) = F;
      PSIs.col(tmburn) = 1./psis_inv;
      ROW_OUT.col(tmburn) = r_idx;
      // LIKELIHOOD(tmburn) = marginal_like;
      // MH_PROBS.col(tmburn) = metropolis_probs;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("LAMBDA",LAMBDA),
                            Rcpp::Named("PSIs",PSIs),
                            Rcpp::Named("ROW_OUT", ROW_OUT),
                            Rcpp::Named("F_OUT", F_OUT),
                            // Rcpp::Named("MH_PROBS", MH_PROBS),
                            Rcpp::Named("ACCEPTED", ACCEPTED)//,Rcpp::Named("LIKELIHOOD", LIKELIHOOD)
                              );
}

