#include <RcppArmadillo.h>
#include <commonEFAfunctions.h>




//unique to IFA and mixed responses
//' @title Sample Bounded Normal Random Variables.
//' @description Sample truncated normal random variables.
//' @param mean location parameter.
//' @param sd scale parameter.
//' @param upper upper bound.
//' @param bound lower bound.
//' @author Steven A Culpepper
//' 
//' @noRd
// [[Rcpp::export]]
double rTruncNorm_bounds(double mean,double sd, double upper, double bound){
  double p0 = R::pnorm(bound,mean,sd, 1, 0);
  double p1 = 1-p0;
  double uZ = R::runif(0,1);
  double Z = R::qnorm( (1-upper)*uZ*p0+upper*(p0 + uZ*p1), mean, sd, 1, 0);
  return(Z);
}

//unique to mixed response
//' @title Sample Unique Factor Variances for Mixed-Type Variables.
//' @description Sample unique factor variables for more general ordinal and bounded model.
//' @param Y A N by J matrix of responses or augmented data.
//' @param Ms model indicator where 0 = "bounded", 1 = "continuous", 2 = "binary", >2 = "ordinal".
//' @param F A N by M matrix of factor scores.
//' @param Lambda A J by M matrix of loadings.
//' @param psis_inv A J vector of inverse uniquenessess.
//' @param continuous_indicator A J vector indicated whether a variable is continuous.
//' @author Albert X Man
//' 
//' @noRd
// [[Rcpp::export]]
void update_uniquenesses_mixed(const arma::mat& Y,
                         const arma::vec& Ms,
                         arma::mat& F,
                         arma::mat& Lambda,
                         arma::vec& psis_inv,
                         const arma::vec& continuous_indicator) {
  // One uniqueness variance per parameter
  // Shared between residuals and loadings
  
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
    if (continuous_indicator(j) == 1) {
      u = R::runif(0.0,1.0);
      n_obs = N;
      
      // qgamma is parameterized as shape, scale rather than shape and rate
      psis_inv(j) = (R::qgamma(u,double(n_obs + 2.0*b1)/2.0,2.0/(2.0*b2 + d2[j]),1,0));
    }
  }
}


//' @title Sample augmented data for mixed-type model
//' @description Sample augmented data for binary, ordinal, and bounded response data 
//' @param Y A N by J matrix of ordinal responses.
//' @param MISS A N by J matrix of missing data indicators.
//' @param Z A N by J matrix of augmented data.
//' @param as A J by M matrix of loadings.
//' @param bs A J vector of intercepts.
//' @param theta A N by M matrix of factor scores.
//' @param Ms A vector of the number of score categories. Note 2 = two parameter normal ogive, >2 ordinal normal ogive.
//' @param Kaps A matrix of category thresholds.
//' @param sdMH A vector of tuning parameters of length J for implementing the Cowles (1996) Metropolis-Hastings threshold sampler.
//' @param psis_inv A J vector of inverse uniquenessess.
//' @param bounds A J by 2 matrix denoting the min and max variable values. Note that bounds are only used for variable j if element j of Ms is zero.
//' @author Steven A Culpepper
//' @return A list containing:
//' 
//' - `Kaps_new`: An updated matrix of category thresholds.
//' - `MHaccept`: A binary vector indicating whether an MH proposal for thresholds was accepted. 
//'  
//' @noRd
// [[Rcpp::export]]
Rcpp::List update_WKappaZ_NA_mixed(const arma::mat& Y,//binary data
                             const arma::mat MISS,//missing indicators
                             arma::mat& Z,
                             const arma::mat& as,//lambda matrix is called "as" here
                             const arma::vec& bs,//intercept is called "bs" here
                             const arma::mat& theta,//F is called "theta" here
                             const arma::vec& Ms,
                             const arma::mat& Kaps,//thresholds is called "Kaps" here
                             const arma::vec& sdMH,
                             const arma::vec& psis_inv,
                             const arma::mat& bounds) {
  
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  // double Yij;
  // arma::vec eta(N);
  arma::vec pk_on(4);//Vector that includes probabilities from Cowles (1996) in order: 
  //P(k_j-1,new|k_j), P(k_j+1|k_j), P(k_j-1|k_j,new), P(k_j+1,new|k_j,new)
  arma::vec KapsMH(max(Ms)+1);
  arma::vec p2pno(2);//,p3pno(2); //P(0|eta_ij) for 2pno draw of truncated normals
  double mpart,ipart; //Ratio for m and i parts of MH step
  // arma::vec W(N);
  arma::vec oneN = arma::ones<arma::vec>(N);
  double uMH, uR, uZ;//, u3pno,u4pno;
  arma::vec pky_on(4); 
  arma::vec pky_new(2); 
  arma::vec MHaccept=arma::zeros<arma::vec>(J);
  arma::mat Kaps_new = Kaps;
  //  double sumz_eta,uNC;
  
  for(unsigned int j=0;j<J;j++){
    // W = arma::zeros<arma::vec>(N);
    //Calculate eta_ij for all of the models and individuals with data
    //old code that computed eta for only nonmissing values
    // for(unsigned int i=0;i<N;i++){
    //   if(MISS(i,j)==1.0){
    //     eta(i) = as(j)*theta(i,area(j)-1) - bs(j);
    //   }  
    // }
    
    arma::vec eta = bs(j)+theta*(as.row(j)).t();
    
    //3PNO and 4PNO code from NAEP_MeasurementModel_050117.cpp used to go here     
    
    ////////////////////////////////////////////////
    //ONO Model
    ////////////////////////////////////////////////
    //Declare items labeled as ONO with 2 categories as 2PNO
    if(Ms(j)>2){
      //Cowles MH sampling for items with > 3 categories
      //    if(Ms(j)>3){
      
      //Sample MH Threshold candidates //
      //Assuming that Kaps is initialized with -Inf,0,kappas,Inf in a matrix
      KapsMH = Kaps.col(j);
      
      for(unsigned int m=0; m<Ms(j)-2; m++){
        //        Rcpp::Rcout << m  << std::endl;
        uMH = R::runif(0.0,1.0);
        pk_on(0) = R::pnorm(KapsMH(m+1),Kaps(m+2,j),sdMH(j), 1, 0);
        pk_on(1) = R::pnorm(Kaps(m+3,j),Kaps(m+2,j),sdMH(j), 1, 0);
        
        KapsMH(m+2) = R::qnorm(pk_on(0) + uMH*(pk_on(1)-pk_on(0) ),Kaps(m+2,j),sdMH(j),1,0) ;
      }
      
      //Step 1a: compute m part related to thresholds
      // mpart = 1.0;
      mpart = .0;
      for(unsigned int m=0; m<Ms(j)-2; m++){
        pk_on(0) = R::pnorm(KapsMH(m+1),Kaps(m+2,j),sdMH(j), 1, 0);
        pk_on(1) = R::pnorm(Kaps(m+3,j),Kaps(m+2,j),sdMH(j), 1, 0);
        pk_on(2) = R::pnorm(Kaps(m+1,j),KapsMH(m+2),sdMH(j), 1, 0);
        pk_on(3) = R::pnorm(KapsMH(m+3),KapsMH(m+2),sdMH(j), 1, 0);
        
        //mpart = (  (pk_on(1)-pk_on(0))/(pk_on(3)-pk_on(2))  )* mpart;
        mpart = log(  (pk_on(1)-pk_on(0))/(pk_on(3)-pk_on(2))  ) + mpart;
      }
      
      //Step 1b: compute i part related to individuals
      // ipart=1.0;
      ipart=.0;
      for(unsigned int i=0;i<N;i++){
        
        if(MISS(i,j)==1.0){
          pky_on(0) = R::pnorm(KapsMH(Y(i,j)      ),eta(i),1.0, 1, 0);
          pky_on(1) = R::pnorm(KapsMH(Y(i,j)+1.0  ),eta(i),1.0, 1, 0);
          pky_on(2) = R::pnorm(Kaps(Y(i,j)      ,j),eta(i),1.0, 1, 0);
          pky_on(3) = R::pnorm(Kaps(Y(i,j)+1.0  ,j),eta(i),1.0, 1, 0);
          
          // ipart = (  (pky_on(1)-pky_on(0))/(pky_on(3)-pky_on(2))  )* ipart;
          ipart = log(  (pky_on(1)-pky_on(0))/(pky_on(3)-pky_on(2))  )+ ipart;
        }
        
      }
      
      //Step 2: Compute R and min prob for rejection & Update Kappas and Z
      //NOte that R_MH = mpart * ipart 
      uR = R::runif(0.0,1.0);
      
      if(log(uR) < mpart + ipart){
        MHaccept(j)=1.0;
        Kaps_new.col(j) = KapsMH;
        
        for(unsigned int i=0;i<N;i++){
          
          uZ = R::runif(0.0,1.0);
          
          if(MISS(i,j)==1.0){
            pky_new(0) = R::pnorm(Kaps_new(Y(i,j)    ,j)  ,eta(i),1.0, 1, 0);
            pky_new(1) = R::pnorm(Kaps_new(Y(i,j)+1.0,j)  ,eta(i),1.0, 1, 0);
            
            Z(i,j) = R::qnorm(pky_new(0) + uZ*(pky_new(1)-pky_new(0) ),eta(i),1.0,1,0) ;          
          }
          if(MISS(i,j)==0){
            Z(i,j) = R::rnorm(eta(i),1.0) ;          
          }
        }
      }
      /*     */
    }
    //    }
    
    ////////////////////////////////////////////////
    //2PNO Model
    ////////////////////////////////////////////////
    if(Ms(j)==2){

      for(unsigned int i=0;i<N;i++){
        
        if(MISS(i,j)==1.0){
          uZ = R::runif(0.0,1.0);
          p2pno(0) = R::pnorm(Kaps(Y(i,j)    ,j)  ,eta(i),1.0, 1, 0);
          p2pno(1) = R::pnorm(Kaps(Y(i,j)+1.0,j)  ,eta(i),1.0, 1, 0);
          Z(i,j) = R::qnorm(p2pno(0) + uZ*( p2pno(1)-p2pno(0) ),eta(i),1.0,1,0) ;
        }
        if(MISS(i,j)==0){
          Z(i,j) = R::rnorm(eta(i),1.0) ;          
        }
      }
    }
    //continuous response model option; only input missing observations from the model
    if(Ms(j)==1){
      double sqrt_psi = sqrt(1./psis_inv(j));
      for(unsigned int i=0;i<N;i++){
        if(MISS(i,j)==0){
          Z(i,j) = R::rnorm(eta(i),sqrt_psi) ;
        }
      }
    }
    //bounded continuous response model
    if(Ms(j)==0){
      double sqrt_psi = sqrt(1./psis_inv(j));
      for(unsigned int i=0;i<N;i++){
        if(MISS(i,j)==1.0){
          double Yij = Y(i,j);
          if(Yij==bounds(j,0)){
            //sample truncated above at bound(j,0)
            Z(i,j) = rTruncNorm_bounds(eta(i),sqrt_psi,0,bounds(j,0));
          }
          if(Yij==bounds(j,1)){
            //sample truncated below at bound(j,1)
            Z(i,j) = rTruncNorm_bounds(eta(i),sqrt_psi,1,bounds(j,1));
          }
        }
        if(MISS(i,j)==0){
          Z(i,j) = R::rnorm(eta(i),sqrt_psi) ;
        }
      }
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("Kaps_new",Kaps_new), 
    Rcpp::Named("MHaccept",MHaccept)  );
}



//shared bw ifa and mixed
//' @title Sample Intercept Parameters for Mixed-Type Model
//' @description Sample intercept parameters from the posterior distribution.
//' @param intercept A J vector of variable intercepts.
//' @param Z A N by J matrix of augmented data.
//' @param F A N by M matrix of factor scores.
//' @param Lambda A J by M matrix of loadings.
//' @param psis_inv A J vector of inverse-uniquenesses.
//' @param intercept_var A hyperparameter for the scale of item intercepts.
//' @author Albert X Man
//' 
//' @noRd
// [[Rcpp::export]]
void update_intercept_mixed(arma::vec& intercept,
                      arma::mat& Z,
                      arma::mat& F,
                      arma::mat& Lambda,
                      arma::vec& psis_inv,
                      double& intercept_var) {
  unsigned int N = Z.n_rows;
  unsigned int J = Z.n_cols;
  
  arma::mat Z_centered = Z - F * Lambda.t();
  
  for (unsigned int j=0; j<J; ++j) {
    double post_var = pow(N*psis_inv(j) + pow(intercept_var, -1), -1);
    double post_mean = psis_inv(j) * sum(Z_centered.col(j)) * post_var;
    
    intercept(j) = R::rnorm(post_mean, pow(post_var, 0.5));
  }
  
  
  // Update intercept variance hyperparameter
  double b1 = 1;
  double b2 = 1;
  
  double d = arma::as_scalar( sum( pow(intercept, 2) ) );

  double u;

  u = R::runif(0.0,1.0);

  // qgamma is parameterized as shape, scale rather than shape and rate
  intercept_var = pow(R::qgamma(u,
                 double(J + 2.0*b1)/2.0,
                 2.0/(2.0*b2 + d),1,0), -1);
  
}


//' @title Exploratory Factor Analysis of Mixed-Type Responses
//' @description Implement the Man and Culpepper (2020) mode-jumping algorithm to factor analyze mixed-type response data. Missing values should be specified as a non-numeric value such as NA.
//' @param Y A N by J matrix of mixed-type item responses. 
//' @param M An interger specifying the number of factors. 
//' @param gamma The value of the mode-jumping tuning parameter. Man and Culpepper (2020) used gamma = 0.5. 
//' @param Ms model indicator where 0 = "bounded", 1 = "continuous", 2 = "binary", >2 = "ordinal".
//' @param sdMH A J vector of tuning parameters for the Cowles (1996) Metropolis-Hastings sampler for ordinal data latent thresholds.
//' @param bounds A J by 2 matrix denoting the min and max variable values. Note that bounds are only used for variable j if element j of Ms is zero.
//' @param burnin Number of burn-in iterations to discard.
//' @param chain_length The total number of iterations (burn-in + post-burn-in).
//' @return A list that contains nsamples = chain_length - burnin array draws from the posterior distribution:
//' 
//' - `LAMBDA`: A J by M by nsamples array of sampled loading matrices on the standardized metric.
//' - `PSIs`: A J by nsamples matrix of vector of variable uniquenesses on the standardized metric.
//' - `ROW_OUT`: A matrix of sampled row indices of founding variables for mode-jumping algorithm.
//' - `THRESHOLDS`: An array of sampled thresholds. 
//' - `INTERCEPTS`: Sampled variable thresholds on the standardized metric.
//' - `ACCEPTED`: Acceptance rates for mode-jumping Metropolis-Hastings (MH) steps.
//' - `MHACCEPT`: A J vector of acceptance rates for item threshold parameters. Note that binary items have an acceptance rate of zero, because MH steps are never performed.
//' - `LAMBDA_unst`: An array of unstandardized loadings.
//' - `PSIs_inv_unst`: A matrix of unstandardized uniquenesses.
//' - `THRESHOLDS_unst`: Unstandardized thresholds.
//' - `INTERCEPTS_unst`: Unstandardized intercepts.
//' 
//' @author Albert X Man, Steven Andrew Culpepper
//' @references
//' Cowles, M. K. (1996), Accelerating Monte Carlo Markov chain convergence for cumulative link generalized linear models," Statistics and Computing, 6, 101-111.
//' 
//' Man, A. & Culpepper, S. A. (2020). A mode-jumping algorithm for Bayesian factor analysis. Journal of the American Statistical Association, doi:10.1080/01621459.2020.1773833.
//' @export
// [[Rcpp::export]]
Rcpp::List IFA_Mode_Jumper_MixedResponses(const arma::mat& Y,
                           unsigned int M,
                           double gamma,
                           const arma::vec& Ms,const arma::vec& sdMH,
                           const arma::mat& bounds,
                           unsigned int burnin,
                           unsigned int chain_length=10000){
  
  //"Ms" : 1 = continuous, 2 = 2pno, >2 = ono
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  double my_gamma = gamma;
  // arma::uvec which_j_is_continuous = find( Ms<2 );
  arma::vec continuous_indicator = arma::zeros<arma::vec>(J);
  // continuous_indicator.elem(which_j_is_continuous) = arma::ones<arma::vec(which_j_is_continuous.n_elem);
  continuous_indicator.elem( find( Ms<2 ) ).fill(1.);

  unsigned int chain_m_burn = chain_length-burnin;
  unsigned int tmburn;
  arma::mat I_K = arma::eye(M,M);
  
  arma::mat MISS = arma::ones<arma::mat>(N,J);
  MISS.elem( arma::find_nonfinite(Y) ).fill(0.0);
  
  //Saving output
  arma::cube LAMBDA(J,M,chain_m_burn);
  // arma::cube F_OUT(N,M,chain_m_burn);
  arma::mat PSIs(J,chain_m_burn);
  arma::umat ROW_OUT(M, chain_m_burn);
  // arma::mat MH_PROBS(4, chain_m_burn);
  arma::vec ACCEPTED(chain_m_burn);
  // arma::vec LIKELIHOOD(chain_m_burn);
  arma::cube THRESHOLDS_OUT(max(Ms)-1,J, chain_m_burn);
  arma::mat INTERCEPT_OUT(J, chain_m_burn);
  
  //saving unstandardized output
  arma::cube LAMBDA_unst(J,M,chain_m_burn);
  arma::mat PSIs_inv_unst(J,chain_m_burn);
  arma::cube THRESHOLDS_unst(max(Ms)+1,J, chain_m_burn);
  arma::mat INTERCEPT_unst(J, chain_m_burn);

  // Initialize r_idx as PLT
  arma::uvec r_idx;
  arma::uvec one_to_J = arma::linspace<arma::uvec>(0, J-1, J);
  r_idx = sort(one_to_J.subvec(0, M-1));

  arma::mat thresholds = kappa_initialize(Ms);
  arma::mat F = arma::zeros<arma::mat>(N,M);//arma::randn<arma::mat>(N,M);
  // arma::mat Cstart=inv(I_K+F.t()*F);

  // arma::mat Z = arma::randn<arma::mat>(N,J);

  arma::mat Lambda = arma::randn<arma::mat>(J,M);
  for(unsigned int j=0;j<M;++j){
    Lambda(j,j) = sim_gamma_type(1.,my_gamma+1.0,0.,2.);
    for(unsigned int k=j+1;k<M;++k){
      Lambda(j,k) = 0.;
    }
  }

  //initialize Z from fc based on starting values. let mean = 0
  // arma::mat Z = arma::randn<arma::mat>(N,J);

  arma::mat Z(N,J);
  arma::rowvec meanY=arma::mean(Y);
  arma::rowvec sdY=arma::stddev(Y);
  for(unsigned int j=0;j<J;j++){

    if(Ms(j)>2){
      arma::vec pky_new(2);

      for(unsigned int i=0;i<N;i++){
        double uZ = R::runif(0.0,1.0);
        if(MISS(i,j)==1.0){
          pky_new(0) = R::pnorm(thresholds(Y(i,j)    ,j)  ,0.0,1.0, 1, 0);
          pky_new(1) = R::pnorm(thresholds(Y(i,j)+1.0,j)  ,0.0,1.0, 1, 0);
          Z(i,j) = R::qnorm(pky_new(0) + uZ*(pky_new(1)-pky_new(0) ),0.0,1.0,1,0) ;
        }
        if(MISS(i,j)==0){
          Z(i,j) = R::rnorm(0.0,1.0) ;
        }
      }
    }

    if(Ms(j)==2){
      for(unsigned int i=0;i<N;i++){
        arma::vec p2pno(2);
        if(MISS(i,j)==1.0){
          double uZ = R::runif(0.0,1.0);
          p2pno(0) = R::pnorm(thresholds(Y(i,j)    ,j)  ,0.0,1.0, 1, 0);
          p2pno(1) = R::pnorm(thresholds(Y(i,j)+1.0,j)  ,0.0,1.0, 1, 0);
          Z(i,j) = R::qnorm(p2pno(0) + uZ*( p2pno(1)-p2pno(0) ),0.0,1.0,1,0) ;
        }
        if(MISS(i,j)==0){
          Z(i,j) = R::rnorm(0.0,1.0) ;
        }
      }
    }

    if(Ms(j)==1){//this is the continuous case
      for(unsigned int i=0;i<N;i++){
        if(MISS(i,j)==1.0){
          Z(i,j) = Y(i,j);
        }
        if(MISS(i,j)==0){
          Z(i,j) = meanY(j);//R::rnorm(0.0,1.0) ;
        }
      }
    }
    if(Ms(j)==0){//this is the continuous case
      for(unsigned int i=0;i<N;i++){
        if(MISS(i,j)==1.0){
          double Yij = Y(i,j);
          Z(i,j) = Yij;
          if(Yij==bounds(j,0)){
            //sample truncated above at bound(j,0)
            Z(i,j) = rTruncNorm_bounds(meanY(j),sdY(j),0,bounds(j,0));
            }
          if(Yij==bounds(j,1)){
            //sample truncated below at bound(j,1)
            Z(i,j) = rTruncNorm_bounds(meanY(j),sdY(j),1,bounds(j,1));
          }
        }
        if(MISS(i,j)==0){
          Z(i,j) = meanY(j);//R::rnorm(0.0,1.0) ;
        }
      }
    }
  }

  arma::vec psis_inv(J, arma::fill::ones);
  arma::rowvec varY=arma::var(Y);
  arma::vec intercept(J, arma::fill::zeros);//=(arma::mean(Z)).t();
  double intercept_var = 1;

  arma::mat invClam(J, M, arma::fill::ones);

  // Main Algorithm Loop
  double accepted;
  Rcpp::List L;


  arma::mat Z_centered = arma::zeros<arma::mat>(N, J);
  arma::vec marginal_like_vec;

  // std::cout << "begin\n";

  //save MH acceptance for thresholds
  arma::vec P_ACCEPT= arma::zeros<arma::vec>(J);
  arma::vec ACCEPT= arma::zeros<arma::vec>(J);

  arma::uvec thresholdindices = arma::linspace<arma::uvec>(1, max(Ms)-1, max(Ms)-1);

  for(unsigned int t = 0; t < chain_length; ++t){

    accepted = 0;

    Rcpp::List step1Z = update_WKappaZ_NA_mixed(Y,//binary data
                                          MISS,//missing indicators
                                          Z,
                                          Lambda,//lambda matrix is called "as" here
                                          intercept,//intercept is called "bs" here
                                          F,//F is called "theta" here
                                          Ms,
                                          thresholds,//thresholds is called "Kaps" here
                                          sdMH,
                                          psis_inv,
                                          bounds);

    thresholds  = Rcpp::as<arma::mat>(step1Z[0]);
    ACCEPT= Rcpp::as<arma::vec>(step1Z[1]);

    update_intercept_mixed(intercept, Z, F, Lambda, psis_inv, intercept_var);
    // std::cout << intercept_var << "\n";
    Z_centered = Z - repmat(intercept, 1, N).t();
    update_F_matrix(Z_centered, I_K, F, Lambda, psis_inv);
    update_Lambda_loadings_hard_zero(Z_centered, r_idx, F, Lambda, psis_inv, invClam, my_gamma);
    update_invClam(Lambda, invClam);
    update_uniquenesses_mixed(Z_centered, Ms, F, Lambda, psis_inv,continuous_indicator);//note changed Z -> Z_centered 10/15/20

    // Mode Jumping Step
    L = mode_jump(Z_centered, Lambda, F, invClam, pow(psis_inv, -0.5), r_idx, my_gamma);
    arma::mat temp = L["lambda"];
    Lambda = temp;
    arma::uvec temp2 = L["r_idx"];
    r_idx = temp2;
    arma::mat temp3 = L["f_mat"];
    F = temp3;
    double temp5 = L["accepted"];
    accepted = temp5;

    // Save output after burn-in
    if(t>(burnin-1)){

      // Correct scale invariance, put loadings, intercept, and residual variance
      // into standardized representation
      arma::mat A1 = diagmat(pow(diagvec(Lambda * Lambda.t()) + pow(psis_inv, -1), -0.5));
      arma::mat new_Lambda = A1 * Lambda;
      arma::vec new_intercept = A1 * intercept;
      arma::mat new_threshold_transposed = A1 * (thresholds.rows(thresholdindices)).t();

      arma::mat Lambda_Lambda_t = new_Lambda * new_Lambda.t();
      arma::vec psis_out = 1 - Lambda_Lambda_t.diag();

      tmburn = t-burnin;
      THRESHOLDS_OUT.slice(tmburn) = new_threshold_transposed.t();//thresholds.rows(thresholdindices);
      INTERCEPT_OUT.col(tmburn) = new_intercept;
      ACCEPTED(tmburn) = accepted;
      LAMBDA.slice(tmburn) = A1 * Lambda;
      // F_OUT.slice(tmburn) = F;
      // PSIs.col(tmburn) = 1./(A*pow(psis_inv, -1));
      PSIs.col(tmburn) = psis_out;
      ROW_OUT.col(tmburn) = r_idx;

      THRESHOLDS_unst.slice(tmburn) = thresholds;
      INTERCEPT_unst.col(tmburn) = intercept;
      LAMBDA_unst.slice(tmburn) = Lambda;
      PSIs_inv_unst.col(tmburn) = psis_inv;

      P_ACCEPT  = (tmburn*P_ACCEPT + ACCEPT)/(tmburn + 1);

    }
  }

  return Rcpp::List::create(Rcpp::Named("LAMBDA",LAMBDA),
                            Rcpp::Named("PSIs",PSIs),
                            Rcpp::Named("ROW_OUT", ROW_OUT),
                            Rcpp::Named("THRESHOLDS", THRESHOLDS_OUT),
                            Rcpp::Named("INTERCEPTS", INTERCEPT_OUT),
                            Rcpp::Named("ACCEPTED", ACCEPTED),
                            Rcpp::Named("MHACCEPT",P_ACCEPT),
                            Rcpp::Named("LAMBDA_unst",LAMBDA_unst),
                            Rcpp::Named("PSIs_inv_unst",PSIs_inv_unst),
                            Rcpp::Named("THRESHOLDS_unst", THRESHOLDS_unst),
                            Rcpp::Named("INTERCEPTS_unst", INTERCEPT_unst)
                              );
}

//' @title Sample Observed Data Given Augmented Data
//' @description Draw observed data given augmented values for the general mixed-type response data for binary, ordinal, and bounded variables.
//' @param threshold A vector of category thresholds.
//' @param Msj The model indicator variable.
//' @param MISS A N by J matrix of missing data indicators.
//' @param Z A scalar augmented value.
//' @param bounds The lower and upper bounds (only used if Msj = 0).
//' @author Steven Andrew Culpepper
//' 
//' @noRd
// [[Rcpp::export]]
double sampleY_given_Z(const arma::vec& threshold,const double& Msj,const double& Z,const arma::vec& bounds){
  
  double Y=0.;
  if(Msj>1){
    Y=0;
    unsigned int thres_index=1;
    while( (Z>threshold(thres_index))&&(thres_index<Msj) ){
      Y+=1.;
      thres_index+=1;
    }
  }
  if(Msj==0){
    Y = Z;
    if(Z<bounds(0)){
      Y = bounds(0);
    }
    if(Z>bounds(1)){
      Y = bounds(1);
    }
  }
  return Y;
}


//' @title Posterior Predictions of Item Factor Analysis with Mixed Response Types
//' @description Generate posterior predictions for new variables using posterior samples.
//' @param OUTPUT A list of output from IFA_Mode_Jumper_MixedResponses.
//' @param Y A N by J matrix of item responses for predictions. Variables to predict are indicated in Y by NAs.
//' @param Ms model indicator where 0 = "bounded", 1 = "continuous", 2 = "binary", >2 = "ordinal".
//' @param variable_predict_flag A J vector. 0 = do not predict the variable; 1 = predict the variable.
//' @param bounds A J by 2 matrix denoting the min and max variable values. Note that bounds are only used for variable j if element j of Ms is zero.
//' @param n_mcmc_iterations The number of Gibbs iterations for sampling posterior predictions for factor scores and missing data. The default is 10 iterations.
//' @return array of predictions for all posterior samples provided in OUTPUT.
//' @author Steven Andrew Culpepper
//' @export
// [[Rcpp::export]]
arma::cube mixedresponse_posterior_prediction(const Rcpp::List& OUTPUT,
                                              const arma::mat& Y,
                                              const arma::vec& Ms,
                                              const arma::vec& variable_predict_flag,
                                              const arma::mat& bounds,
                                              unsigned n_mcmc_iterations=10){
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  arma::mat MISS = arma::ones<arma::mat>(N,J);
  MISS.elem( arma::find_nonfinite(Y) ).fill(0.0);
  
  //extracting samples
  arma::cube LAMBDA_unst = OUTPUT["LAMBDA_unst"];
  arma::mat PSIs_inv_unst = OUTPUT["PSIs_inv_unst"];
  arma::cube THRESHOLDS_unst =  OUTPUT["THRESHOLDS_unst"];
  arma::mat INTERCEPTS_unst =  OUTPUT["INTERCEPTS_unst"];
  
  arma::uvec which_vars_to_predict = find(variable_predict_flag==1);
  unsigned int n_vars_to_predict = which_vars_to_predict.n_elem;
  unsigned int n_mcmc_samples = LAMBDA_unst.n_slices;
  unsigned int M = LAMBDA_unst.n_cols;
  arma::mat I_K = arma::eye(M,M);
  
  arma::cube predicted_Y=arma::zeros<arma::cube>(N,n_vars_to_predict,n_mcmc_samples);
  
  //need a function to generate Y|Z,
  
  for(unsigned int k=0;k<n_mcmc_samples;k++){
    //extract parameters for iteration k
    arma::mat Lambda = LAMBDA_unst.slice(k);
    arma::vec psis_inv = PSIs_inv_unst.col(k);
    arma::mat threshold = THRESHOLDS_unst.slice(k);
    arma::vec intercept = INTERCEPTS_unst.col(k);
    
    arma::mat F=arma::zeros<arma::mat>(N,M);
    arma::mat Z(N,J);
    arma::mat Z_centered(N,J);
    
    arma::vec pky_new(2);
    
    //for loop to iterate between F samples and Z samples
    for(unsigned int Gibbs_its=0;Gibbs_its<n_mcmc_iterations;Gibbs_its++){
      
      //sample all Z|F and parameters
      //loop over i and j
      for(unsigned int j=0;j<J;j++){
        arma::vec eta = intercept(j)+F*(Lambda.row(j)).t();
        if(Ms(j)>1){
          for(unsigned int i=0;i<N;i++){
            double uZ = R::runif(0.0,1.0);
            
            if(MISS(i,j)==1.0){
              pky_new(0) = R::pnorm(threshold(Y(i,j)    ,j)  ,eta(i),1.0, 1, 0);
              pky_new(1) = R::pnorm(threshold(Y(i,j)+1.0,j)  ,eta(i),1.0, 1, 0);
              Z(i,j) = R::qnorm(pky_new(0) + uZ*(pky_new(1)-pky_new(0) ),eta(i),1.0,1,0);
            }
            if(MISS(i,j)==0){
              Z(i,j) = R::rnorm(eta(i),1.0);
            }
          }
          
        }
        
        if(Ms(j)==1){
          double sqrt_psi = sqrt(1./psis_inv(j));
          for(unsigned int i=0;i<N;i++){
            if(MISS(i,j)==1){
              Z(i,j) = Y(i,j);
            }
            if(MISS(i,j)==0){
              Z(i,j) = R::rnorm(eta(i),sqrt_psi) ;
            }
          }
        }
        if(Ms(j)==0){
          double sqrt_psi = sqrt(1./psis_inv(j));
          for(unsigned int i=0;i<N;i++){
            if(MISS(i,j)==1.0){
              double Yij = Y(i,j);
              Z(i,j) = Yij;
              if(Yij==bounds(j,0)){
                //sample truncated above at bound(j,0)
                Z(i,j) = rTruncNorm_bounds(eta(i),sqrt_psi,0,bounds(j,0));
              }
              if(Yij==bounds(j,1)){
                //sample truncated below at bound(j,1)
                Z(i,j) = rTruncNorm_bounds(eta(i),sqrt_psi,1,bounds(j,1));
              }
            }
            if(MISS(i,j)==0){
              Z(i,j) = R::rnorm(eta(i),sqrt_psi) ;
            }
          }
        }
      }
      //sample F|Z and parameters
      Z_centered = Z - repmat(intercept, 1, N).t();
      update_F_matrix(Z_centered,I_K,F,Lambda,psis_inv);
    }
    //sample Y|Z
    for(unsigned int jpred=0;jpred<n_vars_to_predict;jpred++){
      unsigned int j = which_vars_to_predict(jpred);
      for(unsigned int i=0;i<N;i++){
        if(Ms(j)>1){
          double Zij = Z(i,j);
          predicted_Y(i,jpred,k) = sampleY_given_Z(threshold.col(j),Ms(j),Zij,(bounds.row(j)).t() );
        }
        if(Ms(j)==1){
          predicted_Y(i,jpred,k)=Z(i,j);
        }
        if(Ms(j)==0){
          double Zij = Z(i,j);
          predicted_Y(i,jpred,k) = sampleY_given_Z(threshold.col(j),Ms(j),Zij,(bounds.row(j)).t() );
        }
      }
    }
  }
  return predicted_Y;
}

