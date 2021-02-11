#include <RcppArmadillo.h>
#include <commonEFAfunctions.h>





//function that is unique to IFA
//will need a function to sample for ordered categories too
//' @title Sample truncated normal/ 
//' @description Sample a truncated normal random variable that is bounded below.
//' @param mean Location parameter for truncated normal.
//' @param sd Square root of truncated normal scale parameter.
//' @param b_lb The lower bound.
//' @author Albert X Man, Steven Andrew Culpepper
//' @return A value from the truncated normal distribution.
//' @noRd
// [[Rcpp::export]]
double rTruncNorm_lb(double mean,double sd, double b_lb){
  double p0 = R::pnorm(b_lb,mean,sd, 1, 0);
  double p1 = 1-p0;
  double uZ = R::runif(0,1);
  double Z = R::qnorm(p0 + uZ*p1, mean, sd, 1, 0);
  return(Z);
}



//' @title Sample augmented data
//' @description Sample augmented data for binary and ordinal response data 
//' @param Y A N by J matrix of ordinal responses.
//' @param MISS A N by J matrix of missing data indicators.
//' @param Z A N by J matrix of augmented data.
//' @param as A J by M matrix of loadings.
//' @param bs A J vector of intercepts.
//' @param theta A N by M matrix of factor scores.
//' @param Ms A vector of the number of score categories. Note 2 = two parameter normal ogive, >2 ordinal normal ogive.
//' @param Kaps A matrix of category thresholds.
//' @param sdMH A vector of tuning parameters of length J for implementing the Cowles (1996) Metropolis-Hastings threshold sampler.
//' @author Steven Andrew Culpepper
//' @return A list containing:
//' 
//' - `Kaps_new`: An updated matrix of category thresholds.
//' - `MHaccept`: A binary vector indicating whether an MH proposal for thresholds was accepted. 
//'  
//' @noRd
// [[Rcpp::export]]
Rcpp::List update_WKappaZ_NA(const arma::mat& Y,//binary data
                             const arma::mat MISS,//missing indicators
                             arma::mat& Z,
                             const arma::mat& as,//lambda matrix is called "as" here
                             const arma::vec& bs,//intercept is called "bs" here
                             const arma::mat& theta,//F is called "theta" here
                             const arma::vec& Ms,
                             const arma::mat& Kaps,//thresholds is called "Kaps" here
                             const arma::vec& sdMH) {
  
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
  // double Phi_eta;//,sj,tj,n0dot,n1dot,n01,n10;
  // double us,ug,s_temp,pg_temp,ps_temp,g_temp,pg,ps;
  // arma::vec cs_new=arma::zeros<arma::vec>(J);
  // arma::vec ss_new=arma::zeros<arma::vec>(J);
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
      // if(model(j)=="ONO"){
      // if(Ms(j)==2){
      //   goto Model_is_2PNO;
      // }
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
      // if(model(j)=="2PNO"){
      // Model_is_2PNO:
      // 
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
   
  }
   
  return Rcpp::List::create(
    //Rcpp::Named("pk_on",pk_on),
    //Rcpp::Named("uMH",uMH),
    // Rcpp::Named("KapsMH",KapsMH),
    // Rcpp::Named("mpart",mpart),
    // Rcpp::Named("pky_on",pky_on),
    // Rcpp::Named("ipart",ipart),
    Rcpp::Named("Kaps_new",Kaps_new), 
    Rcpp::Named("MHaccept",MHaccept)  );
}



//unique function to IFA
//' @title Sample Intercept Parameters
//' @description Sample intercept parameters from the posterior distribution.
//' @param intercept A J vector of variable intercepts.
//' @param Z A N by J matrix of augmented data.
//' @param F A N by M matrix of factor scores.
//' @param psis_inv A J vector of inverse-uniquenesses.
//' @param intercept_var A hyperparameter for the scale of item intercepts.
//' @author Albert X Man
//' 
//' @noRd
// [[Rcpp::export]]
void update_intercept(arma::vec& intercept,
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
  intercept_var = pow(R::qgamma(u, double(J + 2.0*b1)/2.0, 2.0/(2.0*b2 + d),1,0), -1);
}



//' @title Compute Deviance for Ordinal Probit Factor Model
//' @description Internal function to compute -2LL using unstandardized parameters.
//' @param N The number of observations.  (> 0)
//' @param J The number of items.  (> 0)
//' @param Y A N by J matrix of item responses.
//' @param MISS A N by J matrix of missing data indicators.
//' @param as A matrix of item loadings.
//' @param bs A vector of threshold parameters.
//' @param theta A matrix of factor scores.
//' @param Ms A vector of the number of score categories. Note 2 = two parameter normal ogive, >2 ordinal normal ogive.
//' @param Kaps A matrix of category thresholds.
//' @return -2LL.
//' @author Steven Andrew Culpepper
//' @export
// [[Rcpp::export]]
double min2LL_ono(unsigned int N,unsigned int J,const arma::mat& Y,const arma::mat& MISS,const arma::mat& as,
                   const arma::vec& bs,const arma::mat& theta,
                   const arma::vec& Ms,const arma::mat& Kaps){
  
  double m2ll=0.0;
  double p_eta,p_eta0,yij;
  
  for(unsigned int j=0;j<J;j++){
    
    arma::vec eta = bs(j)+theta*(as.row(j)).t();
    
    // if(Ms(j)==2){
    //   for(unsigned int i=0;i<N;i++){
    //     
    //     if(MISS(i,j)==1.0){
    //       
    //       // eta=as(j)*theta(i,area(j)-1) - bs(j);
    //       p_eta = R::pnorm(eta(i),0.0,1.0, 1, 0);
    //       yij = Y(i,j);
    //       
    //       if(yij==0.0){
    //         m2ll=-2.0*log(1.0-gs(j)-(1.0-ss(j)-gs(j))*p_eta) + m2ll;
    //       }
    //       
    //       if(yij==1.0){
    //         m2ll=-2.0*log(gs(j)+(1.0-ss(j)-gs(j))*p_eta)+ m2ll;
    //       }        
    //     }
    //   }
    // }
    
    if(Ms(j)>1){
      for(unsigned int i=0;i<N;i++){
        
        if(MISS(i,j)==1.0){
          // eta=as(j)*theta(i,area(j)-1) - bs(j);
          yij = Y(i,j);
          
          for(unsigned int m=0;m<Ms(j);m++){
            if(yij==m){
              p_eta  = R::pnorm(Kaps(yij+1.0,j) - eta(i),0.0,1.0, 1, 0);
              p_eta0 = R::pnorm(Kaps(yij    ,j) - eta(i),0.0,1.0, 1, 0);
              m2ll=-2.0*log(p_eta-p_eta0) + m2ll;
            }
          }
        }
      }
    }
  }
  return m2ll;
}



//' @title Exploratory Item Factor Analysis for Ordinal Response Data
//' @description Implement the Man and Culpepper (2020) mode-jumping algorithm to factor analyze ordinal response data. Missing values should be specified as a non-numeric value such as NA.
//' @param Y A N by J matrix of item responses.
//' @param M The number of factors.
//' @param gamma Mode-jumping tuning parameter. Man & Culpepper used 0.5.
//' @param Ms A vector of the number of score categories. Note 2 = two parameter normal ogive, >2 ordinal normal ogive.
//' @param sdMH A vector of tuning parameters of length J for implementing the Cowles (1996) Metropolis-Hastings threshold sampler.
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
//' 
//' @author Albert X Man, Steven Andrew Culpepper
//' @references
//' Cowles, M. K. (1996), Accelerating Monte Carlo Markov chain convergence for cumulative link generalized linear models," Statistics and Computing, 6, 101-111.
//' 
//' Man, A. & Culpepper, S. A. (2020). A mode-jumping algorithm for Bayesian factor analysis. Journal of the American Statistical Association, doi:10.1080/01621459.2020.1773833.
//' 
//' @examples
//' 
//' #Load the psych package to apply the algorithm with the bfi dataset
//' library(psych)
//' data(bfi)
//' 
//' #subtract 1 from all items so scores start at 0
//' ydat<-as.matrix(bfi[,1:25])-1
//' J<-ncol(ydat)
//' N<-nrow(ydat)
//' 
//' #compute the number of item categories
//' Ms<-apply(ydat,2,max)+1
//' 
//' #specify burn-in and chain length
//' #in full application set these to 10000 and 20000
//' burnin = 10 
//' chain_length=20
//' 
//' out5<-IFA_Mode_Jumper(ydat,M=5,gamma=.5,Ms,sdMH=rep(.025,J),burnin,chain_length)
//' 
//' #check the acceptance rates for the threshold tuning parameter fall between .2 and .5
//' mh_acceptance_rates<-t(out5$MHACCEPT)
//' 
//' #compute mean thresholds
//' mthresholds<-apply(out5$THRESHOLDS,c(1,2),mean)
//' 
//' #compute mean intercepts
//' mintercepts<-apply(out5$INTERCEPTS,1,mean)
//' 
//' #rotate loadings to PLT
//' my_lambda_sample = out5$LAMBDA
//' my_M<-5
//' for (j in 1:dim(my_lambda_sample)[3]) {
//'   my_rotate = proposal2(1:my_M,my_lambda_sample[,,j],matrix(0,nrow=N,ncol = my_M))
//'   my_lambda_sample[,,j] = my_rotate$lambda
//' }
//' 
//' #compute mean of PLT loadings
//' mLambda<-apply(my_lambda_sample,c(1,2),mean)
//' 
//' #find promax rotation of posterior mean
//' rotatedLambda<-promax(mLambda)$loadings[1:J,1:my_M]
//' 
//' #save promax rotation matrix
//' promaxrotation<-promax(mLambda)$rotmat
//' 
//' #compute the factor correlation matrix
//' phi<-solve(t(promaxrotation)%*%promaxrotation)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List IFA_Mode_Jumper(const arma::mat& Y,
                           unsigned int M,double gamma,
                           const arma::vec& Ms,const arma::vec& sdMH,
                           unsigned int burnin,unsigned int chain_length=10000){
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  double my_gamma = gamma;
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
  
  // Initialize r_idx as PLT
  arma::uvec r_idx;
  arma::uvec one_to_J = arma::linspace<arma::uvec>(0, J-1, J);
  r_idx = sort(one_to_J.subvec(0, M-1));
  
  arma::mat thresholds = kappa_initialize(Ms);
  arma::mat F = arma::zeros<arma::mat>(N,M);//arma::randn<arma::mat>(N,M);
  arma::mat Cstart=inv(I_K+F.t()*F);
  
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
  }
  
  
  arma::vec psis_inv(J, arma::fill::ones);
  arma::rowvec varY=arma::var(Y);
  for(unsigned int j=0;j<J;++j){
    arma::rowvec lambdaj = Lambda.row(j);
    arma::vec Y_m_Flam = Y.col(j)-F*lambdaj.t();
  }
  
  // arma::vec thresholds(J, arma::fill::zeros);
  arma::vec intercept(J, arma::fill::zeros);
  double intercept_var = 1;
  
  arma::mat invClam(J, M, arma::fill::ones);
  // invClam.fill(1e-5);
  
  // Main Algorithm Loop
  double accepted;
  // arma::vec metropolis_probs;
  Rcpp::List L;
  // double marginal_like;
  
  
  arma::mat Z_centered = arma::zeros<arma::mat>(N, J);
  arma::vec marginal_like_vec;
  
  // std::cout << "begin\n";
  
  //save MH acceptance for thresholds
  arma::vec P_ACCEPT= arma::zeros<arma::vec>(J);
  arma::vec ACCEPT= arma::zeros<arma::vec>(J);
  
  arma::uvec thresholdindices = arma::linspace<arma::uvec>(1, max(Ms)-1, max(Ms)-1);
  
  for(unsigned int t = 0; t < chain_length; ++t){
    accepted = 0;
    
    
    //thresholds need to be in the right format for this function
    //the first and last thresholds are -inf and inf
    //thresholds for items with fewer categories are nan
    Rcpp::List step1Z = update_WKappaZ_NA(Y,//binary data
                      MISS,//missing indicators
                      Z,
                      Lambda,//lambda matrix is called "as" here
                      intercept,//intercept is called "bs" here
                      F,//F is called "theta" here
                      Ms,
                      thresholds,//thresholds is called "Kaps" here
                      sdMH);
    
    thresholds  = Rcpp::as<arma::mat>(step1Z[0]);
    ACCEPT= Rcpp::as<arma::vec>(step1Z[1]);
    update_intercept(intercept, Z, F, Lambda, psis_inv, intercept_var);

    Z_centered = Z - repmat(intercept, 1, N).t();
    update_F_matrix(Z_centered, I_K, F, Lambda, psis_inv);
    update_Lambda_loadings_hard_zero(Z_centered, r_idx, F, Lambda, psis_inv, invClam, my_gamma);

    
    // std::cout << psis_inv(1) << "\n";
    update_invClam(Lambda, invClam);
    

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
    
    // double temp6 = L["marginal_like"];
    // marginal_like = temp6;
    
    /*    */
    
    // Save output after burn-in
    if(t>(burnin-1)){
      
      arma::mat A1 = diagmat(pow(diagvec(Lambda * Lambda.t()) + pow(psis_inv, -1), -0.5));
      // arma::mat A2 = diagmat(pow(diagvec(Lambda * Lambda.t()) + pow(psis_inv, -1), -1));
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
      // LIKELIHOOD(tmburn) = marginal_like;
      // MH_PROBS.col(tmburn) = metropolis_probs;
      
      P_ACCEPT  = (tmburn*P_ACCEPT + ACCEPT)/(tmburn + 1);
      /*  */      
    }
  }
  return Rcpp::List::create(Rcpp::Named("LAMBDA",LAMBDA),
                            // Rcpp::Named("PSIs",PSIs),
                            Rcpp::Named("ROW_OUT", ROW_OUT),
                            // Rcpp::Named("F_OUT", F_OUT),
                            Rcpp::Named("THRESHOLDS", THRESHOLDS_OUT),
                            Rcpp::Named("INTERCEPTS", INTERCEPT_OUT),
                            // Rcpp::Named("MH_PROBS", MH_PROBS),
                            Rcpp::Named("ACCEPTED", ACCEPTED),
                            Rcpp::Named("MHACCEPT",P_ACCEPT)//,
                            // Rcpp::Named("LIKELIHOOD", LIKELIHOOD)
                              );
}










