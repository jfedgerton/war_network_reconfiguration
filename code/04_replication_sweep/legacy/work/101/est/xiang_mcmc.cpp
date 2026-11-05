// Authors' dyadEffect sampler (Xiang & Ling 2016), unchanged except: progress/omp/String.h includes removed, header inlined,
// and an optional alpha_init argument so a chain can be continued from its saved state (segmented runs under a shared lock).
#include <RcppArmadillo.h>
#include <ctime>
#include <vector>
#include <string>
#ifndef MCMC_DYADEFFECT_H
#define MCMC_DYADEFFECT_H



double rtbnorm_combo(const double& m, const double& v, const double& below,
                     const int& iter=10);


void MCMCdyadeffect(
    const arma::mat& Y, const arma::mat& G, const arma::mat& X,
    unsigned int mcmc, unsigned int burnin, unsigned int thin,
    unsigned long seed, unsigned int verbose,
    double alpha_mean0, double alpha_var0, arma::vec& gamma0, arma::vec& beta0,
    arma::mat& alpha_store, arma::mat& alpha_param_store,
    arma::mat& gamma_store, arma::mat& beta_store,
    int cores, const arma::vec& alpha_init);

#endif





// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export(name = "MCMCdyadeffect.matrix")]]
Rcpp::NumericMatrix MCMCdyadeffect_matrix(
    const Rcpp::NumericVector& Y, const Rcpp::NumericVector& G, const Rcpp::NumericVector& X,
    unsigned int mcmc = 1000000, unsigned int burnin = 100000, unsigned int thin = 500,
    unsigned int seed = NA_INTEGER, unsigned int verbose = 0, Rcpp::String outfile = NA_STRING,
    double alpha_mean0 = 0, double alpha_var0 = 1.0,
    Rcpp::NumericVector gamma0 = Rcpp::NumericVector::create(0),
    Rcpp::NumericVector beta0 = Rcpp::NumericVector::create(0),
    int cores = 1,
    Rcpp::NumericVector alpha_init = Rcpp::NumericVector(0)
  ){

   int n_dyad;
   int n_time;
   int n_rel;
   int n_conf;

   {
     Rcpp::IntegerVector d = Y.attr("dim");
     int n_d = std::distance(d.begin(), d.end());
     if(n_d != 2){
       Rcpp::Rcerr << "Error: Y is not a 2-d array." << std::endl;
       Rcpp::Rcerr << "       dim(Y) should be c(n_dyad, n_time)." << std::endl;
       Rcpp::stop("");
     }
     n_dyad = d[0];
     n_time = d[1];
   }

   {
     Rcpp::IntegerVector d = G.attr("dim");
     int n_d = std::distance(d.begin(), d.end());
     if(n_d != 3){
       Rcpp::Rcerr << "Error: G is not a 3-d array." << std::endl;
       Rcpp::Rcerr << "       dim(G) should be c(n_dyad, n_time, n_rel)." << std::endl;
       Rcpp::stop("");
     }

    int n_dyad_rel = d[0];
    int n_time_rel = d[1];
    if(n_dyad_rel != n_dyad || n_time_rel != n_time){
      Rcpp::Rcerr << "Error: dimension not match." << std::endl;
      Rcpp::Rcerr << "       Y is " << n_dyad << "-by-" << n_time << "." << std::endl;
      Rcpp::Rcerr << "       Each slice of G is " << n_dyad_rel << "-by-" << n_time_rel << "." << std::endl;
      Rcpp::stop("");
    }

    n_rel = d[2];
  }

  {
    Rcpp::IntegerVector d = X.attr("dim");
    int n_d = std::distance(d.begin(), d.end());
    if(n_d != 3){
      Rcpp::Rcerr << "Error: X is not a 3-d array." << std::endl;
      Rcpp::Rcerr << "       dim(X) should be c(n_dyad, n_time, n_conf)." << std::endl;
      Rcpp::stop("");
    }

    int n_dyad_conf = d[0];
    int n_time_conf = d[1];
    if(n_dyad_conf != n_dyad || n_time_conf != n_time){
      Rcpp::Rcerr << "Error: dimension not match." << std::endl;
      Rcpp::Rcerr << "       Y is " << n_dyad << "-by-" << n_time << "." << std::endl;
      Rcpp::Rcerr << "       Each slice of X is " << n_dyad_conf << "-by-" << n_time_conf << "." << std::endl;
      Rcpp::stop("");
    }

    n_conf = d[2];
  }

  arma::mat Y_smat(&Y[0], n_dyad, n_time);
  arma::mat G_smat(&G[0], n_dyad*n_time, n_rel);
  arma::mat X_smat(&X[0], n_dyad*n_time, n_conf);
  X_smat = join_rows(X_smat, arma::ones<arma::mat>(n_dyad*n_time, 1));
  n_conf += 1;

  if(seed == NA_INTEGER){
    seed = std::time(NULL);
  }

  if(gamma0.size() == 1 && gamma0[0] == 0)
    gamma0 = Rcpp::NumericVector(n_rel, 0.1);
  arma::vec gamma0_smat(&gamma0[0], n_rel);

  if(beta0.size() == 1 && beta0[0] == 0)
    beta0 = Rcpp::NumericVector(n_conf, 0.1);
  arma::vec beta0_smat(&beta0[0], n_conf);

  arma::mat alpha_store(n_dyad, mcmc/thin + 1);
  arma::mat alpha_param_store(2, mcmc/thin + 1);
  arma::mat gamma_store(n_rel, mcmc/thin + 1);
  arma::mat beta_store(n_conf, mcmc/thin + 1);

  arma::vec alpha_init_v(alpha_init.begin(), alpha_init.size());
  MCMCdyadeffect(Y_smat, G_smat, X_smat, mcmc, burnin, thin, seed, verbose,
                 alpha_mean0, alpha_var0, gamma0_smat, beta0_smat,
                 alpha_store, alpha_param_store, gamma_store, beta_store, cores, alpha_init_v);

  arma::mat param_smat = arma::join_cols(alpha_param_store, gamma_store);
  param_smat = arma::join_cols(param_smat, beta_store);
  param_smat = param_smat.t();
  
  alpha_store = alpha_store.t();


  Rcpp::String outfile_param(outfile);
  outfile_param += ".txt";
  Rcpp::String outfile_alpha(outfile);
  outfile_alpha += "_alpha.txt";


  if(outfile != NA_STRING){
    param_smat.save(outfile_param, arma::raw_ascii);
    alpha_store.save(outfile_alpha, arma::raw_ascii);
  }

  arma::mat out_smat = arma::join_rows(alpha_store, param_smat);

  Rcpp::NumericMatrix out(mcmc/thin + 1, (n_dyad + 2 + n_rel + n_conf), out_smat.begin());

  return out;
}

void MCMCdyadeffect(
    const arma::mat& Y, const arma::mat& G, const arma::mat& X,
    unsigned int mcmc, unsigned int burnin, unsigned int thin,
    unsigned long seed, unsigned int verbose,
    double alpha_mean0, double alpha_var0, arma::vec& gamma0, arma::vec& beta0,
    arma::mat& alpha_store, arma::mat& alpha_param_store,
    arma::mat& gamma_store, arma::mat& beta_store,
    int cores, const arma::vec& alpha_init){

  (void)cores; // omp removed (single core, as authors: n_cores <- 1)
  
  unsigned int n_iter = mcmc + burnin;
  unsigned int n_dyad = Y.n_rows;
  unsigned int n_time = Y.n_cols;
  unsigned int n_rel  = G.n_cols;
  unsigned int n_conf = X.n_cols;

  // create some index to label missing data
  arma::mat Y_ID(n_dyad, n_time);
  arma::mat G_NEW = G;
  arma::mat X_NEW = X;

  int count_missing = 0;
  for(unsigned int j = 0; j < n_time; j++){
    for(unsigned int i = 0; i < n_dyad; i++){
      if(Y(i,j) < 0){ // Y(d,t) unknown
        G_NEW.row(count_missing).fill(0);
        X_NEW.row(count_missing).fill(0);
        Y_ID(i,j) = 0;
      }else{
        Y_ID(i,j) = 1; // Y(d,t) known
      }
      count_missing++;
    }
  }

  arma::vec n_time_new = sum(Y_ID, 1);


  ///////////////////////////////
  // Initialization of Sampler //
  ///////////////////////////////

  // initialize relevance coefficients
  arma::vec b0_rel(n_rel, arma::fill::zeros);
  arma::mat B0_rel = .001 * arma::eye<arma::mat>(n_rel, n_rel);
  arma::vec B0b0_rel = B0_rel * b0_rel;

  arma::vec Zstar(n_dyad*n_time);
  double alpha_mean = alpha_mean0;
  double alpha_var = alpha_var0;
  arma::vec alpha = sqrt(alpha_var) * arma::randn<arma::vec>(n_dyad) + alpha_mean;
  if (alpha_init.n_elem == n_dyad) alpha = alpha_init;   // ADDED (estaudit): continue a chain from its saved dyad effects
  arma::mat alpha_extend = arma::kron(arma::vec(n_time, arma::fill::ones), alpha);

  arma::vec gamma = gamma0;
  arma::mat gamma_var = arma::inv_sympd(G_NEW.t() * G_NEW + B0_rel);
  arma::mat gamma_var_chol = (chol(gamma_var)).t();


  // initialize conflict coefficients
  arma::vec b0_conf(n_conf, arma::fill::zeros);
  arma::mat B0_conf = .001 * arma::eye<arma::mat>(n_conf, n_conf);
  arma::vec B0b0_conf = B0_conf * b0_conf;

  arma::vec Ystar(n_dyad*n_time);
  arma::vec beta = beta0;
  arma::mat beta_var = arma::inv_sympd(X_NEW.t() * X_NEW + B0_conf);
  arma::mat beta_var_chol = (chol(beta_var)).t();



  ///////////////////
  // Gibbs Sampler //
  ///////////////////
  for(unsigned int iter = 0; iter <= n_iter; iter++){
    if(iter % 100 == 0){
      Rcpp::checkUserInterrupt();
    }

    // sample Z
    for(unsigned int j = 0; j < n_time; j++){
      for(unsigned int i = 0; i < n_dyad; i++){
        int row = i + j*n_dyad;
        if(Y(i,j) < 0){ // Y(i,j) unknown
          Zstar(row) = 0;
          Ystar(row) = 0;
        }else{ // Y(i,j) known
          double Zstar_mean = dot(G.row(row),gamma) + alpha(i);
          double Ystar_mean = dot(X.row(row),beta);

          if(Y(i,j) == 1){ // Zstar >= 0, Ystar >= 0
            Zstar(row) = rtbnorm_combo(Zstar_mean, 1.0, 0, 10);
            Ystar(row) = rtbnorm_combo(Ystar_mean, 1.0, 0, 10);
          }else{ // Y(i,j) == 0
            Zstar(row) = R::rnorm (Zstar_mean, 1.0);
			Ystar(row) = R::rnorm (Ystar_mean, 1.0);
			while (Zstar(row)>= 0.0 && Ystar(row)>= 0.0){
				   Zstar(row) = R::rnorm (Zstar_mean, 1.0);
				   Ystar(row) = R::rnorm (Ystar_mean, 1.0);
			}		
          }
        }
      }
    }

    // sample alpha_mean
    double alpha_mean_mean = mean(alpha);
    double alpha_mean_var  = alpha_var / n_dyad;
    alpha_mean = R::rnorm(alpha_mean_mean, sqrt(alpha_mean_var));

    // sample alpha_var
    double alpha_var_a = n_dyad/2 - 1;
    double alpha_var_b = 0.5 * sum(pow((alpha - alpha_mean), 2));
    alpha_var = 1 / R::rgamma(alpha_var_a, 1/alpha_var_b);

    // sample alpha
    {
      arma::mat tmp = Zstar - G*gamma;
      tmp = reshape(tmp, n_dyad, n_time) % Y_ID;

      arma::vec alpha_post_var = alpha_var/(1+n_time_new*alpha_var);
      arma::vec alpha_post_mean = alpha_post_var % (sum(tmp, 1) + alpha_mean/alpha_var);
      alpha = sqrt(alpha_post_var) % arma::randn(n_dyad) + alpha_post_mean;
      alpha_extend = kron(arma::mat(n_time, 1, arma::fill::ones), alpha);
    }

    // sample gamma, beta
      arma::vec gamma_mean = gamma_var * (G_NEW.t()*(Zstar - alpha_extend) + B0b0_rel);
      arma::vec beta_mean = beta_var * (X_NEW.t()*Ystar + B0b0_conf);

      gamma = gamma_var_chol * arma::randn(n_rel) + gamma_mean;
      beta = beta_var_chol * arma::randn(n_conf) + beta_mean;


    // print results to screen
    if(verbose > 0 && (iter % verbose) == 0){
      Rcpp::Rcout << "Gibbs Iteration = " << iter << std::endl;
      Rcpp::Rcout << "gamma = " << gamma.t();
      Rcpp::Rcout << "beta = " << beta.t();
      Rcpp::Rcout << "alpha_mean = " << alpha_mean << " ~ N(" << alpha_mean_mean << ", " << alpha_mean_var << ")" << std::endl;
      Rcpp::Rcout << "alpha_var = " << alpha_var << " ~ InvG(" << alpha_var_a << ", " << 1/alpha_var_b << ")" << std::endl;
      Rcpp::Rcout << std::endl;
    }

    // store results
    if((iter >= burnin) && ((iter-burnin) % thin) == 0){
      int col = (iter - burnin) / thin;
      alpha_store.col(col) = alpha;
      alpha_param_store(0,col) = alpha_mean;
      alpha_param_store(1,col) = alpha_var;
      gamma_store.col(col) = gamma;
      beta_store.col(col) = beta;
    }

  } // end Gibbs loop


  ////////////////////////
  // Write Output Files //
  ////////////////////////

  if(verbose > 0){
    arma::mat est_gamma = arma::join_rows(arma::mean(gamma_store, 1), arma::stddev(gamma_store, 0, 1));
    arma::mat est_beta = arma::join_rows(arma::mean(beta_store, 1), arma::stddev(beta_store, 0, 1));
    arma::mat est_alpha_param = arma::join_rows(arma::mean(alpha_param_store, 1), arma::stddev(alpha_param_store, 0 ,1));

    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "@@@@@@@@@@@@@@@@@@@@@" << std::endl;
    Rcpp::Rcout << "  Posterior Summary  " << std::endl;
    Rcpp::Rcout << "@@@@@@@@@@@@@@@@@@@@@" << std::endl;
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "Variable\t PostMean\t PostStdDev" << std::endl;
    for(unsigned int i = 0; i < n_rel; i++)
      Rcpp::Rcout << "gamma[" << i+1 << "]\t " << est_gamma(i,0) << "\t " << est_gamma(i,1) << std::endl;
    for(unsigned int i = 0; i < n_conf; i++)
      Rcpp::Rcout << "beta[" << i+1 << "]\t\t " << est_beta(i,0) << "\t " << est_beta(i,1) << std::endl;
    Rcpp::Rcout << "mu_alpha\t " << est_alpha_param(0,0) << "\t " << est_alpha_param(0,1) << std::endl;
    Rcpp::Rcout << "var_alpha\t " << est_alpha_param(1,0) << "\t " << est_alpha_param(1,1) << std::endl;
  }

  return;
}



//////////////////////////////////////////////////////////////////////
//rtbnorm_combo is taken from Martin and Quinn's Scythe library 0.2 //
//////////////////////////////////////////////////////////////////////

double rtbnorm_combo(const double& m, const double& v, const double& below,
                     const int& iter){
  // do rejection sampling and return value
  double sd = sqrt(v);
  if (m >= below){
    double x = R::rnorm(m, sd);
    while (x < below)
      x = R::rnorm(m, sd);
    return x;
  }else{ // do slice sampling and return value
    double z = 0;
    double x = below + .00001;
    for (int i=0; i<iter; ++i){
      z = R::runif(0,1)*exp(-1*pow((x-m),2)/(2*v));
      x = R::runif(0,1)*( (m + sqrt(-2*v*log(z))) - below) + below;
    }
    if (std::isinf(x)){
      return below;
    }
    return x;
  }
}




