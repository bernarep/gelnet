// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include "Rcpp.h"
// using namespace arma;  // use the Armadillo library for matrix computations
using namespace Rcpp;

// [[Rcpp::export(stage2solver)]]
List stage2solver(arma::mat S, arma::mat Winit, arma::mat E0, double toler = 1e-5, int maxit = 500){
  
  unsigned int d = S.n_rows;
  arma::mat E0c(d,d);
  E0c.ones();
  E0c = E0c - E0;
  arma::mat W = Winit;
  int matrix_iter = 0;
  bool ConvThr = false;
  arma::mat W_old = W;
  arma::mat W11;
  arma::mat s12;
  arma::mat gamma12;
  arma::mat W11s;
  arma::mat s12s;
  arma::mat b_opt_s;
  arma::mat b_opt((d-1),1);
  arma::mat w12((d-1),1);
  arma::uvec b_idx;
  arma::uvec tmp_idx;
  arma::uvec uj;
  arma::uvec rel_vars;
  arma::vec den;
  
  
  while( (ConvThr == false) && (matrix_iter <= maxit) ){
    W_old = W;
    for(unsigned int j = 0; j < d; j++){
      uj = {j};
      tmp_idx = arma::regspace<arma::uvec>(0,  1,  (d-1));
      tmp_idx.shed_row(j);
      b_idx = find(E0c(tmp_idx,uj) == 0);
      if(b_idx.n_elem != 0){
        W11 = W(tmp_idx, tmp_idx);
        gamma12 = E0c(tmp_idx,uj);
        s12 = S(tmp_idx,uj);
        rel_vars = find(gamma12 == 0);
        W11s = W11(rel_vars, rel_vars);
        s12s = s12.rows(rel_vars);
        b_opt_s = inv(W11s)*s12s;
        b_opt.zeros();
        b_opt.rows(rel_vars) = b_opt_s;
        w12 = W11*b_opt;
      }else{
        w12.zeros();
      }
      W(tmp_idx, uj) = w12;
      W(uj, tmp_idx) = trans(w12);
    }
    
    matrix_iter = matrix_iter + 1;
    
    if(abs(W_old - W).max() <= toler){
      ConvThr = true;
    }
    
  }
  
  arma::mat O = S;
  
  for(unsigned int j = 0; j < d; j++){
    uj = {j};
    tmp_idx = arma::regspace<arma::uvec>(0,  1,  (d-1));
    tmp_idx.shed_row(j);
    b_idx = find(E0c(tmp_idx,uj) == 0);
    if(b_idx.n_elem != 0){
      W11 = W(tmp_idx, tmp_idx);
      gamma12 = E0c(tmp_idx,uj);
      s12 = S(tmp_idx,uj);
      rel_vars = find(gamma12 == 0);
      W11s = W11(rel_vars, rel_vars);
      s12s = s12.rows(rel_vars);
      b_opt_s = inv(W11s)*s12s;
      b_opt.zeros();
      b_opt.rows(rel_vars) = b_opt_s;
      w12 = W11*b_opt;
    }else{
      b_opt.zeros();
      w12.zeros();
    }
    W(tmp_idx, uj) = w12;
    W(uj, tmp_idx) = trans(w12);
    den = S(uj,uj)-trans(w12)*b_opt;
    O(j,j) = 1.0/den(0);
    O(uj, tmp_idx) = trans(-1*b_opt * O(uj,uj));
    O(tmp_idx, uj) = -1*b_opt * O(uj,uj);
  }
  
  List out;
  out["O"] = O;
  out["W"] = W;
  out["iters"] = matrix_iter;
  return(out);
}

// In order to compile the code you have to run sourceCpp("D:/DatiDavide/Downloads/test_gelasticnet.cpp") in console
// END