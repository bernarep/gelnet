// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
using namespace Rcpp;
// [[Rcpp::export]]
List gelnet_mio_cpp(arma::mat sampleCov, arma::mat Winit, arma::mat Oinit, double alpha, double lambda, double toler, int maxiters, bool penalizediagonal)
{
  // Set initial values and additional ancillary variables
  arma::mat S = sampleCov;
  unsigned int d = S.n_rows;
  // Convergence of block descent algorithm
  bool convBD = false;
  // Convergence of the inner cycle of coordinate descents
  bool convCD;
  // Number of iterations for block descent algorithm and internal coordinate descent
  int itersBD = 0;
  int itersCD = 0;
  // Max numbers of iterations
  int maxitersBD = maxiters;
  int maxitersCD = 500;
  // Numerical convergence thresholds
  double tolerBD = toler;
  double tolerCD = 1e-5;
  // Additional combined penalized parameters
  double l1 = (1-alpha)*lambda;
  double l2 = alpha*lambda;
  // Matrix to store estimated b
  arma::mat B((d-1), d, arma::fill::zeros);
  // Ancillary additional variables
  arma::mat W11;
  arma::mat s12;
  arma::mat w12;
  arma::mat o12;
  double o22;
  double w22;
  arma::uvec miss_index_out;
  arma::uvec miss_index_in;
  arma::uvec ui;
  arma::uvec uj;
  arma::mat b_old((d-1),1);
  arma::vec cj;
  arma::vec hj;
  arma::vec den;
  // Set initial values
  // arma::mat O = arma::diagmat(S + l2*arma::eye(d,d));
  // arma::mat O = arma::inv(arma::diagmat(S + l2*arma::eye(d,d)));
  // arma::mat W = S + l2*arma::eye(d,d) + 2*l1*O;
  arma::mat O = Oinit;
  arma::mat W = Winit;
  arma::mat W_old = W;
  
  // Main cycle (BD)
  while((convBD == false) && (itersBD <= maxitersBD))
  {
    W_old = W;
    for(unsigned int i = 0; i <= (d-1); i++)
    {
      ui = {i};
      miss_index_out = arma::regspace<arma::uvec>(0,  1,  (d-1));
      miss_index_out.shed_row(i);
      W11 = W(miss_index_out,miss_index_out);
      s12 = S(miss_index_out,ui);
      o22 = O(i,i);
      w22 = W(i,i);
      
      itersCD = 0;
      convCD = false;
      while((convCD == false) && (itersCD <= maxitersCD))
      {
        b_old = B.col(i);
        for(unsigned int j = 0; j <= (d-2); j++)
        {
          uj = {j};
          miss_index_in = arma::regspace<arma::uvec>(0,  1,  (d-2));
          miss_index_in.shed_row(j);
          cj = W11(uj, miss_index_in)*B(miss_index_in, ui) - s12(j);
          hj = W11(j,j);
          if(cj(0) > l2)
          {
            B(j,i) = (-cj(0) + l2)/(hj(0) + 2*o22*l1);
          }else if(cj(0) < -l2){
            B(j,i) = (-cj(0) - l2)/(hj(0) + 2*o22*l1);
          }else{
            B(j,i) = 0;
          }
        }
        // Check convergence CD (inner cycle)
        if( max(abs(b_old - B.col(i))) <= tolerCD )
        {
          convCD = true;
        }
        itersCD = itersCD+1;
      }
      // Update values of matrices
      w12 = W11*B.col(i);
      den = trans(B.col(i))*w12;
      o22 = 1.0/(w22 - den(0));
      o12 = -B.col(i)*o22;
      if(penalizediagonal==true){
        w22 = S(i,i) + l2 + 2*l1*o22;
      }else{
        w22 = S(i,i);
      }
      O(miss_index_out,ui) = o12;
      O(ui,miss_index_out) = trans(o12);
      O(i,i)  = o22;
      W(miss_index_out,ui) = w12;
      W(ui,miss_index_out) = trans(w12);
      W(i,i) = w22;
    }
    // Check convergence BD (outer cycle)
    if(abs(W_old-W).max() <= tolerBD)
    {
      convBD = true;
    }
    itersBD = itersBD+1;
  }
  // Create output of function
  List out;
  out["O"] = O;
  out["W"] = W;
  out["iters"] = itersBD;
  //out["iters"] = abs(W_old-W).max();
  return(out);
  
  
}
