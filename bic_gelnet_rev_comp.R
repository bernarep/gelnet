# Load functions for algorithm 1
library(Rcpp)
library(RcppArmadillo)
sourceCpp("gelnet_mio_cpp.cpp")

# Load functions for algorithm 2
source("stage1estimator.R")
sourceCpp("stage2solver.cpp")

bic_gelnet = function(X, ebic_par = 0, algorithm = 1, alpha_seq = NULL, lambda_seq = NULL, toler = 1e-5, maxit = 500){

  S = cov(X)
  n = dim(X)[1]
  p = dim(X)[2]

  if(is.null(alpha_seq)){
    alpha_seq = seq(0,1,length.out = 11)
  }
  if(is.null(lambda_seq)){
    lambda_seq = seq(0,2,length.out = 21)
  }

  covariance.penalized = array(0,c(p,p,length(alpha_seq)*length(lambda_seq)))
  precision.penalized = array(0,c(p,p,length(alpha_seq)*length(lambda_seq)))
  AdjMat = array(0,c(p,p,length(alpha_seq)*length(lambda_seq)))
  ebic_data = matrix(0,length(alpha_seq)*length(lambda_seq),4)
  
  count = 1

  for(ii in alpha_seq){

    for(jj in lambda_seq){
      
      # Estimation of the extended BIC measure for the actual set of penalization parameters (ii, jj)
      if(algorithm==1){
        O_pp_init = solve(diag(diag(S+jj*ii*diag(p))))
        # To ensure a symmetric matrix
        O_pp_init = (O_pp_init+t(O_pp_init))/2
        W_pp_init = S + jj*ii*diag(p) + 2*jj*(1-ii)*O_pp_init
        full_out = gelnet_mio_cpp(sampleCov=S, Winit=W_pp_init, Oinit=O_pp_init, alpha=ii, lambda=jj, toler=toler, maxiters=maxit, penalizediagonal=TRUE)
      }
      if(algorithm==2){
        E0 = E0estimator_elasticnet(data = X, alpha = ii, lambda = jj, rule = "AND", toler = toler)
        W_init = S
        full_out = stage2solver(S = S, Winit = W_init, E0 = E0, toler = toler, maxit = maxit)
      }


      covariance.penalized[,,count] = full_out$W
      precision.penalized[,,count] = full_out$O
      AdjMat[,,count] = ((full_out$O!=0)*!diag(p))

      # Calculate eBIC criterion
      nedges_active = sum(AdjMat[,,count])/2
      logLik = (-n*p*log(2*pi))/2 + (n/2)*log(det(full_out$O)) - (n/2)*sum(diag(S%*%full_out$O))
      # logLik = log(det(full_out$O)) - sum(diag(S%*%full_out$O))
      ebic = -2*logLik + nedges_active*log(n) + 4*nedges_active*ebic_par*log(p)

      # Save data for the set (ii, jj)
      ebic_data[count, 1] = ii
      ebic_data[count, 2] = jj
      ebic_data[count, 3] = ebic
      ebic_data[count, 4] = logLik

      count = count + 1
    }
  }
  opt_idx = which.min(ebic_data[,3])
  colnames(ebic_data) = c("alpha","lambda","ebic","logLik")
  optimal = list(opt_idx = opt_idx, opt_ebic_data = ebic_data[opt_idx,], opt_O = precision.penalized[,,opt_idx], opt_W = covariance.penalized[,,opt_idx], opt_Adj = AdjMat[,,opt_idx])
  output = list(ebic_data = ebic_data, precision.penalized_path = precision.penalized, covariance.penalized_path = covariance.penalized, AdjMat_path = AdjMat, optimal_model = optimal)

  return(output)
}
