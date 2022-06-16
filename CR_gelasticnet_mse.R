library(glmnet)

# La suddivisione in k-fold è DIVERSA per ogni cross-validazione di ciascuna regressione conditionale
CRgelnetCV = function(X, kcv=5, lambdas=NULL, alphas=NULL, thr=1e-5, maxit=500, sym_mtd = 1){
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  cv_data = array(0, c(length(lambdas)*length(alphas),4))
  
  count = 0
  
  for(j in alphas){
    tmp_mse = array(0, c(p,length(lambdas)) )
    tmp_stdmse = array(0, c(p,length(lambdas)) )
    for(i in 1:p){
      cv_out = cv.glmnet(X[,-i], X[,i], family="gaussian", alpha=j, lambda=rev(lambdas), nfolds=kcv, type.measure = "mse", intercept=FALSE, thresh=thr, maxit=maxit)
      tmp_mse[i,] = rev(cv_out$cvm)
      tmp_stdmse[i,] = rev(cv_out$cvsd)
    }
    
    tmp_avg_mse = apply(tmp_mse, 2, mean)
    tmp_avg_stdmse = apply(tmp_stdmse, 2, mean)
    
    cv_data[( count*length(lambdas)+1 ):( (count+1)*length(lambdas) ), 1] = rep(j,length(lambdas))
    cv_data[( count*length(lambdas)+1 ):( (count+1)*length(lambdas) ), 2] = lambdas
    cv_data[( count*length(lambdas)+1 ):( (count+1)*length(lambdas) ), 3] = tmp_avg_mse
    cv_data[( count*length(lambdas)+1 ):( (count+1)*length(lambdas) ), 4] = tmp_avg_stdmse
    
    count = count+1
  }
  
  idx_min_avg_mse = which( cv_data[,3] == min(cv_data[,3]) )
  nopts = length(idx_min_avg_mse)
  idx_min = idx_min_avg_mse[1]
  alpha_opt = cv_data[idx_min,1]
  lambda_opt = cv_data[idx_min,2]
  
  O = array( 0, c(p,p) )
  S = cov(X)
  
  for(i in 1:p){
    est = glmnet(X[,-i], X[,i], family="gaussian", alpha=alpha_opt, lambda=lambda_opt, intercept=FALSE, thresh=thr, maxit=maxit)
    b_est = matrix( est$beta, c(p-1,1) )
    O[i,i] = solve( S[i,i] - 2*t(b_est)%*%S[-i,i] + t(b_est)%*%S[-i,-i]%*%b_est )
    O[-i,i] = (-1)*O[i,i]*b_est
  }
  
  # L2-symmetric
  OL2 = ( O + t(O) )/2
  
  # MinEl-symmetric
  OMinEl = array( 0, c(p,p) )
  for(i in 1:p){
    for(j in 1:p){
      if(abs(O[i,j]) <= abs(O[j,i])){
        OMinEl[i,j] = O[i,j]
        OMinEl[j,i] = O[i,j]
      }else{
        OMinEl[i,j] = O[j,i]
        OMinEl[j,i] = O[j,i]
      }
    }
  }
  
  WL2 = solve(OL2)
  WL2 = ( WL2+t(WL2) )/2
  AdjMatL2 = (OL2!=0)*(OL2!=0)
  diag(AdjMatL2) = 0
  
  WMinEl = solve(OMinEl)
  WMinEl = ( WMinEl+t(WMinEl) )/2
  AdjMatMinEl = (OMinEl!=0)*(OMinEl!=0)
  diag(AdjMatMinEl) = 0
  
  # L2-symmetric - choice
  if(sym_mtd==1){
    O = ( O + t(O) )/2
  }
  
  # MinEl-symmetric - choice
  if(sym_mtd==2){
    Oa = O
    O = array( 0, c(p,p) )
    for(i in 1:p){
      for(j in 1:p){
        if(abs(Oa[i,j]) <= abs(Oa[j,i])){
          O[i,j] = Oa[i,j]
          O[j,i] = Oa[i,j]
        }else{
          O[i,j] = Oa[j,i]
          O[j,i] = Oa[j,i]
        }
      }
    }
  }
  
  W = solve(O)
  W = ( W+t(W) )/2
  AdjMat = (O!=0)*(O!=0)
  diag(AdjMat) = 0
  
  optimal = list(opt_idx = idx_min,
                 opt_cv_data = cv_data[idx_min,],
                 opt_OL2 = OL2,
                 opt_WL2 = WL2,
                 opt_AdjL2 = AdjMatL2,
                 opt_OMinEl = OMinEl,
                 opt_WMinEl = WMinEl,
                 opt_AdjMinEl = AdjMatMinEl,
                 opt_O = O,
                 opt_W = W,
                 opt_Adj = AdjMat,
                 nopts = nopts)
  output = list(cv_data = cv_data, optimal_model = optimal)
  return(output)
}

# La suddivisione in k-fold è UGUALE per ogni cross-validazione di ciascuna regressione conditionale
# Setto manualmente per tutte in anticipo i folds da utilizzare!
CRgelnetCV_fixed = function(X, kcv=5, lambdas=NULL, alphas=NULL, thr=1e-5, maxit=500, sym_mtd = 1){
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  subsample_idxs = sample(cut(1:n, breaks=kcv, labels=FALSE), n,replace=FALSE)
  
  cv_data = array(0, c(length(lambdas)*length(alphas),4))
  
  count = 0
  
  for(j in alphas){
    tmp_mse = array(0, c(p,length(lambdas)) )
    tmp_stdmse = array(0, c(p,length(lambdas)) )
    for(i in 1:p){
      cv_out = cv.glmnet(X[,-i], X[,i], family="gaussian", alpha=j, lambda=rev(lambdas), foldid=subsample_idxs, type.measure = "mse", intercept=FALSE, thresh=thr, maxit=maxit)
      tmp_mse[i,] = rev(cv_out$cvm)
      tmp_stdmse[i,] = rev(cv_out$cvsd)
    }
    
    tmp_avg_mse = apply(tmp_mse, 2, mean)
    tmp_avg_stdmse = apply(tmp_stdmse, 2, mean)
    
    cv_data[( count*length(lambdas)+1 ):( (count+1)*length(lambdas) ), 1] = rep(j,length(lambdas))
    cv_data[( count*length(lambdas)+1 ):( (count+1)*length(lambdas) ), 2] = lambdas
    cv_data[( count*length(lambdas)+1 ):( (count+1)*length(lambdas) ), 3] = tmp_avg_mse
    cv_data[( count*length(lambdas)+1 ):( (count+1)*length(lambdas) ), 4] = tmp_avg_stdmse
    
    count = count+1
  }
  
  idx_min_avg_mse = which( cv_data[,3] == min(cv_data[,3]) )
  nopts = length(idx_min_avg_mse)
  idx_min = idx_min_avg_mse[1]
  alpha_opt = cv_data[idx_min,1]
  lambda_opt = cv_data[idx_min,2]
  
  O = array( 0, c(p,p) )
  S = cov(X)
  
  for(i in 1:p){
    est = glmnet(X[,-i], X[,i], family="gaussian", alpha=alpha_opt, lambda=lambda_opt, intercept=FALSE, thresh=thr, maxit=maxit)
    b_est = matrix( est$beta, c(p-1,1) )
    O[i,i] = solve( S[i,i] - 2*t(b_est)%*%S[-i,i] + t(b_est)%*%S[-i,-i]%*%b_est )
    O[-i,i] = (-1)*O[i,i]*b_est
  }
  
  # L2-symmetric
  OL2 = ( O + t(O) )/2
  
  # MinEl-symmetric
  OMinEl = array( 0, c(p,p) )
  for(i in 1:p){
    for(j in 1:p){
      if(abs(O[i,j]) <= abs(O[j,i])){
        OMinEl[i,j] = O[i,j]
        OMinEl[j,i] = O[i,j]
      }else{
        OMinEl[i,j] = O[j,i]
        OMinEl[j,i] = O[j,i]
      }
    }
  }
  
  WL2 = solve(OL2)
  WL2 = ( WL2+t(WL2) )/2
  AdjMatL2 = (OL2!=0)*(OL2!=0)
  diag(AdjMatL2) = 0
  
  WMinEl = solve(OMinEl)
  WMinEl = ( WMinEl+t(WMinEl) )/2
  AdjMatMinEl = (OMinEl!=0)*(OMinEl!=0)
  diag(AdjMatMinEl) = 0
  
  # L2-symmetric - choice
  if(sym_mtd==1){
    O = ( O + t(O) )/2
  }
  
  # MinEl-symmetric - choice
  if(sym_mtd==2){
    Oa = O
    O = array( 0, c(p,p) )
    for(i in 1:p){
      for(j in 1:p){
        if(abs(Oa[i,j]) <= abs(Oa[j,i])){
          O[i,j] = Oa[i,j]
          O[j,i] = Oa[i,j]
        }else{
          O[i,j] = Oa[j,i]
          O[j,i] = Oa[j,i]
        }
      }
    }
  }
  
  W = solve(O)
  W = ( W+t(W) )/2
  AdjMat = (O!=0)*(O!=0)
  diag(AdjMat) = 0
  
  optimal = list(opt_idx = idx_min,
                 opt_cv_data = cv_data[idx_min,],
                 opt_OL2 = OL2,
                 opt_WL2 = WL2,
                 opt_AdjL2 = AdjMatL2,
                 opt_OMinEl = OMinEl,
                 opt_WMinEl = WMinEl,
                 opt_AdjMinEl = AdjMatMinEl,
                 opt_O = O,
                 opt_W = W,
                 opt_Adj = AdjMat,
                 nopts = nopts)
  output = list(cv_data = cv_data, optimal_model = optimal)
  return(output)
}
