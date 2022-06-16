library(glmnet)

E0estimator_elasticnet = function(data, alpha, lambda, rule = "OR", toler = 1e-5){
  
  # This function runs penalized regressions for each component j on the remaining components
  # of a multivariate distribution (generally gaussian, but not needed here). If coefficient
  # b_jk = 0 means that partial correlation between node j and k is 0 (which implies that 
  # the 2 variables are conditionally independet if the distribution is Gaussian). Thus, this
  # is represented with a missing esge in the conditional independence graph.
  # Unfortunately b_jk = 0 doesn't imply that b_kj = 0 with empirical data, so an AND or OR
  # rule is used to decide if the undirected edge jk must be included.
  
  # Inputs:
  # data = n x p matrix with n observations of p variables in the rows. !! Better if centered and normalized !!
  # alpha = parameter that balance the penalties
  #         - If alpha = 1 lasso regression
  #         - If alpha = 0 ridge regression
  # lambda = parameter that regulates the strenght of the penalty
  # rule = rule used to decide how to include and edge in the E0 estimator
  #        - OR : edge included if b_jk = 0 OR b_kj = 0
  #        - AND : edge included if b_jk = 0 AND b_kj = 0
  # toler = convergence parameter for the estimation of the penalized regressions
  
  # Outputs:
  # E0estimator : symmetric matrix with 1s in the jk (and so kj) entries when the undirected edge jk is estimated
  #               to be in E0.
  #               Note: the jj entries are set to be 1 even if this doesn't mean that there is an edge. This is
  #                     done because it is useful to the following estimation of precision matrix
  
  n = dim(data)[1]
  p = dim(data)[2]
  precision = matrix(rep(1,p*p),c(p,p))
  for(i in 1:p){
    y = data[,i]
    X = data[,-i]
    est = glmnet(X, y, alpha = alpha, lambda = lambda, thresh = toler)
    precision[i,-i] = matrix(est$beta, c(1,p))
  }
  
  precision_bin = precision != 0
  if(rule=="OR"){
    E0estimator = precision_bin+t(precision_bin)
    E0estimator = E0estimator != 0
    E0estimator = E0estimator*E0estimator
  }
  if(rule=="AND"){
    E0estimator = precision_bin*t(precision_bin)
  }
  
  # Set diagonal elements equal to 1, useful for following procedure, but it doesn't mean it has a node has a self loop!
  diag(E0estimator) = 1
  return(E0estimator)
}


