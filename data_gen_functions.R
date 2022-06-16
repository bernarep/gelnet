# Additional data generation procedures

library(mvtnorm) # Needed to genertate data from a multivariate normal, given the covariance matrix generated

# Core-Periphery data gen
coreperiphery_graph_gen = function(n, p, k, b, v, u){
  
  # This function generates random vectors from a multivariate normal distribution whose conditional independence
  # is represented by a core-periphery topology graph
  # Comments: this function follows generation procedure used by Gabriele Torri (2018)
  
  # Inputs:
  # n = number of random samples/vectors that must be generated
  # p = dimension of each random vector
  # k = number of nodes in the fully connected core
  # b = probability that a edge is created between a node in the core and a node outside the core
  # v = value of off-diagonal elements in the precision matrix that influences the partial corrlation
  # u = value added to the diagonal elements of the precision matrix that, along with v influences
  #     the partial correlation
  
  # Outputs:
  # A list with : - data: n x p matrix whose rows are observation of random vector
  #               - theta: adjaciency matrix of the graphical model underling the multivariate normal distribution
  #               - omega: precision matrix of the multivariate normal distribution
  #               - sigma: covariance matrix of the multivariate normal distribution
  
  
  theta = matrix(rep(0,p*p), c(p,p))
  omega = matrix(rep(0,p*p), c(p,p))
  diag(omega) = 1
  sigma = matrix(rep(0,p*p), c(p,p))
  data = matrix(rep(0,n*p), c(n,p))
  
  # generate core
  theta[1:k,1:k] = matrix(rep(1,k*k),c(k,k))
  theta[1:k,(k+1):p] = matrix(rbinom(k*(p-k), 1, b), c(k,p-k))
  theta[(k+1):p,1:k] = t(theta[1:k,(k+1):p])
  diag(theta) = 0
  
  tmp1 = theta*v
  tmp2 = eigen(tmp1)$values
  e = which.min(tmp2)
  e = tmp2[e]
  omega = theta*v+(abs(e)+0.1+u)*omega
  sigma = solve(omega)
  data = rmvnorm(n=n, mean=rep(0,p), sigma=sigma)
  
  return(list(theta = theta, omega = omega, sigma = sigma, data = data))
}

# Small-World data gen
smallworld_graph_gen = function(n, p, k, b, v, u){
  
  # This function generates random vectors from a multivariate normal distribution whose conditional independence
  # is represented by a small-world topology graph
  # Comments: this function follows the small-world generation procedure proposed by Watts/Strogatz (1998)
  
  # Inputs:
  # n = number of random samples/vectors that must be generated
  # p = dimension of each random vector
  # k = number of left/right neighbors (total of 2k) in the band-ring topology before rewiring
  # b = probability that a edge in the band-ring topology is rewired, that is a old link is deleted
  #     and a new link with a different node is created
  # v = value of off-diagonal elements in the precision matrix that influences the partial corrlation
  # u = value added to the diagonal elements of the precision matrix that, along with v influences
  #     the partial correlation
  
  # Outputs:
  # A list with : - data: n x p matrix whose rows are observation of random vector
  #               - theta: adjaciency matrix of the graphical model underling the multivariate normal distribution
  #               - omega: precision matrix of the multivariate normal distribution
  #               - sigma: covariance matrix of the multivariate normal distribution
  
  theta = matrix(rep(0,p*p), c(p,p))
  omega = matrix(rep(0,p*p), c(p,p))
  diag(omega) = 1
  # sigma = matrix(rep(0,p*p), c(p,p))
  # data = matrix(rep(0,n*p), c(n,p))
  
  # k-ring generation
  for(i in 1:p){
    if(i-k >= 1){
      theta[i,(i-k):(i-1)] = rep(1,k)
    }else{
      if(i == 1){
        theta[i,(p-(k-i)):p] = rep(1,k-(i-1))
      }else{
        theta[i,1:(i-1)] = rep(1,i-1)
        theta[i,(p-(k-i)):p] = rep(1,k-(i-1))
      }
    }
    if(i+k <= p){
      theta[i,(i+1):(i+k)] = rep(1,k)
    }else{
      if(i == p){
        theta[i,1:(1+(k-p+i-1))] = rep(1,k-p+i)
      }else{
        theta[i,(i+1):p] = rep(1,p-i)
        theta[i,1:(1+(k-p+i-1))] = rep(1,k-p+i)
      }
    }
  }
  
  # print(theta)
  
  # rewiring phase
  tmp0 = sample(1:p,p,replace=FALSE)
  for(j in tmp0){
    idx0 = which(theta[j,]==0)
    tmp_nj = idx0!=j
    idx0 = idx0[tmp_nj]
    idx1 = which(theta[j,]==1)
    tmp1 = rbinom(length(idx1), size=1, prob=1-b)
    
    # this IF statement is needed because if the number of new links exceed the number of free nodes, full rewiring it is not possible 
    if(length(idx0) >= sum(tmp1==0)){
      theta[j,idx1] = tmp1
      if(sum(tmp1) < 2*k){
        tmp2a = sample(idx0, (length(idx1)-sum(tmp1)))
        theta[j,tmp2a] = 1
      }
    }else{
      tmp2b = which(tmp1==0)
      tmp2b = tmp2b[1:length(idx0)]
      tmp2c = rep(1, length(idx1))
      tmp2c[tmp2b] = 0
      theta[j,idx1] = tmp2c
      theta[j,idx0] = 1
    }
    
    
    theta[,j] = t(theta[j,])
  }
  
  tmp3 = theta*v
  tmp4 = eigen(tmp3)$values
  e = which.min(tmp4)
  e = tmp4[e]
  omega = theta*v+(abs(e)+0.1+u)*omega
  sigma = solve(omega)
  data = rmvnorm(n=n, mean=rep(0,p), sigma=sigma)
  return(list(theta = theta, omega = omega, sigma = sigma, data = data))
}

# Data for low/high-correlation scenario
createBDmat = function(d, bd, nb, voff){
  M = array(0, c(d,d))
  diag(M) = 0.5
  for(b in 1:nb){
    for(i in 1:(bd-1)){
      M[i+(b-1)*bd,(i+(b-1)*bd+1):(b*bd)] = voff
    }
  }
  M = M+t(M)
  return(M)
}

M_hc = createBDmat(50,10,5,0.9)
M_lc = createBDmat(50,10,5,0.1)

data_hc = rmvnorm(n=100, mean=rep(0,50),sigma=M_hc)
data_lc = rmvnorm(n=100, mean=rep(0,50),sigma=M_lc)